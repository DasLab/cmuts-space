"""Runs one job's pipeline as cmuts subprocess calls.

On-disk job layout::

    {job_dir}/
        {condition}.h5    final normalized output, one per condition
        {condition}.csv   the same output as CSV
        meta.json         status and condition names
        settings.json     every option the run used
        log.txt           the streaming log
        uploads/          the staged inputs (skipped by the zip download)
        work/             intermediate BAMs and rates (deleted on success)
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import traceback
from dataclasses import dataclass, field

MAX_UPLOAD_MB = int(os.environ.get("CMUTS_MAX_UPLOAD_MB", "500"))
MAX_CONDITIONS = int(os.environ.get("CMUTS_MAX_CONDITIONS", "5"))
RESULTS_TTL_HOURS = int(os.environ.get("CMUTS_RESULTS_TTL_HOURS", "72"))
STEP_TIMEOUT_SEC = int(os.environ.get("CMUTS_STEP_TIMEOUT_SEC", "600"))
THREADS = int(os.environ.get("CMUTS_THREADS", "2"))
CLEANUP_INTERVAL_SEC = int(os.environ.get("CMUTS_CLEANUP_INTERVAL_SEC", "3600"))

# The app tries these directories in order and uses the first one it can
# write to. Hugging Face mounts persistent storage at /data when the space
# has it.
RESULTS_DIR_CANDIDATES = ("/data/results", "/tmp/cmuts-space-results")


def _accepts_a_write(path: str) -> bool:
    """Tests whether the process can write to the directory. Creates the
    directory if it is missing."""
    probe = os.path.join(path, ".write-probe")
    try:
        os.makedirs(path, exist_ok=True)
        with open(probe, "w"):
            pass
        os.unlink(probe)
        return True
    except OSError:
        return False


def _first_writable(candidates: tuple[str, ...]) -> str:
    """Returns the first candidate that the process can write to. Raises
    RuntimeError if it can write to none of them."""
    for path in candidates:
        if _accepts_a_write(path):
            return path
    raise RuntimeError(
        "the app cannot write to any of these results directories: "
        + ", ".join(candidates)
    )


def _results_dir() -> str:
    """Returns the writable results directory. Uses the directory that
    CMUTS_RESULTS_DIR names, and the candidate list when that variable is not
    set. The write test lets the app run whether or not persistent storage is
    attached."""
    configured = os.environ.get("CMUTS_RESULTS_DIR")
    return _first_writable((configured,) if configured else RESULTS_DIR_CANDIDATES)


RESULTS_DIR = _results_dir()

# The settings one job used, written at submission and offered as a download.
SETTINGS_FILE = "settings.json"

# The performance arguments the server sets on each subcommand.
SERVER_ARGS = {
    "align": ["--threads", str(THREADS)],
    "hmm": ["--workers", str(THREADS)],
}


# The roles a read file can play within one condition.
CONDITION_ROLES = ("treated", "untreated", "denatured")


@dataclass
class ConditionInput:
    """One condition: a display name and the read files per role. The app
    stages the files before the pipeline reads them."""

    name: str
    treated: list[str]
    untreated: list[str] = field(default_factory=list)
    denatured: list[str] = field(default_factory=list)


@dataclass
class JobState:
    """Mutable in-memory state for an in-flight job."""

    job_id: str
    status: str = "running"  # "running" | "done" | "error"
    log_lines: list[str] = field(default_factory=list)
    error: str | None = None

    def log(self, msg: str) -> None:
        self.log_lines.append(msg)


# --- Job directories, metadata, cleanup ---


def job_dir_for(job_id: str) -> str:
    return os.path.join(RESULTS_DIR, job_id)


def write_meta(job_dir: str, meta: dict) -> None:
    with open(os.path.join(job_dir, "meta.json"), "w") as f:
        json.dump(meta, f)


def write_settings(job_dir: str, settings: dict) -> None:
    with open(os.path.join(job_dir, SETTINGS_FILE), "w") as f:
        json.dump(settings, f, indent=2)


def read_meta(job_dir: str) -> dict | None:
    """Returns the job's metadata, or None if it is missing or unreadable."""
    path = os.path.join(job_dir, "meta.json")
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError):
        return None


def write_log(job_dir: str, log_lines: list[str]) -> None:
    with open(os.path.join(job_dir, "log.txt"), "w") as f:
        f.write("\n".join(log_lines))


def cleanup_old_results() -> None:
    """Deletes result directories older than RESULTS_TTL_HOURS."""
    cutoff = time.time() - RESULTS_TTL_HOURS * 3600
    if not os.path.isdir(RESULTS_DIR):
        return
    for entry in os.scandir(RESULTS_DIR):
        if not entry.is_dir():
            continue
        meta = read_meta(entry.path)
        created = meta.get("created_at", 0) if meta else entry.stat().st_mtime
        if created < cutoff:
            shutil.rmtree(entry.path, ignore_errors=True)


def start_cleaner() -> None:
    """Starts a thread that prunes old results. The thread prunes once at
    startup and then every CLEANUP_INTERVAL_SEC seconds."""
    def loop() -> None:
        while True:
            try:
                cleanup_old_results()
            except OSError:
                pass
            time.sleep(CLEANUP_INTERVAL_SEC)

    threading.Thread(target=loop, daemon=True).start()


def safe_name(raw: str | None, fallback: str = "condition") -> str:
    """Returns a name that is safe to use in a file name and as an HDF5
    label."""
    name = re.sub(r"[^\w\-]", "_", (raw or "").strip())
    return name or fallback


# --- Subprocess steps ---


class StepFailed(RuntimeError):
    """Signals that one step of the pipeline failed. The message is written
    for the user, and names the step and the reason."""


def cmuts_errors(cmd: list[str], stderr: str) -> list[str]:
    """Returns the lines of the error output that the cmuts subcommand wrote
    itself. Each of these lines starts with the name of the subcommand."""
    prefix = f"{cmd[0]} {cmd[1]}: "
    return [line for line in stderr.splitlines() if line.startswith(prefix)]


def failure_message(step: str, cmd: list[str], result) -> str:
    """Returns the message for a step whose command exited with an error. The
    message holds the errors that cmuts reported, or the exit code if cmuts
    reported none."""
    reported = cmuts_errors(cmd, result.stderr or "")
    if not reported:
        return f"{step} failed: {cmd[0]} {cmd[1]} exited with code {result.returncode}."
    return "\n".join([f"{step} failed:", *reported])


def make_runner(state: JobState, cwd: str):
    """Returns a run(cmd, step) function that logs the command and its output.
    The function raises StepFailed if the command fails or times out. The step
    describes the command in words for the error message."""

    def run(cmd: list[str], step: str, stdout_path: str | None = None) -> None:
        state.log("$ " + " ".join(cmd))
        out = open(stdout_path, "w") if stdout_path else subprocess.PIPE
        try:
            result = subprocess.run(
                cmd, cwd=cwd, stdout=out, stderr=subprocess.PIPE,
                text=True, timeout=STEP_TIMEOUT_SEC,
            )
        except subprocess.TimeoutExpired as error:
            raise StepFailed(
                f"{step} did not finish within {STEP_TIMEOUT_SEC} seconds."
            ) from error
        finally:
            if stdout_path:
                out.close()
        if not stdout_path and result.stdout:
            state.log(result.stdout.rstrip())
        if result.stderr:
            state.log(result.stderr.rstrip())
        if result.returncode != 0:
            raise StepFailed(failure_message(step, cmd, result))

    return run


def align_reads(run, fasta: str, reads: list[str], bam: str, extra: list[str],
                label: str) -> None:
    run(["cmuts", "align", "-f", fasta, "-o", bam,
         *SERVER_ARGS["align"], *extra, *reads],
        f"Aligning {label}")


def count_mutations(run, fasta: str, bam: str, h5: str, extra: list[str],
                    label: str) -> None:
    run(["cmuts", "hmm", "-f", fasta, "-o", h5,
         *SERVER_ARGS["hmm"], *extra, bam],
        f"Counting the mutations in {label}")


def subtract_background(run, treated: str, untreated: str, h5: str, extra: list[str],
                        name: str) -> None:
    run(["cmuts", "sub", "-o", h5, *extra, treated, untreated],
        f"Subtracting the untreated background of {name}")


def divide_by_control(run, rates: str, control: str, h5: str, extra: list[str],
                      name: str) -> None:
    run(["cmuts", "div", "-o", h5, *extra, rates, control],
        f"Dividing {name} by its denatured control")


def normalize_conditions(run, inputs: list[str], outputs: list[str], extra: list[str]) -> None:
    cmd = ["cmuts", "norm", *extra]
    for output in outputs:
        cmd.extend(["-o", output])
    run([*cmd, *inputs], "Normalizing the conditions")


def write_csv(run, fasta: str, h5: str, csv_path: str, name: str) -> None:
    run(["cmuts", "csv", "-f", fasta, h5], f"Writing the CSV file of {name}",
        stdout_path=csv_path)


# --- Per-condition assembly ---


def sample_rates(run, fasta: str, reads: list[str], tag: str,
                 workdir: str, extra: dict[str, list[str]], label: str) -> str:
    """Aligns one sample's reads and counts its mutations. Returns the path of
    the rates file. The label names the reads in an error message."""
    bam = os.path.join(workdir, f"{tag}.bam")
    h5 = os.path.join(workdir, f"{tag}.h5")
    align_reads(run, fasta, reads, bam, extra["align"], label)
    count_mutations(run, fasta, bam, h5, extra["hmm"], label)
    return h5


def reads_label(role: str, condition: ConditionInput) -> str:
    """Returns the words that name one role's reads of a condition."""
    return f"the {role} reads of {condition.name}"


def condition_rates(run, fasta: str, condition: ConditionInput, tag: str,
                    workdir: str, extra: dict[str, list[str]]) -> str:
    """Returns the rates of one condition before normalization: the treated
    sample, less the untreated background, divided by the denatured control.
    Each step runs only where the condition has those reads."""
    rates = sample_rates(run, fasta, condition.treated, f"{tag}-treated",
                         workdir, extra, reads_label("treated", condition))
    if condition.untreated:
        untreated = sample_rates(run, fasta, condition.untreated,
                                 f"{tag}-untreated", workdir, extra,
                                 reads_label("untreated", condition))
        subtracted = os.path.join(workdir, f"{tag}-subtracted.h5")
        subtract_background(run, rates, untreated, subtracted, extra["sub"],
                            condition.name)
        rates = subtracted
    if condition.denatured:
        denatured = sample_rates(run, fasta, condition.denatured,
                                 f"{tag}-denatured", workdir, extra,
                                 reads_label("denatured", condition))
        divided = os.path.join(workdir, f"{tag}-divided.h5")
        divide_by_control(run, rates, denatured, divided, extra["div"],
                          condition.name)
        rates = divided
    return rates


# --- Main pipeline ---


def run_pipeline(
    job_id: str,
    job_dir: str,
    fasta_path: str,
    conditions: list[ConditionInput],
    extra: dict[str, list[str]],
    state: JobState,
) -> None:
    """Runs the whole pipeline for one job. Status and logging go through
    state. The metadata and the log are written to the job directory whether
    the run succeeds or fails."""
    workdir = os.path.join(job_dir, "work")
    os.makedirs(workdir, exist_ok=True)

    try:
        tags = unique_tags(conditions)

        state.log("=== Computing per-condition rates ===")
        rates = [
            condition_rates(make_runner(state, workdir), fasta_path,
                            condition, tag, workdir, extra)
            for condition, tag in zip(conditions, tags)
        ]

        state.log("\n=== Normalizing across conditions ===")
        finals = [os.path.join(job_dir, f"{tag}.h5") for tag in tags]
        run = make_runner(state, workdir)
        normalize_conditions(run, rates, finals, extra["norm"])

        state.log("\n=== Writing CSVs ===")
        for condition, tag, final in zip(conditions, tags, finals):
            write_csv(run, fasta_path, final, os.path.join(job_dir, f"{tag}.csv"),
                      condition.name)

        state.log(f"\nDone. Computed profiles for {len(conditions)} condition(s).")
        write_meta(job_dir, {
            "job_id": job_id,
            "conditions": [
                {"name": c.name, "tag": tag}
                for c, tag in zip(conditions, tags)
            ],
            "created_at": time.time(),
        })
        write_log(job_dir, state.log_lines)
        shutil.rmtree(workdir, ignore_errors=True)
        state.status = "done"

    except StepFailed as error:
        record_failure(state, job_dir, str(error))
    except Exception:  # noqa: BLE001
        report_server_error(job_id)
        record_failure(state, job_dir, UNEXPECTED_ERROR)


# The results page shows this message when the server fails for a reason other
# than a pipeline step. The details go only to the server log.
UNEXPECTED_ERROR = "The server had an unexpected error. Please try again."

INTERRUPTED_ERROR = "The server restarted during this job; please resubmit."


def report_server_error(job_id: str) -> None:
    """Writes the traceback of the current exception to the server log."""
    print(f"cmuts: job {job_id} failed with a server error", file=sys.stderr)
    traceback.print_exc(file=sys.stderr)


def record_failure(state: JobState, job_dir: str, message: str) -> None:
    """Marks a job as failed, and writes the message to its log, its state
    and its metadata."""
    state.log(f"Error: {message}")
    state.error = message
    state.status = "error"
    write_meta(job_dir, {
        "job_id": state.job_id,
        "status": "error",
        "error": message,
        "created_at": time.time(),
    })
    write_log(job_dir, state.log_lines)


def mark_interrupted(state: JobState) -> None:
    """Records a running job as interrupted, so its results page reports
    what happened instead of finding nothing. Called from the SIGTERM
    handler when the container is replaced."""
    if state.status != "running":
        return
    try:
        record_failure(state, job_dir_for(state.job_id), INTERRUPTED_ERROR)
    except OSError:
        pass


def unique_tags(conditions: list[ConditionInput]) -> list[str]:
    """Returns one tag per condition, safe to use in a file name. Adds a
    number where two conditions would otherwise get the same tag."""
    tags: list[str] = []
    for i, condition in enumerate(conditions):
        tag = safe_name(condition.name, f"condition_{i + 1}")
        base, k = tag, 2
        while tag in tags:
            tag = f"{base}_{k}"
            k += 1
        tags.append(tag)
    return tags
