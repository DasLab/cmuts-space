"""Runs one job's pipeline as cmuts subprocess calls.

Framework-agnostic. Takes staged input paths plus the per-subcommand extra
arguments the option layer built, and writes everything for one job into a
single directory under RESULTS_DIR.

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
import time
import traceback
from dataclasses import dataclass, field

MAX_UPLOAD_MB = int(os.environ.get("CMUTS_MAX_UPLOAD_MB", "500"))
MAX_CONDITIONS = int(os.environ.get("CMUTS_MAX_CONDITIONS", "5"))
RESULTS_TTL_HOURS = int(os.environ.get("CMUTS_RESULTS_TTL_HOURS", "72"))
STEP_TIMEOUT_SEC = int(os.environ.get("CMUTS_STEP_TIMEOUT_SEC", "600"))
THREADS = int(os.environ.get("CMUTS_THREADS", "2"))

_default_results_dir = (
    "/data/results" if os.path.isdir("/data") else "/tmp/cmuts-space-results"
)
RESULTS_DIR = os.environ.get("CMUTS_RESULTS_DIR", _default_results_dir)

# The settings one job used, written at submission and offered as a download.
SETTINGS_FILE = "settings.json"

# The performance arguments the server sets on each subcommand.
SERVER_ARGS = {
    "align": ["--threads", str(THREADS)],
    "hmm": ["--workers", str(THREADS)],
}


@dataclass
class ConditionInput:
    """One condition: a display name and the staged read files per role."""

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
    """The job's metadata, or None where it is missing or unreadable."""
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


def safe_name(raw: str | None, fallback: str = "condition") -> str:
    """A filesystem-safe name for output files and HDF5 labels."""
    name = re.sub(r"[^\w\-]", "_", (raw or "").strip())
    return name or fallback


# --- Subprocess steps ---


def make_runner(state: JobState, cwd: str):
    """A run(cmd) closure that logs the command and its output, and raises
    on failure or timeout."""

    def run(cmd: list[str], stdout_path: str | None = None) -> None:
        state.log("$ " + " ".join(cmd))
        out = open(stdout_path, "w") if stdout_path else subprocess.PIPE
        try:
            result = subprocess.run(
                cmd, cwd=cwd, stdout=out, stderr=subprocess.PIPE,
                text=True, timeout=STEP_TIMEOUT_SEC,
            )
        except subprocess.TimeoutExpired as error:
            raise RuntimeError(
                f"{cmd[0]} {cmd[1]} timed out after {STEP_TIMEOUT_SEC}s"
            ) from error
        finally:
            if stdout_path:
                out.close()
        if not stdout_path and result.stdout:
            state.log(result.stdout.rstrip())
        if result.stderr:
            state.log(result.stderr.rstrip())
        if result.returncode != 0:
            raise RuntimeError(
                f"{cmd[0]} {cmd[1]} failed with exit code {result.returncode}"
            )

    return run


def align_reads(run, fasta: str, reads: list[str], bam: str, extra: list[str]) -> None:
    run(["cmuts", "align", "-f", fasta, "-o", bam,
         *SERVER_ARGS["align"], *extra, *reads])


def count_mutations(run, fasta: str, bam: str, h5: str, extra: list[str]) -> None:
    run(["cmuts", "hmm", "-f", fasta, "-o", h5,
         *SERVER_ARGS["hmm"], *extra, bam])


def subtract_background(run, treated: str, untreated: str, h5: str, extra: list[str]) -> None:
    run(["cmuts", "sub", "-o", h5, *extra, treated, untreated])


def divide_by_control(run, rates: str, control: str, h5: str, extra: list[str]) -> None:
    run(["cmuts", "div", "-o", h5, *extra, rates, control])


def normalize_conditions(run, inputs: list[str], outputs: list[str], extra: list[str]) -> None:
    cmd = ["cmuts", "norm", *extra]
    for output in outputs:
        cmd.extend(["-o", output])
    run([*cmd, *inputs])


def write_csv(run, fasta: str, h5: str, csv_path: str) -> None:
    run(["cmuts", "csv", "-f", fasta, h5], stdout_path=csv_path)


# --- Per-condition assembly ---


def sample_rates(run, fasta: str, reads: list[str], tag: str,
                 workdir: str, extra: dict[str, list[str]]) -> str:
    """Aligns one sample's reads and counts its mutations; returns the rates
    HDF5."""
    bam = os.path.join(workdir, f"{tag}.bam")
    h5 = os.path.join(workdir, f"{tag}.h5")
    align_reads(run, fasta, reads, bam, extra["align"])
    count_mutations(run, fasta, bam, h5, extra["hmm"])
    return h5


def condition_rates(run, fasta: str, condition: ConditionInput, tag: str,
                    workdir: str, extra: dict[str, list[str]]) -> str:
    """One condition's pre-normalization rates: treated, less the untreated
    background, over the denatured control, as far as the inputs go."""
    rates = sample_rates(run, fasta, condition.treated, f"{tag}-treated",
                         workdir, extra)
    if condition.untreated:
        untreated = sample_rates(run, fasta, condition.untreated,
                                 f"{tag}-untreated", workdir, extra)
        subtracted = os.path.join(workdir, f"{tag}-subtracted.h5")
        subtract_background(run, rates, untreated, subtracted, extra["sub"])
        rates = subtracted
    if condition.denatured:
        denatured = sample_rates(run, fasta, condition.denatured,
                                 f"{tag}-denatured", workdir, extra)
        divided = os.path.join(workdir, f"{tag}-divided.h5")
        divide_by_control(run, rates, denatured, divided, extra["div"])
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
    """Runs the full pipeline for one job. All status and logging goes
    through state; on-disk meta and log are written at the end either way."""
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
        for tag, final in zip(tags, finals):
            write_csv(run, fasta_path, final, os.path.join(job_dir, f"{tag}.csv"))

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

    except Exception as error:  # noqa: BLE001
        msg = str(error) or "Unexpected error. See log for details."
        state.log(f"Error: {msg}")
        state.log(traceback.format_exc())
        state.error = msg
        state.status = "error"
        write_meta(job_dir, {
            "job_id": job_id,
            "status": "error",
            "error": msg,
            "created_at": time.time(),
        })
        write_log(job_dir, state.log_lines)


INTERRUPTED_ERROR = "The server restarted during this job; please resubmit."


def mark_interrupted(state: JobState) -> None:
    """Records a running job as interrupted, so its results page reports
    what happened instead of finding nothing. Called from the SIGTERM
    handler when the container is replaced."""
    if state.status != "running":
        return
    state.log(f"Error: {INTERRUPTED_ERROR}")
    state.error = INTERRUPTED_ERROR
    state.status = "error"
    job_dir = job_dir_for(state.job_id)
    try:
        write_meta(job_dir, {
            "job_id": state.job_id,
            "status": "error",
            "error": INTERRUPTED_ERROR,
            "created_at": time.time(),
        })
        write_log(job_dir, state.log_lines)
    except OSError:
        pass


def unique_tags(conditions: list[ConditionInput]) -> list[str]:
    """A filesystem-safe tag per condition, disambiguated on clashes."""
    tags: list[str] = []
    for i, condition in enumerate(conditions):
        tag = safe_name(condition.name, f"condition_{i + 1}")
        base, k = tag, 2
        while tag in tags:
            tag = f"{base}_{k}"
            k += 1
        tags.append(tag)
    return tags
