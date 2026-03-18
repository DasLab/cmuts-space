#!/usr/bin/env python3
"""Gradio app for cmuts: chemical mutation profiling for RNA structure analysis."""

from __future__ import annotations

import glob
import json
import os
import re
import shutil
import subprocess
import tempfile
import time
import traceback
import uuid
from dataclasses import dataclass

import gradio as gr
import h5py
import numpy as np
import plotly.graph_objects as go
import cmuts
from cmuts.visualize.plotly import (
    plot_correlation,
    plot_coverage,
    plot_cumulative_reads,
    plot_examples,
    plot_heatmap,
    plot_mi,
    plot_profile,
    plot_read_hist,
    plot_snr_scaling,
    plot_termination,
)
from fastapi import FastAPI
from fastapi.responses import FileResponse, HTMLResponse, RedirectResponse


# --- Constants and paths ---

EXAMPLES_DIR = os.environ.get("CMUTS_EXAMPLES_DIR", os.path.join(os.path.dirname(__file__), "examples"))
MAX_FASTQ_MB = int(os.environ.get("CMUTS_MAX_FASTQ_MB", "500"))
RESULTS_TTL_HOURS = int(os.environ.get("CMUTS_RESULTS_TTL_HOURS", "48"))
DEFAULT_GROUP_NAME = "profile"
PIPELINE_TIMEOUT_SEC = int(os.environ.get("CMUTS_PIPELINE_TIMEOUT_SEC", "600"))

_default_results_dir = "/data/results" if os.path.isdir("/data") else "/tmp/cmuts_results"
RESULTS_DIR = os.environ.get("CMUTS_RESULTS_DIR", _default_results_dir)
os.makedirs(RESULTS_DIR, exist_ok=True)

_FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


# --- Dataclasses ---


@dataclass
class AlignConfig:
    trim_5: str = ""
    trim_3: str = ""
    local_align: bool = False


@dataclass
class CoreConfig:
    min_mapq: int = 10
    min_phred: int = 10
    min_length: int = 2
    max_length: int = 1024
    no_insertions: bool = True
    no_mismatches: bool = False
    strand: str = "both"


@dataclass
class NormConfig:
    norm_method: str = "ubr"
    no_insertions: bool = True
    no_deletions: bool = False
    clip_low: bool = False
    clip_high: bool = False
    blank_5p: int = 0
    blank_3p: int = 0
    blank_cutoff: int = 10
    norm_cutoff: int = 500
    norm_percentile: int = 90


# --- Input validation ---


def _sanitize_group_name(raw: str | None) -> str:
    """Normalize a user-supplied group name to a safe HDF5 path component."""
    name = re.sub(r"[^\w\-]", "_", (raw or "").strip())
    return name or DEFAULT_GROUP_NAME


def _fastq_stem(path: str) -> str:
    """Return the sample name from a FASTQ path by stripping known extensions.

    Handles multi-dot filenames like 'sample.rep1.fastq.gz' correctly.
    """
    name = os.path.basename(path)
    for suffix in _FASTQ_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return os.path.splitext(name)[0]


def _file_size_mb(path: str) -> float:
    return os.path.getsize(path) / (1024 * 1024)


# --- CLI command builders ---


def _build_align_cmd(
    fasta_path: str,
    output_dir: str,
    fastq_files: list[str],
    cfg: AlignConfig,
) -> list[str]:
    cmd = ["cmuts", "align", "--fasta", fasta_path, "--output", output_dir]
    if cfg.trim_5.strip():
        cmd.extend(["--trim-5", cfg.trim_5.strip()])
    if cfg.trim_3.strip():
        cmd.extend(["--trim-3", cfg.trim_3.strip()])
    if cfg.local_align:
        cmd.append("--local")
    cmd.extend(fastq_files)
    return cmd


def _build_core_cmd(
    fasta_path: str,
    output_h5: str,
    bam_files: list[str],
    cfg: CoreConfig,
) -> list[str]:
    cmd = [
        "cmuts", "core",
        "-f", fasta_path,
        "-o", output_h5,
        "--min-mapq", str(cfg.min_mapq),
        "--min-phred", str(cfg.min_phred),
        "--min-length", str(cfg.min_length),
        "--max-length", str(cfg.max_length),
    ]
    if cfg.no_insertions:
        cmd.append("--no-insertions")
    if cfg.no_mismatches:
        cmd.append("--no-mismatches")
    if cfg.strand == "forward":
        cmd.append("--no-reverse")
    elif cfg.strand == "reverse":
        cmd.append("--only-reverse")
    cmd.extend(bam_files)
    return cmd



# --- Intermediate checks ---


def _check_bam_files(alignments_dir: str) -> list[str]:
    """Return sorted relative BAM paths, or raise if none found."""
    bam_files = sorted(glob.glob(os.path.join(alignments_dir, "*.bam")))
    if not bam_files:
        raise RuntimeError(
            f"Alignment produced no BAM files in {alignments_dir}. "
            "Check the log for bowtie2 errors — the reference FASTA may not "
            "match the reads, or the FASTQ may be empty."
        )
    return bam_files


def _check_output_h5(path: str, step: str) -> None:
    """Raise if an expected HDF5 output file is missing."""
    if not os.path.isfile(path):
        raise RuntimeError(
            f"{step} did not produce output file: {os.path.basename(path)}. "
            "Check the log for errors."
        )


# --- HDF5 reading and plotting ---


def _build_plots(
    mod: cmuts.ProbingData,
    nomod: cmuts.ProbingData | None,
    combined: cmuts.ProbingData,
    name: str,
) -> dict[str, go.Figure | None]:
    """Build all plots from in-memory ProbingData objects."""
    plots: dict[str, go.Figure | None] = {}

    plots["profile"] = plot_examples(
        np.asarray(combined.reactivity), np.asarray(combined.error), name,
    )
    plots["mod_heatmap"] = plot_heatmap(np.asarray(combined.heatmap), name)
    plots["termination"] = plot_termination(np.asarray(combined.terminations), name)
    plots["coverage"] = plot_coverage(
        np.asarray(combined.coverage), np.asarray(combined.reads), name,
    )

    is_multi = not combined.single()
    reads = np.asarray(combined.reads)
    plots["read_hist"] = plot_read_hist(reads, name) if is_multi else None
    plots["cumulative_reads"] = plot_cumulative_reads(reads, name) if is_multi else None

    plots["snr_scaling"] = plot_snr_scaling(mod, nomod, combined, name)

    if combined.mi is not None:
        plots["mi"] = plot_mi(np.asarray(combined.mi)[0], name)
    else:
        plots["mi"] = None

    if combined.covariance is not None:
        plots["correlation"] = plot_correlation(np.asarray(combined.covariance)[0], name)
    else:
        plots["correlation"] = None

    return plots


def _read_profiles(h5_path: str, group_name: str) -> tuple[np.ndarray, list[str]]:
    """Read reactivity profiles and sequence names from an HDF5 file."""
    with h5py.File(h5_path, "r") as f:
        grp = f[group_name] if group_name in f else f
        reactivity = np.array(grp["reactivity"])
        sequences = None
        if "sequence" in f:
            sequences = [
                s.decode() if isinstance(s, bytes) else s for s in f["sequence"]
            ]
    names = []
    for i in range(reactivity.shape[0]):
        seq = sequences[i] if sequences and i < len(sequences) else None
        if seq and len(seq) > 50:
            names.append(seq[:50] + "...")
        elif seq:
            names.append(seq)
        else:
            names.append(f"Sequence {i + 1}")
    return reactivity, names


def _build_stats_table(h5_path: str, group_name: str) -> list[list[str]]:
    """Build a stats table as a list of [Statistic, Value] rows."""
    with h5py.File(h5_path, "r") as f:
        grp = f[group_name] if group_name in f else f
        reactivity = np.array(grp["reactivity"])
        reads = np.array(grp["reads"])
        error = np.array(grp["error"])
        snr = np.array(grp["SNR"])

    n_refs = reactivity.shape[0]
    seq_len = reactivity.shape[1]
    total_reads = int(reads.sum())
    valid = np.isfinite(reactivity)

    rows = [
        ["References", f"{n_refs:,}"],
        ["Reference length", f"{seq_len:,}"],
        ["Total reads", f"{total_reads:,}"],
        ["Mean reads per reference", f"{np.mean(reads):,.1f}"],
        ["Median reads per reference", f"{int(np.median(reads)):,}"],
    ]

    if valid.any():
        rows.extend([
            ["Mean reactivity", f"{np.mean(reactivity[valid]):.3f}"],
            ["Mean error", f"{np.mean(error[valid]):.3f}"],
            ["Mean SNR", f"{np.mean(snr):.2f}"],
            ["SNR > 1", f"{np.mean(snr > 1):.1%}"],
        ])

    dropout = float(np.mean(reads == 0))
    if dropout > 0:
        rows.append(["Dropout fraction", f"{dropout:.1%}"])

    return rows


# --- Results persistence ---


def cleanup_old_results() -> None:
    """Delete result directories older than RESULTS_TTL_HOURS."""
    cutoff = time.time() - RESULTS_TTL_HOURS * 3600
    if not os.path.isdir(RESULTS_DIR):
        return
    for entry in os.scandir(RESULTS_DIR):
        if entry.is_dir():
            meta_path = os.path.join(entry.path, "meta.json")
            try:
                if os.path.exists(meta_path):
                    with open(meta_path) as f:
                        created = json.load(f).get("created_at", 0)
                else:
                    created = entry.stat().st_mtime
                if created < cutoff:
                    shutil.rmtree(entry.path, ignore_errors=True)
            except Exception:
                pass


_PLOT_KEYS = [
    "profile", "mod_heatmap", "termination", "coverage",
    "read_hist", "cumulative_reads", "snr_scaling", "mi", "correlation",
]


def save_results(
    h5_path: str,
    group_name: str,
    plots: dict[str, go.Figure | None],
    stats_rows: list[list[str]],
    names: list[str],
    job_id: str | None = None,
) -> str:
    """Save pipeline results to persistent storage. Returns the job ID."""
    if job_id is None:
        job_id = uuid.uuid4().hex[:12]
    job_dir = os.path.join(RESULTS_DIR, job_id)
    os.makedirs(job_dir)

    shutil.copy(h5_path, os.path.join(job_dir, "profiles.h5"))

    meta = {
        "group_name": group_name,
        "created_at": time.time(),
        "names": names,
        "stats_rows": stats_rows,
    }
    with open(os.path.join(job_dir, "meta.json"), "w") as f:
        json.dump(meta, f)

    for key in _PLOT_KEYS:
        fig = plots.get(key)
        if fig is not None:
            with open(os.path.join(job_dir, f"{key}.json"), "w") as f:
                f.write(fig.to_json())

    return job_id


# --- Pipeline orchestration ---


def _empty_plot() -> go.Figure:
    fig = go.Figure()
    fig.update_layout(template="plotly_white", height=400)
    return fig




def _plot_update(fig: go.Figure | None):
    """Wrap a figure in a gr.update that shows/hides the component."""
    if fig is None:
        return gr.update(visible=False, value=None)
    return gr.update(visible=True, value=fig)


def _progress_yield(result_url: str, log_lines: list[str]) -> tuple:
    """Build the in-progress yield tuple. Single source of truth for the
    yield shape: (file, profile, dropdown, stats, url,
    mod_heatmap, termination, coverage, read_hist, cumulative_reads,
    snr_scaling, mi, correlation, log)."""
    hidden = gr.update(visible=False, value=None)
    return (hidden, hidden, None, hidden, hidden) + (hidden,) * (len(_PLOT_KEYS) - 1) + ("\n".join(log_lines),)


def run_pipeline(
    fasta_file: str,
    mod_fastq: str,
    nomod_fastq: str | None,
    group_name: str,
    align_cfg: AlignConfig,
    core_cfg: CoreConfig,
    norm_cfg: NormConfig,
):
    """Run the full cmuts pipeline: align -> core -> normalize.

    Yields (output_file, profile, dropdown, stats, url,
            mod_heatmap, termination, coverage, read_hist, cumulative_reads,
            snr_scaling, mi, correlation, log)
    so the log updates in real time and the interactive plots appear at the end.
    """
    if fasta_file is None or mod_fastq is None:
        raise gr.Error("Please upload a FASTA file and at least one modified FASTQ file.")

    for path, label in [(mod_fastq, "Modified FASTQ"), (nomod_fastq, "Control FASTQ")]:
        if path is not None and _file_size_mb(path) > MAX_FASTQ_MB:
            raise gr.Error(
                f"{label} is {_file_size_mb(path):.0f} MB. "
                f"The free tier has limited RAM (16 GB); files over {MAX_FASTQ_MB} MB "
                f"may cause out-of-memory errors. Consider downsampling first."
            )

    cleanup_old_results()

    workdir = tempfile.mkdtemp(prefix="cmuts_")
    outdir = os.path.join(workdir, "outputs")
    os.makedirs(outdir)

    group_name = _sanitize_group_name(group_name)

    # Generate job ID and result URL upfront so the link appears immediately
    job_id = uuid.uuid4().hex[:12]
    space_host = os.environ.get("SPACE_HOST", "")
    base = f"https://{space_host}" if space_host else ""
    result_url = f"{base}/results/{job_id}"

    log_lines: list[str] = []

    def log(msg: str) -> None:
        log_lines.append(msg)

    def run(cmd: list[str], cwd: str = outdir) -> bool:
        log(f"$ {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            cwd=cwd,
            capture_output=True,
            text=True,
            timeout=PIPELINE_TIMEOUT_SEC,
        )
        if result.stdout:
            log(result.stdout.rstrip())
        if result.stderr:
            log(result.stderr.rstrip())
        if result.returncode != 0:
            log(f"Command failed with exit code {result.returncode}")
            return False
        return True

    try:
        fasta_path = os.path.join(workdir, "ref.fasta")
        shutil.copy(fasta_file, fasta_path)

        fastq_dir = os.path.join(workdir, "fastq")
        os.makedirs(fastq_dir)

        mod_basename = os.path.basename(mod_fastq)
        mod_name = _fastq_stem(mod_fastq)
        shutil.copy(mod_fastq, os.path.join(fastq_dir, mod_basename))

        nomod_name = None
        if nomod_fastq is not None:
            nomod_basename = os.path.basename(nomod_fastq)
            nomod_name = _fastq_stem(nomod_fastq)
            shutil.copy(nomod_fastq, os.path.join(fastq_dir, nomod_basename))

        # Step 1: Align
        log("=== Step 1: Aligning reads ===")
        yield _progress_yield(result_url, log_lines)

        fastq_files = sorted(glob.glob(os.path.join(fastq_dir, "*")))
        alignments_dir = os.path.join(outdir, "alignments")
        align_cmd = _build_align_cmd(fasta_path, alignments_dir, fastq_files, align_cfg)
        if not run(align_cmd, cwd=outdir):
            yield _progress_yield(result_url, log_lines)
            return
        yield _progress_yield(result_url, log_lines)

        # Step 2: Count mutations
        log("\n=== Step 2: Counting mutations ===")
        yield _progress_yield(result_url, log_lines)

        bam_abs = _check_bam_files(alignments_dir)
        bam_files = sorted(os.path.relpath(p, outdir) for p in bam_abs)
        counts_h5 = "counts.h5"
        core_cmd = _build_core_cmd(fasta_path, counts_h5, bam_files, core_cfg)
        if not run(core_cmd, cwd=outdir):
            yield _progress_yield(result_url, log_lines)
            return
        _check_output_h5(os.path.join(outdir, counts_h5), "cmuts core")
        yield _progress_yield(result_url, log_lines)

        # Step 3: Normalize (in-process via cmuts Python API)
        log("\n=== Step 3: Normalizing reactivities ===")
        yield _progress_yield(result_url, log_lines)

        mod_group = f"alignments/{mod_name}"
        nomod_group = f"alignments/{nomod_name}" if nomod_name else None

        opts = cmuts.Opts(
            cmuts.DataGroups([mod_group]),
            cmuts.DataGroups([nomod_group] if nomod_group else None),
            norm_cfg.blank_cutoff,
            not norm_cfg.no_insertions,
            not norm_cfg.no_deletions,
            norm_cfg.norm_method,
            (norm_cfg.blank_5p, norm_cfg.blank_3p),
            (norm_cfg.clip_low, norm_cfg.clip_high),
            0.05,  # significance level
        )

        counts_path = os.path.join(outdir, counts_h5)
        with h5py.File(counts_path, "r") as f:
            mod_data, nomod_data, combined = cmuts.compute_reactivity(
                f, fasta_path, opts,
            )

        # Save HDF5 output
        final_name = f"{group_name}-profiles.h5"
        final_path = os.path.join(outdir, final_name)
        combined.save(group_name, final_path)
        log("Normalization complete.")

        # Build all plots directly from in-memory ProbingData
        plots = _build_plots(mod_data, nomod_data, combined, group_name)

        reactivity, names = _read_profiles(final_path, group_name)
        dropdown_update = gr.Dropdown(
            choices=names, value=names[0],
            visible=(len(names) > 1),
        )
        stats_md = _build_stats_table(final_path, group_name)

        save_results(final_path, group_name, plots, stats_md, names, job_id=job_id)

        log(f"\nDone. Generated {len(names)} profile(s).")
        log(f"Results available at: {result_url} (expires in {RESULTS_TTL_HOURS}h)")
        yield (
            gr.update(visible=True, value=final_path),
            _plot_update(plots["profile"]), dropdown_update,
            gr.update(visible=True, value=stats_md),
            gr.update(visible=True, value=result_url),
            *[_plot_update(plots[k]) for k in _PLOT_KEYS[1:]],
            "\n".join(log_lines),
        )

    except subprocess.TimeoutExpired:
        log(f"Pipeline timed out ({PIPELINE_TIMEOUT_SEC // 60} minute limit).")
        yield _progress_yield(result_url, log_lines)
    except Exception as e:
        log(f"Error: {e}")
        log(traceback.format_exc())
        yield _progress_yield(result_url, log_lines)


# --- Gradio callbacks ---


def _run_pipeline_gradio(
    fasta_file: str,
    mod_fastq: str,
    nomod_fastq: str | None,
    group_name: str,
    norm_method: str,
    no_insertions: bool,
    no_deletions: bool,
    clip_low: bool,
    clip_high: bool,
    trim_5: str,
    trim_3: str,
    local_align: bool,
    min_mapq: int,
    min_phred: int,
    min_length: int,
    max_length: int,
    no_mismatches: bool,
    strand: str,
    blank_5p: int,
    blank_3p: int,
    blank_cutoff: int,
    norm_cutoff: int,
    norm_percentile: int,
):
    """Gradio-facing wrapper: packs flat args into dataclasses."""
    yield from run_pipeline(
        fasta_file=fasta_file,
        mod_fastq=mod_fastq,
        nomod_fastq=nomod_fastq,
        group_name=group_name,
        align_cfg=AlignConfig(
            trim_5=trim_5 or "",
            trim_3=trim_3 or "",
            local_align=local_align,
        ),
        core_cfg=CoreConfig(
            min_mapq=int(min_mapq or 10),
            min_phred=int(min_phred or 10),
            min_length=int(min_length or 2),
            max_length=int(max_length or 1024),
            no_insertions=no_insertions,
            no_mismatches=no_mismatches,
            strand=strand or "both",
        ),
        norm_cfg=NormConfig(
            norm_method=norm_method or "ubr",
            no_insertions=no_insertions,
            no_deletions=no_deletions,
            clip_low=clip_low,
            clip_high=clip_high,
            blank_5p=int(blank_5p or 0),
            blank_3p=int(blank_3p or 0),
            blank_cutoff=int(blank_cutoff or 10),
            norm_cutoff=int(norm_cutoff or 500),
            norm_percentile=int(norm_percentile or 90),
        ),
    )


def select_profile(
    seq_name: str,
    output_file: str,
    group_name: str,
) -> go.Figure:
    """Switch the displayed profile when the user picks a different sequence."""
    from cmuts.visualize.plotly import plot_profile

    if not output_file or not seq_name:
        return go.Figure()
    group_name = _sanitize_group_name(group_name)
    reactivity, names = _read_profiles(output_file, group_name)
    try:
        idx = names.index(seq_name)
    except ValueError:
        idx = 0
    with h5py.File(output_file, "r") as f:
        grp = f[group_name] if group_name in f else f
        error = np.array(grp["error"])
    return plot_profile(reactivity[idx], error[idx], names[idx])


def load_example():
    """Load bundled example files by convention from the examples directory.

    Expected layout:
        examples/
            ref.fasta (or .fa)          — required
            treated.fastq.gz (or .fq*)  — required (modified condition)
            untreated.fastq.gz          — optional (control condition)
            group.txt                   — optional (single line: group name)
    """
    fasta = None
    treated = None
    untreated = None

    for f in os.listdir(EXAMPLES_DIR):
        path = os.path.join(EXAMPLES_DIR, f)
        lower = f.lower()
        if lower.endswith((".fasta", ".fa")):
            fasta = path
        elif "untreated" in lower or "nomod" in lower or "control" in lower:
            untreated = path
        elif lower.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
            treated = path

    group_file = os.path.join(EXAMPLES_DIR, "group.txt")
    if os.path.isfile(group_file):
        with open(group_file) as f:
            group_name = f.read().strip() or "example"
    else:
        group_name = "example"

    return fasta, treated, untreated, group_name


def load_saved_result(job_id: str):
    """Load a previously saved result by job ID."""
    job_id = (job_id or "").strip()
    hidden = gr.update(visible=False, value=None)
    n_extra = len(_PLOT_KEYS) - 1  # all plot keys except "profile"
    if not job_id:
        return (hidden, None, hidden, "") + (hidden,) * n_extra + ("",)

    job_dir = os.path.join(RESULTS_DIR, job_id)
    meta_path = os.path.join(job_dir, "meta.json")

    if not os.path.isdir(job_dir):
        return (hidden, None, hidden, (
            f"Result not found. It may have expired "
            f"(results are kept for {RESULTS_TTL_HOURS} hours)."
        )) + (hidden,) * n_extra + ("",)

    with open(meta_path) as f:
        meta = json.load(f)

    loaded: dict[str, go.Figure | None] = {}
    for key in _PLOT_KEYS:
        path = os.path.join(job_dir, f"{key}.json")
        if os.path.isfile(path):
            with open(path) as f:
                loaded[key] = go.Figure(json.load(f))
        else:
            loaded[key] = None

    names = meta.get("names", [])
    dropdown_update = gr.Dropdown(choices=names, value=names[0] if names else None, visible=len(names) > 1)

    stats = meta.get("stats_rows", meta.get("stats_md", []))
    return (
        _plot_update(loaded.get("profile")), dropdown_update,
        gr.update(visible=True, value=stats), "",
        *[_plot_update(loaded.get(k)) for k in _PLOT_KEYS[1:]],
        "",
    )




# --- Gradio UI ---

with gr.Blocks(title="cmuts — RNA Chemical Probing Analysis") as demo:
    gr.Markdown(
        """
        # cmuts — RNA Chemical Probing Analysis

        Upload a FASTA reference and FASTQ file(s) from a MaP-seq experiment
        to compute normalized reactivity profiles.

        **Pipeline:** `cmuts align` &rarr; `cmuts core` &rarr; `cmuts normalize`
        &ensp;|&ensp;
        [GitHub](https://github.com/hmblair/cmuts)
        &ensp;|&ensp;
        [Documentation](https://hmblair.github.io/cmuts)
        &ensp;|&ensp;
        Free and open-source under the
        [MIT License](https://github.com/hmblair/cmuts/blob/main/LICENSE)
        """
    )

    with gr.Tab("Run"):
        gr.Markdown("### Input data")
        with gr.Row():
            with gr.Column():
                fasta_input = gr.File(label="Reference FASTA", file_types=[".fasta", ".fa"])
            with gr.Column():
                mod_input = gr.File(label="Modified FASTQ (required)", file_types=[".fastq", ".fq", ".gz"])
            with gr.Column():
                nomod_input = gr.File(label="Control FASTQ (optional)", file_types=[".fastq", ".fq", ".gz"])
        with gr.Row():
            group_name = gr.Textbox(label="Group name", value="experiment", placeholder="e.g. DMS, 2A3", scale=3)
            example_btn = gr.Button("Load example data", variant="secondary", size="sm", scale=1)

        gr.Markdown("### Options")
        with gr.Accordion("Alignment", open=False):
            with gr.Row():
                trim_5 = gr.Textbox(label="5' adapter to trim", placeholder="e.g. AGATCGGAAGAG")
                trim_3 = gr.Textbox(label="3' adapter to trim", placeholder="e.g. AGATCGGAAGAG")
                local_align = gr.Checkbox(label="Local alignment", value=False)

        with gr.Accordion("Read filtering", open=False):
            with gr.Row():
                min_mapq = gr.Slider(minimum=0, maximum=60, step=1, value=10, label="Min mapping quality")
                min_phred = gr.Slider(minimum=0, maximum=40, step=1, value=10, label="Min PHRED score")
            with gr.Row():
                min_length = gr.Number(value=2, label="Min read length", precision=0)
                max_length = gr.Number(value=1024, label="Max read length", precision=0)
            with gr.Row():
                no_mismatches = gr.Checkbox(label="Exclude mismatches", value=False)
                strand = gr.Radio(choices=["both", "forward", "reverse"], value="both", label="Strand")

        with gr.Accordion("Normalization", open=False):
            norm_method = gr.Radio(
                choices=["ubr", "outlier", "raw"],
                value="ubr",
                label="Normalization method",
            )
            with gr.Row():
                no_insertions = gr.Checkbox(label="Exclude insertions", value=True)
                no_deletions = gr.Checkbox(label="Exclude deletions", value=False)
            with gr.Row():
                clip_low = gr.Checkbox(label="Clip negative reactivities", value=False)
                clip_high = gr.Checkbox(label="Clip reactivities above 1", value=False)
            with gr.Row():
                blank_5p = gr.Number(value=0, label="Blank 5' bases", precision=0)
                blank_3p = gr.Number(value=0, label="Blank 3' bases", precision=0)
                blank_cutoff = gr.Number(value=10, label="Min reads for position", precision=0)
            with gr.Row():
                norm_cutoff = gr.Number(value=500, label="Min reads for normalization", precision=0)
                norm_percentile = gr.Slider(minimum=50, maximum=100, step=1, value=90, label="Normalization percentile")

        run_btn = gr.Button("Run Pipeline", variant="primary")

        gr.Markdown("### Results")
        output_file = gr.File(label="Output HDF5", visible=False)
        result_url = gr.Textbox(
            label=f"Result link (bookmark this — expires in {RESULTS_TTL_HOURS}h)",
            interactive=False,
            visible=False,
        )
        seq_dropdown = gr.Dropdown(label="Sequence", visible=False, interactive=True)
        output_plot = gr.Plot(label="Reactivity Profile", visible=False)
        with gr.Row():
            with gr.Column(scale=1):
                mod_heatmap_plot = gr.Plot(label="Modification Heatmap", visible=False)
            with gr.Column(scale=1):
                pass
        with gr.Row():
            termination_plot = gr.Plot(label="Termination by Position", visible=False)
            coverage_plot = gr.Plot(label="Coverage by Position", visible=False)
        with gr.Row():
            read_hist_plot = gr.Plot(label="Read Depth Distribution", visible=False)
            cumulative_reads_plot = gr.Plot(label="Cumulative Reads", visible=False)
        snr_scaling_plot = gr.Plot(label="SNR vs Read Depth", visible=False)
        with gr.Row():
            mi_plot = gr.Plot(label="Mutual Information", visible=False)
            correlation_plot = gr.Plot(label="Correlation", visible=False)
        output_stats = gr.Dataframe(
            label="Summary Statistics",
            headers=["Statistic", "Value"],
            interactive=False,
            visible=False,
        )
        with gr.Accordion("Log", open=False):
            output_log = gr.Textbox(label="Log", lines=15, max_lines=30, show_label=False)

        with gr.Accordion("Load previous results", open=False):
            with gr.Row():
                prev_job_id = gr.Textbox(label="Job ID", placeholder="e.g. a3f2b1c4d5e6", scale=3)
                load_btn = gr.Button("Load", variant="secondary", scale=1)
            load_status = gr.Textbox(label="Status", interactive=False)

        example_btn.click(
            fn=load_example,
            outputs=[fasta_input, mod_input, nomod_input, group_name],
        )

        run_btn.click(
            fn=_run_pipeline_gradio,
            inputs=[
                fasta_input,
                mod_input,
                nomod_input,
                group_name,
                norm_method,
                no_insertions,
                no_deletions,
                clip_low,
                clip_high,
                # Alignment options
                trim_5,
                trim_3,
                local_align,
                # Read filtering
                min_mapq,
                min_phred,
                min_length,
                max_length,
                no_mismatches,
                strand,
                # Normalization
                blank_5p,
                blank_3p,
                blank_cutoff,
                norm_cutoff,
                norm_percentile,
            ],
            outputs=[output_file, output_plot, seq_dropdown, output_stats, result_url,
                     mod_heatmap_plot, termination_plot, coverage_plot,
                     read_hist_plot, cumulative_reads_plot, snr_scaling_plot,
                     mi_plot, correlation_plot, output_log],
        )

        seq_dropdown.change(
            fn=select_profile,
            inputs=[seq_dropdown, output_file, group_name],
            outputs=[output_plot],
        )

        load_btn.click(
            fn=load_saved_result,
            inputs=[prev_job_id],
            outputs=[output_plot, seq_dropdown, output_stats, load_status,
                     mod_heatmap_plot, termination_plot, coverage_plot,
                     read_hist_plot, cumulative_reads_plot, snr_scaling_plot,
                     mi_plot, correlation_plot, output_log],
        )

    _load_outputs = [
        output_plot, seq_dropdown, output_stats, load_status,
        mod_heatmap_plot, termination_plot, coverage_plot,
        read_hist_plot, cumulative_reads_plot, snr_scaling_plot,
        mi_plot, correlation_plot, output_log,
    ]

    def _load_from_query(request: gr.Request):
        """Auto-load results when ?job_id= is present in the URL."""
        job_id = (request.query_params.get("job_id") or "").strip()
        if not job_id:
            return [gr.update()] * len(_load_outputs)
        return load_saved_result(job_id)

    demo.load(
        fn=_load_from_query,
        outputs=_load_outputs,
    )

    with gr.Tab("About"):
        gr.Markdown(
            """
            ## Why cmuts?

            Existing tools for analyzing MaP-seq data — such as ShapeMapper2 and
            RNAframework — were designed for single-RNA experiments and do not
            scale to modern high-throughput libraries with thousands or millions
            of reference sequences.

            **cmuts** is a ground-up rewrite in C/C++ that addresses these
            limitations:

            - **100–200x faster** than ShapeMapper2 and RNAframework. A dataset
              of 100 billion aligned reads across 24 million references was
              processed in under 24 hours on 32 cores — a task that would take
              RNAframework approximately 3 months.
            - **Constant memory footprint** regardless of library size, thanks to
              streamed single-pass I/O. Competing tools either grow linearly in
              memory or require processing one reference at a time.
            - **More accurate deletion handling.** cmuts uses a depth-first
              search algorithm to enumerate all possible positions of ambiguous
              deletions and weights them probabilistically using observed mutation
              rates. Prior tools arbitrarily assign deletions to the 3'-most
              position, which can misplace reactivity signals — particularly in
              homopolymer regions and structurally important motifs like
              kink-turns.
            - **HDF5 output** for compact storage, fast random access, and direct
              compatibility with Python and machine-learning pipelines.
            """
        )

    with gr.Tab("Help"):
        gr.Markdown(
            """
            ## Quick Start

            1. Click **Load example data** on the Run tab to populate the inputs
               with a bundled dataset.
            2. Leave the default settings and click **Run Pipeline**.
            3. The pipeline runs three steps — alignment, mutation counting, and
               normalization — and streams its progress to the log.
            4. When finished, download the output HDF5 file and explore the
               interactive reactivity profile. If multiple reference sequences
               are present, use the dropdown to switch between them.

            ## Inputs

            | Field | Description |
            |-------|-------------|
            | **Reference FASTA** | One or more RNA sequences in FASTA format. Each sequence is treated as a separate reference for alignment. |
            | **Modified FASTQ** | Reads from the chemically treated condition (e.g., DMS, 2A3, SHAPE). Compressed `.fastq.gz` is accepted. |
            | **Control FASTQ** | *(Optional)* Reads from the untreated/DMSO condition. Providing a control enables background subtraction during normalization. |
            | **Group name** | A label for the experiment, used to name the output file and the HDF5 group. |

            ## Settings

            | Setting | Description |
            |---------|-------------|
            | **Normalization method** | `ubr` (default): divides by the 90th percentile of reactivities at positions with >100 reads, giving a robust upper-bound reference. `outlier`: from the top 10% of reactivities, discards the top 2% as outliers and divides by the mean of the remaining 8%. `raw`: no normalization — returns raw mutation rates. |
            | **Exclude insertions** | Do not count inserted bases as mutations (recommended for most protocols). |
            | **Exclude deletions** | Do not count deleted bases as mutations. |
            | **Clip negative reactivities** | Set negative normalized reactivities to zero. |
            | **Clip reactivities above 1** | Cap normalized reactivities at 1.0. |

            ## Interpreting the Output

            ### Reactivity profile

            The interactive bar chart shows per-nucleotide reactivity values.
            Hover over any bar to see the exact position, nucleotide identity,
            and reactivity value. Peaks correspond to unpaired or flexible
            nucleotides; low/near-zero regions correspond to base-paired or
            otherwise protected positions.

            ### HDF5 file

            The output file (`<group>-profiles.h5`) contains normalized
            per-nucleotide reactivity profiles. Open it in Python with:

            ```python
            import h5py

            with h5py.File("experiment-profiles.h5") as f:
                for name in f:
                    reactivities = f[name][:]
                    print(name, reactivities.shape)
            ```

            Each dataset is a 1-D array of floats, one value per nucleotide.
            Higher values indicate more flexible (unpaired) positions; lower
            values indicate structured (paired) regions.

            ### Result link

            After the pipeline completes, a **result link** is displayed that
            you can bookmark or share. Opening the link loads the results
            directly into the app with the interactive plots, summary
            statistics, and HDF5 download. Results are stored for
            **{RESULTS_TTL_HOURS} hours** and automatically deleted after that.

            To reload previous results manually, expand
            **Load previous results** on the Run tab and enter the job ID.

            ## Limits

            This server has limited resources (16 GB RAM, CPU only). Very
            large FASTQ files may cause out-of-memory errors. For larger
            datasets, install cmuts locally:

            ```bash
            pip install cmuts
            ```

            See the [full documentation](https://hmblair.github.io/cmuts) for
            CLI usage and advanced options.

            ## Privacy

            All uploaded data is processed in ephemeral temporary directories.
            Pipeline results (reactivity profiles, plots, and summary
            statistics) are stored for **48 hours** to provide bookmarkable
            result links, then automatically deleted. No user accounts,
            tracking, or cookies are used.
            """
        )


# --- FastAPI routes ---

app = FastAPI()


@app.get("/results/{job_id}")
async def results_page(job_id: str):
    """Redirect to the Gradio app with the job ID as a query parameter."""
    return RedirectResponse(url=f"/?job_id={job_id}")


@app.get("/results/{job_id}/download")
async def results_download(job_id: str):
    h5_path = os.path.join(RESULTS_DIR, job_id, "profiles.h5")
    if not os.path.isfile(h5_path):
        return HTMLResponse(
            "<h1>File not found</h1><p>This result may have expired.</p>",
            status_code=404,
        )

    with open(os.path.join(RESULTS_DIR, job_id, "meta.json")) as f:
        meta = json.load(f)
    group_name = meta.get("group_name", "profiles")

    return FileResponse(
        h5_path,
        media_type="application/x-hdf5",
        filename=f"{group_name}-profiles.h5",
    )


# Mount Gradio onto the FastAPI app
app = gr.mount_gradio_app(app, demo, path="")

# Run cleanup on startup
cleanup_old_results()

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=7860)
