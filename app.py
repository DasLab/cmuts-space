#!/usr/bin/env python3
"""Gradio app for cmuts: chemical mutation profiling for RNA structure analysis."""

from __future__ import annotations

import csv
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
    plot_pairwise_coverage,
    plot_profile,
    plot_profiles,
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
MAX_GROUPS = 5

_default_results_dir = "/data/results" if os.path.isdir("/data") else "/tmp/cmuts_results"
RESULTS_DIR = os.environ.get("CMUTS_RESULTS_DIR", _default_results_dir)
os.makedirs(RESULTS_DIR, exist_ok=True)

_FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


# --- Dataclasses ---


@dataclass
class GroupInput:
    """One experiment group with modified and optional control FASTQ files."""
    name: str
    mod_fastq: str
    nomod_fastq: str | None = None


@dataclass
class ResultUpdate:
    """Collected Gradio output updates for the results section."""
    result_url: object = None
    output_file: object = None
    csv_file: object = None
    profile_plot: object = None
    seq_dropdown: object = None
    stats: object = None
    load_status: str = ""
    mod_heatmap: object = None
    termination: object = None
    coverage: object = None
    read_hist: object = None
    cumulative_reads: object = None
    snr_scaling: object = None
    mi: object = None
    correlation: object = None
    pairwise_coverage: object = None
    log: str = ""

    @classmethod
    def hidden(cls) -> ResultUpdate:
        h = gr.update(visible=False, value=None)
        return cls(
            result_url=h, output_file=h, csv_file=h, profile_plot=h,
            seq_dropdown=None, stats=h, load_status="",
            mod_heatmap=h, termination=h, coverage=h,
            read_hist=h, cumulative_reads=h, snr_scaling=h,
            mi=h, correlation=h, pairwise_coverage=h,
            log="",
        )

    def to_tuple(self) -> tuple:
        return (
            self.result_url, self.output_file, self.csv_file,
            self.profile_plot, self.seq_dropdown,
            self.stats, self.load_status,
            self.mod_heatmap, self.termination, self.coverage,
            self.read_hist, self.cumulative_reads, self.snr_scaling,
            self.mi, self.correlation, self.pairwise_coverage,
            self.log,
        )


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
    compute_pairwise: bool = False


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
    sig: float = 0.05


# --- FASTA parsing ---


def _parse_fasta(fasta_path: str) -> list[tuple[str, str]]:
    """Parse FASTA file, returning list of (name, sequence) tuples."""
    entries: list[tuple[str, str]] = []
    name = ""
    seq_parts: list[str] = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq_parts:
                    entries.append((name, "".join(seq_parts)))
                    seq_parts = []
                name = line[1:].split()[0]
            elif line:
                seq_parts.append(line)
    if seq_parts:
        entries.append((name, "".join(seq_parts)))
    return entries


# --- Input validation ---


def _sanitize_group_name(raw: str | None) -> str:
    """Normalize a user-supplied group name to a safe HDF5 path component."""
    name = re.sub(r"[^\w\-]", "_", (raw or "").strip())
    return name or DEFAULT_GROUP_NAME


def _fastq_stem(path: str) -> str:
    """Return the sample name from a FASTQ path by stripping known extensions."""
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
    if cfg.compute_pairwise:
        cmd.append("--pairwise")
    cmd.extend(bam_files)
    return cmd


# --- Intermediate checks ---


def _check_bam_files(alignments_dir: str) -> list[str]:
    """Return sorted absolute BAM paths, or raise if none found."""
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


# --- CSV generation ---


def _generate_csv(
    h5_path: str,
    fasta_path: str,
    group_names: list[str],
) -> str:
    """Generate a CSV from HDF5 profiles with columns for each group."""
    fasta_entries = _parse_fasta(fasta_path)
    csv_path = h5_path.rsplit(".", 1)[0] + ".csv"

    with h5py.File(h5_path, "r") as f, open(csv_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)

        first_grp = f[group_names[0]]
        n_refs = first_grp["reactivity"].shape[0]
        seq_len = first_grp["reactivity"].shape[1]
        multi_ref = n_refs > 1

        header: list[str] = []
        if multi_ref:
            header.append("Reference")
        header.extend(["Position", "Nucleotide"])
        for gn in group_names:
            header.extend([gn, f"{gn}_error"])
        writer.writerow(header)

        for ref_idx in range(n_refs):
            ref_name = fasta_entries[ref_idx][0] if ref_idx < len(fasta_entries) else f"ref_{ref_idx + 1}"
            ref_seq = fasta_entries[ref_idx][1] if ref_idx < len(fasta_entries) else ""

            group_data = {}
            for gn in group_names:
                group_data[gn] = {
                    "reactivity": np.array(f[gn]["reactivity"])[ref_idx],
                    "error": np.array(f[gn]["error"])[ref_idx],
                }

            for pos in range(seq_len):
                row: list[str] = []
                if multi_ref:
                    row.append(ref_name)
                row.append(str(pos + 1))
                row.append(ref_seq[pos] if pos < len(ref_seq) else "")
                for gn in group_names:
                    r = group_data[gn]["reactivity"][pos]
                    e = group_data[gn]["error"][pos]
                    row.append(f"{r:.6f}" if np.isfinite(r) else "")
                    row.append(f"{e:.6f}" if np.isfinite(e) else "")
                writer.writerow(row)

    return csv_path


# --- HDF5 reading and plotting ---


def _build_plots(
    mod: cmuts.ProbingData,
    nomod: cmuts.ProbingData | None,
    combined: cmuts.ProbingData,
    name: str,
) -> dict[str, go.Figure | None]:
    """Build all diagnostic plots from in-memory ProbingData objects."""
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

    if combined.probability is not None:
        prob = np.asarray(combined.probability)
        plots["pairwise_coverage"] = plot_pairwise_coverage(prob[0, :, :, 1, 1], name)
    else:
        plots["pairwise_coverage"] = None

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


def _build_stats_table(h5_path: str, group_names: list[str]) -> list[list[str]]:
    """Build a stats table for one or more groups."""
    rows: list[list[str]] = []
    multi = len(group_names) > 1
    with h5py.File(h5_path, "r") as f:
        for gn in group_names:
            grp = f[gn] if gn in f else f
            reactivity = np.array(grp["reactivity"])
            reads = np.array(grp["reads"])
            error = np.array(grp["error"])
            snr = np.array(grp["SNR"])

            n_refs = reactivity.shape[0]
            seq_len = reactivity.shape[1]
            total_reads = int(reads.sum())
            valid = np.isfinite(reactivity)

            if multi:
                rows.append([f"--- {gn} ---", ""])

            rows.extend([
                ["References", f"{n_refs:,}"],
                ["Reference length", f"{seq_len:,}"],
                ["Total reads", f"{total_reads:,}"],
                ["Mean reads per reference", f"{np.mean(reads):,.1f}"],
                ["Median reads per reference", f"{int(np.median(reads)):,}"],
            ])

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
    "pairwise_coverage",
]


def save_results(
    h5_path: str,
    csv_path: str | None,
    group_names: list[str],
    plots: dict[str, go.Figure | None],
    stats_rows: list[list[str]],
    names: list[str],
    job_id: str | None = None,
) -> str:
    """Save pipeline results to persistent storage. Returns the job ID."""
    if job_id is None:
        job_id = uuid.uuid4().hex[:12]
    job_dir = os.path.join(RESULTS_DIR, job_id)
    os.makedirs(job_dir, exist_ok=True)

    shutil.copy(h5_path, os.path.join(job_dir, "profiles.h5"))
    if csv_path and os.path.isfile(csv_path):
        shutil.copy(csv_path, os.path.join(job_dir, "profiles.csv"))

    meta = {
        "group_names": group_names,
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


def _plot_update(fig: go.Figure | None):
    """Wrap a figure in a gr.update that shows/hides the component."""
    if fig is None:
        return gr.update(visible=False, value=None)
    return gr.update(visible=True, value=fig)


def _error_yield(msg: str) -> tuple:
    """Build a yield tuple that displays an error without breaking UI state."""
    r = ResultUpdate.hidden()
    r.log = f"Error: {msg}"
    r.load_status = msg
    return r.to_tuple()


def _progress_yield(result_url: str, log_lines: list[str]) -> tuple:
    """Build the in-progress yield tuple."""
    r = ResultUpdate.hidden()
    r.result_url = gr.update(visible=True, value=result_url)
    r.log = "\n".join(log_lines)
    return r.to_tuple()


def run_pipeline(
    fasta_file: str,
    groups: list[GroupInput],
    align_cfg: AlignConfig,
    core_cfg: CoreConfig,
    norm_cfg: NormConfig,
):
    """Run the full cmuts pipeline for one or more experiment groups.

    All groups share a single alignment and mutation-counting pass.
    Normalization is computed from pooled data across all groups so that
    reactivity values are directly comparable.
    """
    if fasta_file is None:
        yield _error_yield("Please upload a FASTA file.")
        return
    if not groups:
        yield _error_yield("Please add at least one group with a modified FASTQ file.")
        return

    for g in groups:
        for path, label in [(g.mod_fastq, f"{g.name} Modified FASTQ"),
                            (g.nomod_fastq, f"{g.name} Control FASTQ")]:
            if path is not None and _file_size_mb(path) > MAX_FASTQ_MB:
                yield _error_yield(
                    f"{label} is {_file_size_mb(path):.0f} MB. "
                    f"The free tier has limited RAM (16 GB); files over "
                    f"{MAX_FASTQ_MB} MB may cause out-of-memory errors. "
                    f"Consider downsampling first."
                )
                return

    cleanup_old_results()

    workdir = tempfile.mkdtemp(prefix="cmuts_")
    outdir = os.path.join(workdir, "outputs")
    os.makedirs(outdir)

    job_id = uuid.uuid4().hex[:12]
    space_host = os.environ.get("SPACE_HOST", "")
    base = f"https://{space_host}" if space_host else ""
    result_url = f"{base}/results/{job_id}"

    job_dir = os.path.join(RESULTS_DIR, job_id)
    os.makedirs(job_dir, exist_ok=True)
    meta = {
        "group_names": [g.name for g in groups],
        "created_at": time.time(),
        "names": [],
        "stats_rows": [],
    }
    with open(os.path.join(job_dir, "meta.json"), "w") as f:
        json.dump(meta, f)

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

        # Copy all FASTQs, prefixed by group name to avoid collisions
        fastq_dir = os.path.join(workdir, "fastq")
        os.makedirs(fastq_dir)

        group_fastq_info: list[tuple[str, str, str | None]] = []
        for g in groups:
            prefix = g.name + "__"
            mod_basename = prefix + os.path.basename(g.mod_fastq)
            mod_stem = _fastq_stem(mod_basename)
            shutil.copy(g.mod_fastq, os.path.join(fastq_dir, mod_basename))

            nomod_stem = None
            if g.nomod_fastq is not None:
                nomod_basename = prefix + os.path.basename(g.nomod_fastq)
                nomod_stem = _fastq_stem(nomod_basename)
                shutil.copy(g.nomod_fastq, os.path.join(fastq_dir, nomod_basename))

            group_fastq_info.append((g.name, mod_stem, nomod_stem))

        # Step 1: Align all FASTQs together
        log("=== Step 1: Aligning reads ===")
        yield _progress_yield(result_url, log_lines)

        fastq_files = sorted(glob.glob(os.path.join(fastq_dir, "*")))
        alignments_dir = os.path.join(outdir, "alignments")
        align_cmd = _build_align_cmd(fasta_path, alignments_dir, fastq_files, align_cfg)
        if not run(align_cmd, cwd=outdir):
            yield _progress_yield(result_url, log_lines)
            return
        yield _progress_yield(result_url, log_lines)

        # Step 2: Count mutations for all BAMs
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

        # Step 3: Normalize reactivities with shared normalization factor
        log("\n=== Step 3: Normalizing reactivities ===")
        yield _progress_yield(result_url, log_lines)

        counts_path = os.path.join(outdir, counts_h5)

        cmuts_groups = [
            cmuts.Group(
                name=name,
                mod=[f"alignments/{mod_stem}"],
                nomod=[f"alignments/{nomod_stem}"] if nomod_stem else None,
            )
            for name, mod_stem, nomod_stem in group_fastq_info
        ]

        norm_opts = cmuts.Opts(
            cmuts.DataGroups([]),
            cmuts.DataGroups(None),
            norm_cfg.blank_cutoff,
            not norm_cfg.no_insertions,
            not norm_cfg.no_deletions,
            norm_cfg.norm_method,
            (norm_cfg.blank_5p, norm_cfg.blank_3p),
            (norm_cfg.clip_low, norm_cfg.clip_high),
            norm_cfg.sig,
        )

        with h5py.File(counts_path, "r") as f:
            results = cmuts.compute_reactivities(
                f, fasta_path, cmuts_groups, norm_opts, shared_norm=True,
            )

        if len(results) > 1:
            log(f"  Pooled {norm_cfg.norm_method} normalization across {len(results)} groups.")
        log("Normalization complete.")

        # Save multi-group HDF5
        final_path = os.path.join(outdir, "profiles.h5")
        cmuts.save_groups(final_path, [(r.group.name, r.combined) for r in results])

        # Generate CSV
        all_group_names = [r.group.name for r in results]
        csv_path = _generate_csv(final_path, fasta_path, all_group_names)
        log("Generated CSV output.")

        # Build plots
        first = results[0]
        first_name = first.group.name

        # Profile: overlay all groups for single-reference data
        if len(results) > 1 and first.combined.single():
            reactivities = [np.asarray(r.combined.reactivity)[0] for r in results]
            profile_fig = plot_profiles(reactivities, all_group_names)
        else:
            profile_fig = plot_examples(
                np.asarray(first.combined.reactivity),
                np.asarray(first.combined.error),
                first_name,
            )

        # Diagnostic plots from first group
        diag_plots = _build_plots(first.mod, first.nomod, first.combined, first_name)
        diag_plots["profile"] = profile_fig

        # Sequence names for dropdown (shared across groups)
        reactivity, names = _read_profiles(final_path, first_name)
        dropdown_update = gr.Dropdown(
            choices=names, value=names[0],
            visible=(len(names) > 1),
        )

        # Stats for all groups
        stats_rows = _build_stats_table(final_path, all_group_names)

        save_results(
            final_path, csv_path, all_group_names,
            diag_plots, stats_rows, names, job_id=job_id,
        )

        log(f"\nDone. Generated profiles for {len(results)} group(s).")
        log(f"Results available at: {result_url} (expires in {RESULTS_TTL_HOURS}h)")
        yield ResultUpdate(
            result_url=gr.update(visible=True, value=result_url),
            output_file=gr.update(visible=True, value=final_path),
            csv_file=gr.update(visible=True, value=csv_path),
            profile_plot=_plot_update(diag_plots["profile"]),
            seq_dropdown=dropdown_update,
            stats=gr.update(visible=True, value=stats_rows),
            mod_heatmap=_plot_update(diag_plots["mod_heatmap"]),
            termination=_plot_update(diag_plots["termination"]),
            coverage=_plot_update(diag_plots["coverage"]),
            read_hist=_plot_update(diag_plots["read_hist"]),
            cumulative_reads=_plot_update(diag_plots["cumulative_reads"]),
            snr_scaling=_plot_update(diag_plots["snr_scaling"]),
            mi=_plot_update(diag_plots["mi"]),
            correlation=_plot_update(diag_plots["correlation"]),
            pairwise_coverage=_plot_update(diag_plots["pairwise_coverage"]),
            log="\n".join(log_lines),
        ).to_tuple()

    except subprocess.TimeoutExpired:
        log(f"Pipeline timed out ({PIPELINE_TIMEOUT_SEC // 60} minute limit).")
        yield _progress_yield(result_url, log_lines)
    except Exception as e:
        log(f"Error: {e}")
        log(traceback.format_exc())
        yield _progress_yield(result_url, log_lines)


# --- Gradio callbacks ---


def _run_pipeline_gradio(fasta_file, n_groups, *all_args):
    """Gradio-facing wrapper: unpacks flat args into structured input."""
    gn_vals = all_args[:MAX_GROUPS]
    mod_vals = all_args[MAX_GROUPS : 2 * MAX_GROUPS]
    nomod_vals = all_args[2 * MAX_GROUPS : 3 * MAX_GROUPS]
    cfg = all_args[3 * MAX_GROUPS :]

    (norm_method, no_insertions, no_deletions, clip_low, clip_high,
     trim_5, trim_3, local_align,
     min_mapq, min_phred, min_length, max_length, no_mismatches, strand,
     blank_5p, blank_3p, blank_cutoff, norm_cutoff, norm_percentile,
     compute_pairwise, sig) = cfg

    groups: list[GroupInput] = []
    for i in range(int(n_groups)):
        if mod_vals[i] is not None:
            gn = _sanitize_group_name(gn_vals[i]) or f"group_{i + 1}"
            groups.append(GroupInput(gn, mod_vals[i], nomod_vals[i]))

    yield from run_pipeline(
        fasta_file=fasta_file,
        groups=groups,
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
            compute_pairwise=compute_pairwise,
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
            sig=float(sig or 0.05),
        ),
    )


def select_profile(seq_name: str, output_file: str) -> tuple:
    """Switch the displayed profile and pairwise plots for the selected sequence."""
    if not output_file or not seq_name:
        empty = go.Figure()
        hidden = gr.update(visible=False, value=None)
        return empty, hidden, hidden, hidden

    with h5py.File(output_file, "r") as f:
        group_names = sorted(k for k in f.keys() if k != "sequence")
        if not group_names:
            hidden = gr.update(visible=False, value=None)
            return go.Figure(), hidden, hidden, hidden

        first_grp = f[group_names[0]]
        reactivity = np.array(first_grp["reactivity"])
        sequences = None
        if "sequence" in f:
            sequences = [
                s.decode() if isinstance(s, bytes) else s for s in f["sequence"]
            ]

        names: list[str] = []
        for i in range(reactivity.shape[0]):
            seq = sequences[i] if sequences and i < len(sequences) else None
            if seq and len(seq) > 50:
                names.append(seq[:50] + "...")
            elif seq:
                names.append(seq)
            else:
                names.append(f"Sequence {i + 1}")

        try:
            idx = names.index(seq_name)
        except ValueError:
            idx = 0

        if len(group_names) > 1:
            reactivities = [np.array(f[gn]["reactivity"])[idx] for gn in group_names]
            profile_fig = plot_profiles(reactivities, group_names)
        else:
            error = np.array(first_grp["error"])
            profile_fig = plot_profile(reactivity[idx], error[idx], names[idx])

        mi_fig = None
        corr_fig = None
        if "mutual-information" in first_grp:
            mi_fig = plot_mi(np.array(first_grp["mutual-information"])[idx], names[idx])
        if "covariance" in first_grp:
            corr_fig = plot_correlation(np.array(first_grp["covariance"])[idx], names[idx])

    return profile_fig, _plot_update(mi_fig), _plot_update(corr_fig), _plot_update(None)


def load_example():
    """Load bundled example files into group 1."""
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
    if not job_id:
        return ResultUpdate.hidden().to_tuple()

    job_dir = os.path.join(RESULTS_DIR, job_id)
    meta_path = os.path.join(job_dir, "meta.json")

    if not os.path.isdir(job_dir):
        r = ResultUpdate.hidden()
        r.load_status = (
            f"Result not found. It may have expired "
            f"(results are kept for {RESULTS_TTL_HOURS} hours)."
        )
        return r.to_tuple()

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

    # Handle both old (group_name) and new (group_names) meta formats
    group_names = meta.get("group_names")
    if group_names is None:
        group_names = [meta.get("group_name", DEFAULT_GROUP_NAME)]
    names = meta.get("names", [])

    space_host = os.environ.get("SPACE_HOST", "")
    base = f"https://{space_host}" if space_host else ""

    h5_path = os.path.join(job_dir, "profiles.h5")
    csv_path = os.path.join(job_dir, "profiles.csv")
    hidden = gr.update(visible=False, value=None)

    return ResultUpdate(
        result_url=gr.update(visible=True, value=f"{base}/results/{job_id}"),
        output_file=gr.update(visible=True, value=h5_path) if os.path.isfile(h5_path) else hidden,
        csv_file=gr.update(visible=True, value=csv_path) if os.path.isfile(csv_path) else hidden,
        profile_plot=_plot_update(loaded.get("profile")),
        seq_dropdown=gr.Dropdown(choices=names, value=names[0] if names else None, visible=len(names) > 1),
        stats=gr.update(visible=True, value=meta.get("stats_rows", [])),
        mod_heatmap=_plot_update(loaded.get("mod_heatmap")),
        termination=_plot_update(loaded.get("termination")),
        coverage=_plot_update(loaded.get("coverage")),
        read_hist=_plot_update(loaded.get("read_hist")),
        cumulative_reads=_plot_update(loaded.get("cumulative_reads")),
        snr_scaling=_plot_update(loaded.get("snr_scaling")),
        mi=_plot_update(loaded.get("mi")),
        correlation=_plot_update(loaded.get("correlation")),
        pairwise_coverage=_plot_update(loaded.get("pairwise_coverage")),
    ).to_tuple()


# --- Gradio UI ---

with gr.Blocks(title="cmuts — RNA Chemical Probing Analysis") as demo:
    gr.Markdown(
        """
        # cmuts — RNA Chemical Probing Analysis

        Upload a FASTA reference and FASTQ file(s) from a MaP-seq experiment
        to compute normalized reactivity profiles.  Use **+ Add group** to
        compare multiple conditions (e.g. with/without ligand) — normalization
        is applied across all groups so values are directly comparable.

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
        with gr.Column() as input_section:
            gr.Markdown("### Input data")
            fasta_input = gr.File(label="Reference FASTA", file_types=[".fasta", ".fa"])

            n_groups_state = gr.State(1)

            # Dynamic group inputs (up to MAX_GROUPS)
            group_rows: list[gr.Group] = []
            group_names_inputs: list[gr.Textbox] = []
            mod_inputs: list[gr.File] = []
            nomod_inputs: list[gr.File] = []

            for _i in range(MAX_GROUPS):
                with gr.Group(visible=(_i == 0)) as row:
                    with gr.Row():
                        gn = gr.Textbox(
                            label=f"Group {_i + 1} name",
                            value="experiment" if _i == 0 else "",
                            placeholder="e.g. 2A3_with_cdiGMP",
                            scale=1,
                        )
                        mod = gr.File(
                            label=f"Group {_i + 1} Modified FASTQ (required)",
                            file_types=[".fastq", ".fq", ".gz"],
                            scale=2,
                        )
                        nomod = gr.File(
                            label=f"Group {_i + 1} Control FASTQ (optional)",
                            file_types=[".fastq", ".fq", ".gz"],
                            scale=2,
                        )
                group_rows.append(row)
                group_names_inputs.append(gn)
                mod_inputs.append(mod)
                nomod_inputs.append(nomod)

            with gr.Row():
                add_group_btn = gr.Button("+ Add group", variant="secondary", size="sm")
                remove_group_btn = gr.Button("- Remove last group", variant="secondary", size="sm")
                example_btn = gr.Button("Load example data", variant="secondary", size="sm")

            def _add_group(n):
                n = min(n + 1, MAX_GROUPS)
                return (n,) + tuple(gr.update(visible=(i < n)) for i in range(MAX_GROUPS))

            def _remove_group(n):
                n = max(n - 1, 1)
                return (n,) + tuple(gr.update(visible=(i < n)) for i in range(MAX_GROUPS))

            add_group_btn.click(
                fn=_add_group,
                inputs=[n_groups_state],
                outputs=[n_groups_state] + group_rows,
            )
            remove_group_btn.click(
                fn=_remove_group,
                inputs=[n_groups_state],
                outputs=[n_groups_state] + group_rows,
            )

            gr.Markdown("### Options")
            with gr.Accordion("Alignment", open=False):
                gr.Markdown(
                    "*If no adapter sequences are provided, cmuts will attempt to "
                    "recognize adapters outside the reference sequence and auto-trim.*"
                )
                with gr.Row():
                    trim_5 = gr.Textbox(label="5' adapter to trim (optional)", placeholder="auto-detected if blank")
                    trim_3 = gr.Textbox(label="3' adapter to trim (optional)", placeholder="auto-detected if blank")
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

            with gr.Accordion("Pairwise analysis", open=False):
                gr.Markdown(
                    "Compute pairwise modification correlations (mutual information "
                    "and Pearson correlation). Cost is O(L&sup2;) in sequence length, "
                    "so this is slow for long references."
                )
                compute_pairwise = gr.Checkbox(label="Compute pairwise correlations", value=False)
                sig = gr.Slider(
                    minimum=0.001, maximum=0.1, step=0.001, value=0.05,
                    label="Significance threshold (Bonferroni-corrected)",
                )

            run_btn = gr.Button("Run Pipeline", variant="primary")

        gr.Markdown("### Results")
        result_url = gr.Textbox(
            label=f"Result link (bookmark this — expires in {RESULTS_TTL_HOURS}h)",
            interactive=False,
            visible=False,
        )
        output_file = gr.File(label="Output HDF5", visible=False)
        csv_file = gr.File(label="Output CSV", visible=False)
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
        pairwise_coverage_plot = gr.Plot(label="Pairwise Coverage", visible=False)
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
            outputs=[fasta_input, mod_inputs[0], nomod_inputs[0], group_names_inputs[0]],
        )

        # Shared outputs list matching ResultUpdate.to_tuple() field order
        _result_outputs = [
            result_url, output_file, csv_file,
            output_plot, seq_dropdown,
            output_stats, load_status,
            mod_heatmap_plot, termination_plot, coverage_plot,
            read_hist_plot, cumulative_reads_plot, snr_scaling_plot,
            mi_plot, correlation_plot, pairwise_coverage_plot,
            output_log,
        ]

        run_btn.click(
            fn=_run_pipeline_gradio,
            inputs=[
                fasta_input,
                n_groups_state,
                # Group inputs: names, then mod FASTQs, then nomod FASTQs
                *group_names_inputs,
                *mod_inputs,
                *nomod_inputs,
                # Config options
                norm_method,
                no_insertions,
                no_deletions,
                clip_low,
                clip_high,
                trim_5,
                trim_3,
                local_align,
                min_mapq,
                min_phred,
                min_length,
                max_length,
                no_mismatches,
                strand,
                blank_5p,
                blank_3p,
                blank_cutoff,
                norm_cutoff,
                norm_percentile,
                compute_pairwise,
                sig,
            ],
            outputs=_result_outputs,
        )

        seq_dropdown.change(
            fn=select_profile,
            inputs=[seq_dropdown, output_file],
            outputs=[output_plot, mi_plot, correlation_plot, pairwise_coverage_plot],
        )

        load_btn.click(
            fn=load_saved_result,
            inputs=[prev_job_id],
            outputs=_result_outputs,
        )

    _load_outputs = [input_section] + _result_outputs

    def _load_from_query(request: gr.Request):
        """Auto-load results when ?job_id= is present in the URL."""
        job_id = (request.query_params.get("job_id") or "").strip()
        if not job_id:
            return [gr.update()] * len(_load_outputs)
        result = load_saved_result(job_id)
        return (gr.update(visible=False),) + result

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

            - **100-200x faster** than ShapeMapper2 and RNAframework. A dataset
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
            f"""
            ## Quick Start

            1. Click **Load example data** on the Run tab to populate the inputs
               with a bundled dataset.
            2. Leave the default settings and click **Run Pipeline**.
            3. The pipeline runs three steps — alignment, mutation counting, and
               normalization — and streams its progress to the log.
            4. When finished, download the output HDF5 or CSV file and explore
               the interactive reactivity profile. If multiple reference
               sequences are present, use the dropdown to switch between them.

            ## Inputs

            | Field | Description |
            |-------|-------------|
            | **Reference FASTA** | One or more RNA sequences in FASTA format. Each sequence is treated as a separate reference for alignment. |
            | **Modified FASTQ** | Reads from the chemically treated condition (e.g., DMS, 2A3, SHAPE). Compressed `.fastq.gz` is accepted. |
            | **Control FASTQ** | *(Optional)* Reads from the untreated/DMSO condition. Providing a control enables background subtraction during normalization. |
            | **Group name** | A label for the experiment group, used in the output files. |

            ### Multiple groups

            Use **+ Add group** to compare multiple experimental conditions
            (e.g. with/without a ligand). Each group has its own modified and
            control FASTQ files. All groups share the same reference FASTA.

            Normalization is applied **across all groups** so that reactivity
            values are directly comparable. The profile plot overlays all
            groups for easy visual comparison.

            ## Adapter trimming

            If no adapter sequences are provided in the Alignment settings,
            cmuts will attempt to recognize adapter sequences outside the
            reference and auto-trim them.

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

            The interactive chart shows per-nucleotide reactivity values.
            Hover over any position to see the exact position, nucleotide
            identity, and reactivity value. Peaks correspond to unpaired or
            flexible nucleotides; low/near-zero regions correspond to
            base-paired or otherwise protected positions.

            When multiple groups are present, the profile overlays all groups
            for direct comparison.

            ### Output files

            **HDF5** (`profiles.h5`) contains normalized per-nucleotide
            reactivity profiles organized by group. Open it in Python with:

            ```python
            import h5py

            with h5py.File("profiles.h5") as f:
                for group in f:
                    if group == "sequence":
                        continue
                    reactivities = f[group]["reactivity"][:]
                    print(group, reactivities.shape)
            ```

            **CSV** (`profiles.csv`) contains the same data in a flat table
            with columns: Position, Nucleotide, and for each group:
            group_name, group_name_error. This is convenient for spreadsheet
            analysis or cross-checking with other tools.

            ### Result link

            After the pipeline completes, a **result link** is displayed that
            you can bookmark or share. Opening the link loads the results
            directly into the app with the interactive plots, summary
            statistics, and file downloads. Results are stored for
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
            statistics) are stored for **{RESULTS_TTL_HOURS} hours** to provide
            bookmarkable result links, then automatically deleted. No user
            accounts, tracking, or cookies are used.
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
    group_names = meta.get("group_names", [meta.get("group_name", "profiles")])
    filename = "-".join(group_names) + "-profiles.h5"

    return FileResponse(
        h5_path,
        media_type="application/x-hdf5",
        filename=filename,
    )


@app.get("/results/{job_id}/download/csv")
async def results_download_csv(job_id: str):
    csv_path = os.path.join(RESULTS_DIR, job_id, "profiles.csv")
    if not os.path.isfile(csv_path):
        return HTMLResponse(
            "<h1>File not found</h1><p>CSV not available for this result.</p>",
            status_code=404,
        )

    with open(os.path.join(RESULTS_DIR, job_id, "meta.json")) as f:
        meta = json.load(f)
    group_names = meta.get("group_names", [meta.get("group_name", "profiles")])
    filename = "-".join(group_names) + "-profiles.csv"

    return FileResponse(
        csv_path,
        media_type="text/csv",
        filename=filename,
    )


# Mount Gradio onto the FastAPI app
app = gr.mount_gradio_app(app, demo, path="")

# Run cleanup on startup
cleanup_old_results()

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=7860)
