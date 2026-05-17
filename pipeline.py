"""Pipeline orchestration for the cmuts web app.

Framework-agnostic. Takes input paths + configs, runs ``cmuts align`` /
``core`` / ``compute_reactivities``, and writes everything for one job
into a single directory under ``RESULTS_DIR``.

On-disk job layout::

    {job_dir}/
        profiles.h5               final per-group HDF5
        profiles.csv              flat per-position table
        meta.json                 status, group/sequence names, stats, log
        log.txt                   raw streaming log
        defattr/{group}.defattr   (if a CIF was supplied)
        groups/{group}/{key}.json plotly figure JSON, first reference
        combined_profile.json     overlay across groups (single-ref + multi-group)
"""

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
from dataclasses import dataclass, field
from typing import Callable

import h5py
import numpy as np
import plotly.graph_objects as go


# --- Constants ---

MAX_GROUPS = int(os.environ.get("CMUTS_MAX_GROUPS", "5"))
MAX_FASTQ_MB = int(os.environ.get("CMUTS_MAX_FASTQ_MB", "500"))
RESULTS_TTL_HOURS = int(os.environ.get("CMUTS_RESULTS_TTL_HOURS", "48"))
PIPELINE_TIMEOUT_SEC = int(os.environ.get("CMUTS_PIPELINE_TIMEOUT_SEC", "600"))
DEFAULT_GROUP_NAME = "profile"

_default_results_dir = "/data/results" if os.path.isdir("/data") else "/tmp/cmuts_results"
RESULTS_DIR = os.environ.get("CMUTS_RESULTS_DIR", _default_results_dir)
os.makedirs(RESULTS_DIR, exist_ok=True)

_FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")

PLOT_KEYS = [
    "profile", "mod_heatmap", "termination", "coverage",
    "read_hist", "cumulative_reads", "snr_scaling",
    "mi", "correlation", "pairwise_coverage",
]


# --- Dataclasses ---


@dataclass
class GroupInput:
    name: str
    mod_fastq: str
    nomod_fastq: str | None = None


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
    clip_below: float | None = None
    clip_above: float | None = None
    blank_5p: int = 0
    blank_3p: int = 0
    blank_cutoff: int = 10
    norm_cutoff: int = 500
    norm_percentile: int = 90
    sig: float = 0.05


@dataclass
class JobState:
    """Mutable in-memory state for an in-flight job."""
    job_id: str
    status: str = "running"  # "running" | "done" | "error"
    log_lines: list[str] = field(default_factory=list)
    error: str | None = None

    def log(self, msg: str) -> None:
        self.log_lines.append(msg)


# --- FASTA / FASTQ helpers ---


def parse_fasta(fasta_path: str) -> list[tuple[str, str]]:
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


def sanitize_group_name(raw: str | None) -> str:
    name = re.sub(r"[^\w\-]", "_", (raw or "").strip())
    return name or DEFAULT_GROUP_NAME


def fastq_stem(path: str) -> str:
    name = os.path.basename(path)
    for suffix in _FASTQ_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return os.path.splitext(name)[0]


def file_size_mb(path: str) -> float:
    return os.path.getsize(path) / (1024 * 1024)


def build_seq_names(n: int, sequences: list[str] | None) -> list[str]:
    """Build display labels for the sequence selector, disambiguating
    truncated duplicates."""
    raw: list[str] = []
    for i in range(n):
        seq = sequences[i] if sequences and i < len(sequences) else None
        if seq and len(seq) > 50:
            raw.append(seq[:50] + "...")
        elif seq:
            raw.append(seq)
        else:
            raw.append(f"Sequence {i + 1}")
    counts: dict[str, int] = {}
    out: list[str] = []
    for label in raw:
        if raw.count(label) > 1:
            counts[label] = counts.get(label, 0) + 1
            out.append(f"{label} (#{counts[label]})")
        else:
            out.append(label)
    return out


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


# --- CSV generation ---


def _generate_csv(
    h5_path: str,
    fasta_path: str,
    group_names: list[str],
) -> str:
    """Generate a CSV from HDF5 profiles with columns for each group."""
    fasta_entries = parse_fasta(fasta_path)
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

            group_data: dict[str, dict[str, np.ndarray]] = {}
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


# --- Structure visualization (ChimeraX defattr) ---


def _build_defattrs(
    cif_path: str,
    sequence: str,
    results: list,
    out_dir: str,
    chimerax_bin: str = "ChimeraX",
):
    """Generate one defattr per group plus a markdown command snippet.

    Returns (defattr_paths, markdown).
    """
    import cmuts as _cmuts  # imported lazily so app starts without cmuts

    cif_basename = os.path.basename(cif_path)
    aln_seq = sequence.upper().replace("U", "T")

    defattr_paths: list[str] = []
    blocks: list[str] = [
        "### Visualize the structure with ChimeraX",
        "",
        f"Download each `.defattr` file below, place it next to your "
        f"`{cif_basename}` (a copy of your uploaded structure), and run "
        f"the matching command in ChimeraX's command line.",
        "",
    ]

    for r in results:
        name = r.group.name
        reactivity = np.asarray(r.combined.reactivity)
        if reactivity.shape[0] != 1:
            continue
        defattr_path = os.path.join(out_dir, f"{name}.defattr")
        try:
            max_value = _cmuts.visualize.make_defattr(
                reactivity[0], aln_seq, cif_path, defattr_path,
            )
        except Exception as e:  # noqa: BLE001
            blocks.append(f"**{name}:** could not generate defattr — {e}")
            blocks.append("")
            continue
        cmd = _cmuts.visualize.chimerax_command(
            cif_basename, os.path.basename(defattr_path),
            color="indianred", max_value=max_value,
        )
        blocks.append(f"**{name}:**")
        blocks.append("```")
        blocks.append(f"{chimerax_bin} --cmd '{cmd}'")
        blocks.append("```")
        blocks.append("")
        defattr_paths.append(defattr_path)

    if not defattr_paths:
        return [], ""
    return defattr_paths, "\n".join(blocks)


# --- Plot building ---


def _build_plots_for_group(
    mod, nomod, combined,
    group_name: str,
    sequence: str | None,
) -> dict[str, go.Figure | None]:
    """Build every diagnostic plot for one group. Returns key -> Figure or
    None when the underlying data is absent."""
    from cmuts.visualize.plotly import (
        plot_correlation, plot_coverage, plot_cumulative_reads, plot_examples,
        plot_heatmap, plot_mi, plot_pairwise_coverage,
        plot_read_hist, plot_snr_scaling, plot_termination,
    )

    plots: dict[str, go.Figure | None] = {}
    plots["profile"] = plot_examples(
        np.asarray(combined.reactivity), np.asarray(combined.error),
        group_name, sequence=sequence,
    )
    plots["mod_heatmap"] = plot_heatmap(np.asarray(combined.heatmap), group_name)
    plots["termination"] = plot_termination(np.asarray(combined.terminations), group_name)
    plots["coverage"] = plot_coverage(
        np.asarray(combined.coverage), np.asarray(combined.reads), group_name,
    )

    is_multi = not combined.single()
    reads = np.asarray(combined.reads)
    plots["read_hist"] = plot_read_hist(reads, group_name) if is_multi else None
    plots["cumulative_reads"] = plot_cumulative_reads(reads, group_name) if is_multi else None

    plots["snr_scaling"] = plot_snr_scaling(mod, nomod, combined, group_name)

    plots["mi"] = (
        plot_mi(np.asarray(combined.mi)[0], group_name)
        if combined.mi is not None else None
    )
    plots["correlation"] = (
        plot_correlation(np.asarray(combined.covariance)[0], group_name)
        if combined.covariance is not None else None
    )
    if combined.probability is not None:
        prob = np.asarray(combined.probability)
        plots["pairwise_coverage"] = plot_pairwise_coverage(prob[0, :, :, 1, 1], group_name)
    else:
        plots["pairwise_coverage"] = None
    return plots


def _save_plot_json(fig: go.Figure | None, path: str) -> bool:
    if fig is None:
        return False
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(fig.to_json())
    return True


def _build_stats_for_group(grp, ref_count: int) -> list[list[str]]:
    reactivity = np.array(grp["reactivity"])
    reads = np.array(grp["reads"])
    error = np.array(grp["error"])
    snr = np.array(grp["SNR"])

    n_refs = reactivity.shape[0]
    seq_len = reactivity.shape[1]
    total_reads = int(reads.sum())
    valid = np.isfinite(reactivity)

    rows: list[list[str]] = [
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


# --- Results persistence + cleanup ---


def job_dir_for(job_id: str) -> str:
    return os.path.join(RESULTS_DIR, job_id)


def cleanup_old_results() -> None:
    """Delete result directories older than RESULTS_TTL_HOURS."""
    cutoff = time.time() - RESULTS_TTL_HOURS * 3600
    if not os.path.isdir(RESULTS_DIR):
        return
    for entry in os.scandir(RESULTS_DIR):
        if not entry.is_dir():
            continue
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


def write_meta(job_dir: str, meta: dict) -> None:
    with open(os.path.join(job_dir, "meta.json"), "w") as f:
        json.dump(meta, f)


def read_meta(job_dir: str) -> dict | None:
    path = os.path.join(job_dir, "meta.json")
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        return json.load(f)


def write_log(job_dir: str, log_lines: list[str]) -> None:
    with open(os.path.join(job_dir, "log.txt"), "w") as f:
        f.write("\n".join(log_lines))


# --- On-demand plot generation (sequence switching) ---


def build_profile_plot(job_dir: str, group_name: str, seq_idx: int) -> str | None:
    """Build a profile plot JSON on demand for a given group/sequence."""
    from cmuts.visualize.plotly import plot_profile

    h5_path = os.path.join(job_dir, "profiles.h5")
    if not os.path.isfile(h5_path):
        return None
    with h5py.File(h5_path, "r") as f:
        if group_name not in f:
            return None
        grp = f[group_name]
        reactivity = np.array(grp["reactivity"])
        error = np.array(grp["error"])
        sequences = None
        if "sequence" in f:
            sequences = [
                s.decode() if isinstance(s, bytes) else s for s in f["sequence"]
            ]
    if seq_idx < 0 or seq_idx >= reactivity.shape[0]:
        return None
    seq = sequences[seq_idx] if sequences and seq_idx < len(sequences) else None
    names = build_seq_names(reactivity.shape[0], sequences)
    fig = plot_profile(reactivity[seq_idx], error[seq_idx], names[seq_idx], sequence=seq)
    return fig.to_json()


def build_perref_plot(
    job_dir: str, group_name: str, key: str, seq_idx: int,
) -> str | None:
    """Rebuild a per-reference plot (mi / correlation / pairwise_coverage)
    for a different reference. Returns None if not available."""
    from cmuts.visualize.plotly import (
        plot_correlation, plot_mi, plot_pairwise_coverage,
    )

    h5_path = os.path.join(job_dir, "profiles.h5")
    if not os.path.isfile(h5_path):
        return None
    with h5py.File(h5_path, "r") as f:
        if group_name not in f:
            return None
        grp = f[group_name]
        if key == "mi":
            if "mutual-information" not in grp:
                return None
            arr = np.array(grp["mutual-information"])
        elif key == "correlation":
            if "covariance" not in grp:
                return None
            arr = np.array(grp["covariance"])
        elif key == "pairwise_coverage":
            if "probability" not in grp:
                return None
            arr = np.array(grp["probability"])
        else:
            return None
    if seq_idx < 0 or seq_idx >= arr.shape[0]:
        return None
    if key == "mi":
        fig = plot_mi(arr[seq_idx], group_name)
    elif key == "correlation":
        fig = plot_correlation(arr[seq_idx], group_name)
    else:  # pairwise_coverage
        fig = plot_pairwise_coverage(arr[seq_idx, :, :, 1, 1], group_name)
    return fig.to_json()


# --- Main pipeline ---


def run_pipeline(
    job_id: str,
    job_dir: str,
    fasta_path: str,
    groups: list[GroupInput],
    align_cfg: AlignConfig,
    core_cfg: CoreConfig,
    norm_cfg: NormConfig,
    cif_path: str | None,
    state: JobState,
) -> None:
    """Run the full pipeline. All status/log goes through ``state``.

    Writes results into ``job_dir``. On failure, ``state.status`` becomes
    ``"error"`` and ``state.error`` holds a short message. On success,
    ``state.status`` becomes ``"done"``.
    """
    workdir = tempfile.mkdtemp(prefix="cmuts_")
    outdir = os.path.join(workdir, "outputs")
    os.makedirs(outdir)

    def log(msg: str) -> None:
        state.log(msg)

    def run_subprocess(cmd: list[str], cwd: str = outdir) -> bool:
        log(f"$ {' '.join(cmd)}")
        try:
            result = subprocess.run(
                cmd, cwd=cwd, capture_output=True, text=True,
                timeout=PIPELINE_TIMEOUT_SEC,
            )
        except subprocess.TimeoutExpired:
            log(f"Command timed out after {PIPELINE_TIMEOUT_SEC}s")
            return False
        if result.stdout:
            log(result.stdout.rstrip())
        if result.stderr:
            log(result.stderr.rstrip())
        if result.returncode != 0:
            log(f"Command failed with exit code {result.returncode}")
            return False
        return True

    try:
        try:
            import cmuts as _cmuts
        except ImportError as e:
            raise RuntimeError(
                "cmuts is not installed in this environment. "
                "Run the app under Docker (see Dockerfile) or install cmuts."
            ) from e

        # Stage FASTQ files with a group-name prefix to avoid collisions.
        fastq_dir = os.path.join(workdir, "fastq")
        os.makedirs(fastq_dir)

        group_fastq_info: list[tuple[str, str, str | None]] = []
        for g in groups:
            prefix = g.name + "__"
            mod_basename = prefix + os.path.basename(g.mod_fastq)
            mod_stem = fastq_stem(mod_basename)
            shutil.copy(g.mod_fastq, os.path.join(fastq_dir, mod_basename))

            nomod_stem = None
            if g.nomod_fastq is not None:
                nomod_basename = prefix + os.path.basename(g.nomod_fastq)
                nomod_stem = fastq_stem(nomod_basename)
                shutil.copy(g.nomod_fastq, os.path.join(fastq_dir, nomod_basename))

            group_fastq_info.append((g.name, mod_stem, nomod_stem))

        # Step 1: align
        log("=== Step 1: Aligning reads ===")
        fastq_files = sorted(glob.glob(os.path.join(fastq_dir, "*")))
        alignments_dir = os.path.join(outdir, "alignments")
        if not run_subprocess(_build_align_cmd(fasta_path, alignments_dir, fastq_files, align_cfg)):
            raise RuntimeError("Alignment failed. See log for details.")

        bam_files_abs = sorted(glob.glob(os.path.join(alignments_dir, "*.bam")))
        if not bam_files_abs:
            raise RuntimeError(
                "Alignment produced no BAM files. Check log — the reference "
                "may not match the reads, or the FASTQ may be empty."
            )

        # Step 2: count mutations
        log("\n=== Step 2: Counting mutations ===")
        bam_files = sorted(os.path.relpath(p, outdir) for p in bam_files_abs)
        counts_h5 = "counts.h5"
        if not run_subprocess(_build_core_cmd(fasta_path, counts_h5, bam_files, core_cfg)):
            raise RuntimeError("Mutation counting failed. See log for details.")

        counts_path = os.path.join(outdir, counts_h5)
        if not os.path.isfile(counts_path):
            raise RuntimeError("cmuts core did not produce counts.h5.")

        # Step 3: normalize
        log("\n=== Step 3: Normalizing reactivities ===")
        cmuts_groups = [
            _cmuts.Group(
                name=name,
                mod=[f"alignments/{mod_stem}"],
                nomod=[f"alignments/{nomod_stem}"] if nomod_stem else None,
            )
            for name, mod_stem, nomod_stem in group_fastq_info
        ]

        norm_opts = _cmuts.Opts(
            _cmuts.DataGroups([]),
            _cmuts.DataGroups(None),
            norm_cfg.blank_cutoff,
            not norm_cfg.no_insertions,
            not norm_cfg.no_deletions,
            norm_cfg.norm_method,
            (norm_cfg.blank_5p, norm_cfg.blank_3p),
            (norm_cfg.clip_below, norm_cfg.clip_above),
            norm_cfg.sig,
        )

        with h5py.File(counts_path, "r") as f:
            results = _cmuts.compute_reactivities(
                f, fasta_path, cmuts_groups, norm_opts, shared_norm=True,
            )

        if len(results) > 1:
            log(f"  Pooled {norm_cfg.norm_method} normalization across {len(results)} groups.")
        log("Normalization complete.")

        # Save the combined HDF5 + CSV.
        final_h5 = os.path.join(job_dir, "profiles.h5")
        _cmuts.save_groups(final_h5, [(r.group.name, r.combined) for r in results])

        group_names = [r.group.name for r in results]
        csv_path = _generate_csv(final_h5, fasta_path, group_names)
        # Move CSV next to HDF5.
        final_csv = os.path.join(job_dir, "profiles.csv")
        if csv_path != final_csv:
            shutil.move(csv_path, final_csv)
        log("Wrote profiles.h5 and profiles.csv.")

        # Sequence names + single-ref check.
        fasta_entries = parse_fasta(fasta_path)
        first_combined = results[0].combined
        single_ref = bool(first_combined.single()) and len(fasta_entries) == 1
        ref_sequence = fasta_entries[0][1] if single_ref else None
        sequence_names_per_group: dict[str, list[str]] = {}

        # Per-group plots: profile, heatmap, termination, coverage, ... (first ref).
        for r in results:
            gname = r.group.name
            group_plot_dir = os.path.join(job_dir, "groups", gname)
            os.makedirs(group_plot_dir, exist_ok=True)
            plots = _build_plots_for_group(
                r.mod, r.nomod, r.combined, gname, sequence=ref_sequence,
            )
            for key, fig in plots.items():
                _save_plot_json(fig, os.path.join(group_plot_dir, f"{key}.json"))

        # Combined profile across groups (only when single-ref + >1 group).
        has_combined = False
        if single_ref and len(results) > 1:
            from cmuts.visualize.plotly import plot_profiles
            reactivities = [np.asarray(r.combined.reactivity)[0] for r in results]
            combined_fig = plot_profiles(reactivities, group_names, sequence=ref_sequence)
            _save_plot_json(combined_fig, os.path.join(job_dir, "combined_profile.json"))
            has_combined = True

        # Per-group sequence names + stats.
        stats: dict[str, list[list[str]]] = {}
        with h5py.File(final_h5, "r") as f:
            sequences = None
            if "sequence" in f:
                sequences = [
                    s.decode() if isinstance(s, bytes) else s for s in f["sequence"]
                ]
            for gname in group_names:
                reactivity = np.array(f[gname]["reactivity"])
                seq_names = build_seq_names(reactivity.shape[0], sequences)
                sequence_names_per_group[gname] = seq_names
                stats[gname] = _build_stats_for_group(f[gname], reactivity.shape[0])

        # Optional CIF visualization (per-group defattrs).
        defattr_files: list[str] = []
        chimerax_md = ""
        if cif_path is not None and ref_sequence is not None:
            defattr_dir = os.path.join(job_dir, "defattr")
            os.makedirs(defattr_dir, exist_ok=True)
            cif_workdir_path = os.path.join(workdir, os.path.basename(cif_path))
            shutil.copy(cif_path, cif_workdir_path)
            paths, md = _build_defattrs(cif_workdir_path, ref_sequence, results, defattr_dir)
            defattr_files = [os.path.basename(p) for p in paths]
            chimerax_md = md
            if defattr_files:
                log(f"Wrote {len(defattr_files)} defattr file(s).")
        elif cif_path is not None:
            log(
                "Skipping structure visualization: defattr generation requires "
                "a single-reference FASTA."
            )

        meta = {
            "job_id": job_id,
            "group_names": group_names,
            "sequence_names": sequence_names_per_group,
            "single_ref": single_ref,
            "has_combined_profile": has_combined,
            "ref_count": int(np.array(results[0].combined.reactivity).shape[0]),
            "is_pairwise": core_cfg.compute_pairwise,
            "stats": stats,
            "defattr_files": defattr_files,
            "chimerax_md": chimerax_md,
            "created_at": time.time(),
        }
        write_meta(job_dir, meta)
        write_log(job_dir, state.log_lines)

        log(f"\nDone. Generated profiles for {len(results)} group(s).")
        state.status = "done"

    except Exception as e:  # noqa: BLE001
        msg = str(e) or "Unexpected error. See log for details."
        log(f"Error: {msg}")
        log(traceback.format_exc())
        state.error = msg
        state.status = "error"
        # Still persist meta + log so the results page can show what happened.
        write_meta(job_dir, {
            "job_id": job_id,
            "status": "error",
            "error": msg,
            "created_at": time.time(),
        })
        write_log(job_dir, state.log_lines)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
