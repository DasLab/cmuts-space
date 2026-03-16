#!/usr/bin/env python3
"""Gradio app for cmuts: chemical mutation profiling for RNA structure analysis."""

from __future__ import annotations

import glob
import json
import os
import shutil
import subprocess
import tempfile
import time
import uuid

import gradio as gr
import h5py
import numpy as np
import plotly.graph_objects as go
from fastapi.responses import FileResponse, HTMLResponse

EXAMPLES_DIR = os.path.join(os.path.dirname(__file__), "examples")
MAX_FASTQ_MB = 500
RESULTS_TTL_HOURS = 48

# Use HF persistent storage if available, else fall back to /tmp
RESULTS_DIR = "/data/results" if os.path.isdir("/data") else "/tmp/cmuts_results"
os.makedirs(RESULTS_DIR, exist_ok=True)


def _file_size_mb(path: str) -> float:
    return os.path.getsize(path) / (1024 * 1024)


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


def save_results(
    h5_path: str,
    group_name: str,
    fig: go.Figure,
    stats_md: str,
    names: list[str],
) -> str:
    """Save pipeline results to persistent storage. Returns the job ID."""
    job_id = uuid.uuid4().hex[:12]
    job_dir = os.path.join(RESULTS_DIR, job_id)
    os.makedirs(job_dir)

    shutil.copy(h5_path, os.path.join(job_dir, "profiles.h5"))

    meta = {
        "group_name": group_name,
        "created_at": time.time(),
        "names": names,
        "stats_md": stats_md,
    }
    with open(os.path.join(job_dir, "meta.json"), "w") as f:
        json.dump(meta, f)

    with open(os.path.join(job_dir, "plot.json"), "w") as f:
        f.write(fig.to_json())

    return job_id


def _build_profile_plot(
    reactivity: np.ndarray,
    sequence: str | None,
    title: str,
) -> go.Figure:
    """Build an interactive Plotly bar chart for a single reactivity profile."""
    x = np.arange(1, len(reactivity) + 1)
    if sequence:
        hover = [
            f"{nt}{p}<br>Reactivity: {r:.4f}"
            for p, r, nt in zip(x, reactivity, sequence)
        ]
    else:
        hover = [f"Position {p}<br>Reactivity: {r:.4f}" for p, r in zip(x, reactivity)]

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=x,
        y=reactivity,
        hovertext=hover,
        hoverinfo="text",
        marker_color="indianred",
    ))
    fig.update_layout(
        title=title,
        xaxis_title="Position",
        yaxis_title="Reactivity",
        template="plotly_white",
        height=400,
        margin=dict(l=50, r=20, t=40, b=40),
    )
    return fig


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


def _build_stats_table(h5_path: str, group_name: str) -> str:
    """Build a markdown summary table from the output HDF5 file."""
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
        ("References", f"{n_refs:,}"),
        ("Reference length", f"{seq_len:,}"),
        ("Total reads", f"{total_reads:,}"),
        ("Mean reads per reference", f"{np.mean(reads):,.1f}"),
        ("Median reads per reference", f"{int(np.median(reads)):,}"),
    ]

    if valid.any():
        rows.extend([
            ("Mean reactivity", f"{np.mean(reactivity[valid]):.3f}"),
            ("Mean error", f"{np.mean(error[valid]):.3f}"),
            ("Mean SNR", f"{np.mean(snr):.2f}"),
            ("SNR > 1", f"{np.mean(snr > 1):.1%}"),
        ])

    dropout = float(np.mean(reads == 0))
    if dropout > 0:
        rows.append(("Dropout fraction", f"{dropout:.1%}"))

    md = "| Statistic | Value |\n|-----------|-------|\n"
    md += "\n".join(f"| {label} | {value} |" for label, value in rows)
    return md


def run_pipeline(
    fasta_file: str,
    mod_fastq: str,
    nomod_fastq: str | None,
    group_name: str,
    norm_method: str,
    no_insertions: bool,
    no_deletions: bool,
    clip_low: bool,
    clip_high: bool,
    # Alignment options
    trim_5: str,
    trim_3: str,
    local_align: bool,
    # Read filtering
    min_mapq: int,
    min_phred: int,
    min_length: int,
    max_length: int,
    no_mismatches: bool,
    strand: str,
    # Normalization
    blank_5p: int,
    blank_3p: int,
    blank_cutoff: int,
    norm_cutoff: int,
    norm_percentile: int,
):
    """Run the full cmuts pipeline: align -> core -> normalize.

    Yields (output_file, plot, sequence_dropdown_update, stats, result_url, log)
    so the log updates in real time and the interactive plot appears at the end.
    """
    empty_plot = go.Figure()
    empty_plot.update_layout(template="plotly_white", height=400)

    if fasta_file is None or mod_fastq is None:
        yield None, empty_plot, None, "", "", "Please upload a FASTA file and at least one modified FASTQ file."
        return

    for path, label in [(mod_fastq, "Modified FASTQ"), (nomod_fastq, "Control FASTQ")]:
        if path is not None and _file_size_mb(path) > MAX_FASTQ_MB:
            yield None, empty_plot, None, "", "", (
                f"{label} is {_file_size_mb(path):.0f} MB. "
                f"The free tier has limited RAM (16 GB); files over {MAX_FASTQ_MB} MB "
                f"may cause out-of-memory errors. Consider downsampling first."
            )
            return

    cleanup_old_results()

    workdir = tempfile.mkdtemp(prefix="cmuts_")
    outdir = os.path.join(workdir, "outputs")
    os.makedirs(outdir)

    group_name = (group_name or "").strip() or "profile"
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
            timeout=600,
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
        mod_name = mod_basename.split(".")[0]
        shutil.copy(mod_fastq, os.path.join(fastq_dir, mod_basename))

        nomod_name = None
        if nomod_fastq is not None:
            nomod_basename = os.path.basename(nomod_fastq)
            nomod_name = nomod_basename.split(".")[0]
            shutil.copy(nomod_fastq, os.path.join(fastq_dir, nomod_basename))

        # Step 1: Align
        log("=== Step 1: Aligning reads ===")
        yield None, empty_plot, None, "", "", "\n".join(log_lines)

        fastq_files = sorted(glob.glob(os.path.join(fastq_dir, "*")))
        align_cmd = [
            "cmuts", "align",
            "--fasta", fasta_path,
            "--output", os.path.join(outdir, "alignments"),
        ]
        if trim_5 and trim_5.strip():
            align_cmd.extend(["--trim-5", trim_5.strip()])
        if trim_3 and trim_3.strip():
            align_cmd.extend(["--trim-3", trim_3.strip()])
        if local_align:
            align_cmd.append("--local")
        align_cmd.extend(fastq_files)
        if not run(align_cmd, cwd=outdir):
            yield None, empty_plot, None, "", "", "\n".join(log_lines)
            return
        yield None, empty_plot, None, "", "", "\n".join(log_lines)

        # Step 2: Count mutations
        log("\n=== Step 2: Counting mutations ===")
        yield None, empty_plot, None, "", "", "\n".join(log_lines)

        bam_files = sorted(
            os.path.relpath(p, outdir)
            for p in glob.glob(os.path.join(outdir, "alignments", "*.bam"))
        )
        core_cmd = [
            "cmuts", "core",
            "-f", fasta_path,
            "-o", "counts.h5",
            "--min-mapq", str(min_mapq),
            "--min-phred", str(min_phred),
            "--min-length", str(min_length),
            "--max-length", str(max_length),
        ]
        if no_insertions:
            core_cmd.append("--no-insertions")
        if no_mismatches:
            core_cmd.append("--no-mismatches")
        if strand == "forward":
            core_cmd.append("--no-reverse")
        elif strand == "reverse":
            core_cmd.append("--only-reverse")
        core_cmd.extend(bam_files)
        if not run(core_cmd, cwd=outdir):
            yield None, empty_plot, None, "", "", "\n".join(log_lines)
            return
        yield None, empty_plot, None, "", "", "\n".join(log_lines)

        # Step 3: Normalize
        log("\n=== Step 3: Normalizing reactivities ===")
        yield None, empty_plot, None, "", "", "\n".join(log_lines)

        mod_group = f"alignments/{mod_name}"
        norm_cmd = [
            "cmuts", "normalize",
            "-o", "profiles.h5",
            "--mod", mod_group,
            "--fasta", fasta_path,
            "--group", group_name,
            "--norm", norm_method,
            "--blank-5p", str(blank_5p),
            "--blank-3p", str(blank_3p),
            "--blank-cutoff", str(blank_cutoff),
            "--norm-cutoff", str(norm_cutoff),
            "--norm-percentile", str(norm_percentile),
        ]
        if nomod_name:
            nomod_group = f"alignments/{nomod_name}"
            norm_cmd.extend(["--nomod", nomod_group])
        if no_insertions:
            norm_cmd.append("--no-insertions")
        if no_deletions:
            norm_cmd.append("--no-deletions")
        if clip_low:
            norm_cmd.append("--clip-low")
        if clip_high:
            norm_cmd.append("--clip-high")
        norm_cmd.append("counts.h5")

        if not run(norm_cmd, cwd=outdir):
            yield None, empty_plot, None, "", "", "\n".join(log_lines)
            return

        final_name = f"{group_name}-profiles.h5"
        final_path = os.path.join(outdir, final_name)
        os.rename(os.path.join(outdir, "profiles.h5"), final_path)

        reactivity, names = _read_profiles(final_path, group_name)
        fig = _build_profile_plot(reactivity[0], names[0], names[0])
        dropdown_update = gr.Dropdown(choices=names, value=names[0], visible=len(names) > 1)
        stats_md = _build_stats_table(final_path, group_name)

        # Save results for persistent access
        job_id = save_results(final_path, group_name, fig, stats_md, names)
        space_host = os.environ.get("SPACE_HOST", "")
        base = f"https://{space_host}" if space_host else ""
        result_url = f"{base}/results/{job_id}"

        log(f"\nDone. Generated {len(names)} profile(s).")
        log(f"Results available at: {result_url} (expires in {RESULTS_TTL_HOURS}h)")
        yield final_path, fig, dropdown_update, stats_md, result_url, "\n".join(log_lines)

    except subprocess.TimeoutExpired:
        log("Pipeline timed out (10 minute limit).")
        yield None, empty_plot, None, "", "", "\n".join(log_lines)
    except Exception as e:
        import traceback
        log(f"Error: {e}")
        log(traceback.format_exc())
        yield None, empty_plot, None, "", "", "\n".join(log_lines)


def select_profile(
    seq_name: str,
    output_file: str,
    group_name: str,
) -> go.Figure:
    """Switch the displayed profile when the user picks a different sequence."""
    if not output_file or not seq_name:
        return go.Figure()
    group_name = (group_name or "").strip() or "profile"
    reactivity, names = _read_profiles(output_file, group_name)
    try:
        idx = names.index(seq_name)
    except ValueError:
        idx = 0
    return _build_profile_plot(reactivity[idx], names[idx], names[idx])


def load_example():
    """Load bundled example files into the input fields."""
    return (
        os.path.join(EXAMPLES_DIR, "ref.fasta"),
        os.path.join(EXAMPLES_DIR, "treated.fastq.gz"),
        os.path.join(EXAMPLES_DIR, "untreated.fastq.gz"),
        "example",
    )


def load_saved_result(job_id: str):
    """Load a previously saved result by job ID."""
    job_id = (job_id or "").strip()
    if not job_id:
        return go.Figure(), None, "", ""

    job_dir = os.path.join(RESULTS_DIR, job_id)
    meta_path = os.path.join(job_dir, "meta.json")
    h5_path = os.path.join(job_dir, "profiles.h5")
    plot_path = os.path.join(job_dir, "plot.json")

    if not os.path.isdir(job_dir):
        return go.Figure(), None, "", "Result not found. It may have expired (results are kept for 48 hours)."

    with open(meta_path) as f:
        meta = json.load(f)

    with open(plot_path) as f:
        fig = go.Figure(json.load(f))

    names = meta.get("names", [])
    dropdown_update = gr.Dropdown(choices=names, value=names[0] if names else None, visible=len(names) > 1)

    return fig, dropdown_update, meta.get("stats_md", ""), ""


# ---------- Results page HTML template ----------

RESULTS_PAGE_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>cmuts results — {job_id}</title>
    <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
               max-width: 960px; margin: 2rem auto; padding: 0 1rem; color: #333; }}
        h1 {{ font-size: 1.5rem; }}
        h1 a {{ color: inherit; text-decoration: none; }}
        .meta {{ color: #666; margin-bottom: 1.5rem; }}
        table {{ border-collapse: collapse; margin: 1.5rem 0; }}
        th, td {{ border: 1px solid #ddd; padding: 0.5rem 1rem; text-align: left; }}
        th {{ background: #f5f5f5; }}
        .download {{ display: inline-block; margin: 1rem 0; padding: 0.5rem 1.5rem;
                     background: #4a90d9; color: white; text-decoration: none;
                     border-radius: 4px; }}
        .download:hover {{ background: #357abd; }}
        .expiry {{ color: #999; font-size: 0.85rem; margin-top: 2rem; }}
    </style>
</head>
<body>
    <h1><a href="/">cmuts</a> — Results</h1>
    <p class="meta">Job ID: <code>{job_id}</code></p>

    <div id="plot"></div>
    <script>
        var plotData = {plot_json};
        Plotly.newPlot('plot', plotData.data, plotData.layout, {{responsive: true}});
    </script>

    {stats_html}

    <a class="download" href="/results/{job_id}/download">Download HDF5 file</a>

    <p class="expiry">Results are stored for {ttl} hours and will be automatically deleted after that.</p>
</body>
</html>"""


def _md_table_to_html(md: str) -> str:
    """Convert a simple markdown table to HTML."""
    lines = [l.strip() for l in md.strip().split("\n") if l.strip() and not l.strip().startswith("|---")]
    if not lines:
        return ""
    html = "<table>\n"
    for i, line in enumerate(lines):
        cells = [c.strip() for c in line.strip("|").split("|")]
        tag = "th" if i == 0 else "td"
        html += "  <tr>" + "".join(f"<{tag}>{c}</{tag}>" for c in cells) + "</tr>\n"
    html += "</table>"
    return html


# ---------- Gradio UI ----------

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
        output_file = gr.File(label="Output HDF5")
        result_url = gr.Textbox(label="Result link (bookmark this — expires in 48h)", interactive=False)
        seq_dropdown = gr.Dropdown(label="Sequence", visible=False, interactive=True)
        output_plot = gr.Plot(label="Reactivity Profile")
        output_stats = gr.Markdown(label="Summary Statistics")
        with gr.Accordion("Log", open=False):
            output_log = gr.Textbox(label="Log", lines=15, max_lines=30, show_label=False)

        with gr.Accordion("Load previous results", open=False):
            with gr.Row():
                prev_job_id = gr.Textbox(label="Job ID", placeholder="e.g. a3f2b1c4d5e6", scale=3)
                load_btn = gr.Button("Load", variant="secondary", scale=1)
            load_status = gr.Textbox(label="Status", interactive=False, visible=False)

        example_btn.click(
            fn=load_example,
            outputs=[fasta_input, mod_input, nomod_input, group_name],
        )

        run_btn.click(
            fn=run_pipeline,
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
            outputs=[output_file, output_plot, seq_dropdown, output_stats, result_url, output_log],
        )

        seq_dropdown.change(
            fn=select_profile,
            inputs=[seq_dropdown, output_file, group_name],
            outputs=[output_plot],
        )

        load_btn.click(
            fn=load_saved_result,
            inputs=[prev_job_id],
            outputs=[output_plot, seq_dropdown, output_stats, load_status],
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
               with a bundled dataset (a ~615 nt RNA profiled with 2A3).
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
            you can bookmark or share. The link opens a standalone page with
            the interactive plot, summary statistics, and a download button
            for the HDF5 file. Results are stored for **48 hours** and
            automatically deleted after that.

            To reload previous results within the app, expand
            **Load previous results** on the Run tab and enter the job ID.

            ## Limits

            This server runs on Hugging Face Spaces with limited resources
            (16 GB RAM, CPU only). FASTQ files larger than 500 MB may cause
            out-of-memory errors. For larger datasets, install cmuts locally:

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


# ---------- FastAPI app with custom routes + mounted Gradio ----------

from fastapi import FastAPI

app = FastAPI()


@app.get("/results/{job_id}", response_class=HTMLResponse)
async def results_page(job_id: str):
    job_dir = os.path.join(RESULTS_DIR, job_id)
    if not os.path.isdir(job_dir):
        return HTMLResponse(
            "<h1>Result not found</h1><p>This result may have expired. "
            "Results are kept for 48 hours.</p>"
            '<p><a href="/">Return to cmuts</a></p>',
            status_code=404,
        )

    with open(os.path.join(job_dir, "meta.json")) as f:
        meta = json.load(f)
    with open(os.path.join(job_dir, "plot.json")) as f:
        plot_json = f.read()

    stats_html = _md_table_to_html(meta.get("stats_md", ""))

    html = RESULTS_PAGE_TEMPLATE.format(
        job_id=job_id,
        plot_json=plot_json,
        stats_html=stats_html,
        ttl=RESULTS_TTL_HOURS,
    )
    return HTMLResponse(html)


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
