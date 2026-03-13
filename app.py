#!/usr/bin/env python3
"""Gradio app for cmuts: chemical mutation profiling for RNA structure analysis."""

from __future__ import annotations

import glob
import os
import shutil
import subprocess
import tempfile
from collections.abc import Generator

import gradio as gr
import h5py
import numpy as np
import plotly.graph_objects as go

EXAMPLES_DIR = os.path.join(os.path.dirname(__file__), "examples")
MAX_FASTQ_MB = 500


def _file_size_mb(path: str) -> float:
    return os.path.getsize(path) / (1024 * 1024)


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
):
    """Run the full cmuts pipeline: align -> core -> normalize.

    Yields (output_file, plot, sequence_dropdown_update, log) so the log
    updates in real time and the interactive plot appears at the end.
    """
    empty_plot = go.Figure()
    empty_plot.update_layout(template="plotly_white", height=400)

    if fasta_file is None or mod_fastq is None:
        yield None, empty_plot, None, "", "Please upload a FASTA file and at least one modified FASTQ file."
        return

    for path, label in [(mod_fastq, "Modified FASTQ"), (nomod_fastq, "Control FASTQ")]:
        if path is not None and _file_size_mb(path) > MAX_FASTQ_MB:
            yield None, empty_plot, None, "", (
                f"{label} is {_file_size_mb(path):.0f} MB. "
                f"The free tier has limited RAM (16 GB); files over {MAX_FASTQ_MB} MB "
                f"may cause out-of-memory errors. Consider downsampling first."
            )
            return

    workdir = tempfile.mkdtemp(prefix="cmuts_")
    outdir = os.path.join(workdir, "outputs")
    os.makedirs(outdir)

    group_name = group_name.strip() or "profile"
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
        yield None, empty_plot, None, "", "\n".join(log_lines)

        fastq_files = sorted(glob.glob(os.path.join(fastq_dir, "*")))
        align_cmd = [
            "cmuts", "align",
            "--fasta", fasta_path,
            "--output", os.path.join(outdir, "alignments"),
            *fastq_files,
        ]
        if not run(align_cmd, cwd=outdir):
            yield None, empty_plot, None, "", "\n".join(log_lines)
            return
        yield None, empty_plot, None, "", "\n".join(log_lines)

        # Step 2: Count mutations
        log("\n=== Step 2: Counting mutations ===")
        yield None, empty_plot, None, "", "\n".join(log_lines)

        bam_files = sorted(
            os.path.relpath(p, outdir)
            for p in glob.glob(os.path.join(outdir, "alignments", "*.bam"))
        )
        core_cmd = [
            "cmuts", "core",
            "-f", fasta_path,
            "-o", "counts.h5",
        ]
        if no_insertions:
            core_cmd.append("--no-insertions")
        core_cmd.extend(bam_files)
        if not run(core_cmd, cwd=outdir):
            yield None, empty_plot, None, "", "\n".join(log_lines)
            return
        yield None, empty_plot, None, "", "\n".join(log_lines)

        # Step 3: Normalize
        log("\n=== Step 3: Normalizing reactivities ===")
        yield None, empty_plot, None, "", "\n".join(log_lines)

        mod_group = f"alignments/{mod_name}"
        norm_cmd = [
            "cmuts", "normalize",
            "-o", "profiles.h5",
            "--mod", mod_group,
            "--fasta", fasta_path,
            "--group", group_name,
            "--norm", norm_method,
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
            yield None, empty_plot, None, "", "\n".join(log_lines)
            return

        final_name = f"{group_name}-profiles.h5"
        final_path = os.path.join(outdir, final_name)
        os.rename(os.path.join(outdir, "profiles.h5"), final_path)

        reactivity, names = _read_profiles(final_path, group_name)
        fig = _build_profile_plot(reactivity[0], names[0], names[0])
        dropdown_update = gr.Dropdown(choices=names, value=names[0], visible=len(names) > 1)
        stats_md = _build_stats_table(final_path, group_name)

        log(f"\nDone. Generated {len(names)} profile(s).")
        yield final_path, fig, dropdown_update, stats_md, "\n".join(log_lines)

    except subprocess.TimeoutExpired:
        log("Pipeline timed out (10 minute limit).")
        yield None, empty_plot, None, "", "\n".join(log_lines)
    except Exception as e:
        log(f"Error: {e}")
        yield None, empty_plot, None, "", "\n".join(log_lines)


def select_profile(
    seq_name: str,
    output_file: str,
    group_name: str,
) -> go.Figure:
    """Switch the displayed profile when the user picks a different sequence."""
    if not output_file or not seq_name:
        return go.Figure()
    group_name = group_name.strip() or "profile"
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
        with gr.Row():
            with gr.Column():
                fasta_input = gr.File(label="Reference FASTA", file_types=[".fasta", ".fa"])
                mod_input = gr.File(label="Modified FASTQ (required)", file_types=[".fastq", ".fq", ".fastq.gz", ".fq.gz"])
                nomod_input = gr.File(label="Control FASTQ (optional)", file_types=[".fastq", ".fq", ".fastq.gz", ".fq.gz"])
                group_name = gr.Textbox(label="Group name", value="experiment", placeholder="e.g. DMS, 2A3")
                example_btn = gr.Button("Load example data", variant="secondary", size="sm")

            with gr.Column():
                norm_method = gr.Radio(
                    choices=["ubr", "outlier", "raw"],
                    value="ubr",
                    label="Normalization method",
                )
                no_insertions = gr.Checkbox(label="Exclude insertions", value=True)
                no_deletions = gr.Checkbox(label="Exclude deletions", value=False)
                clip_low = gr.Checkbox(label="Clip negative reactivities", value=False)
                clip_high = gr.Checkbox(label="Clip reactivities above 1", value=False)

        run_btn = gr.Button("Run Pipeline", variant="primary")

        output_file = gr.File(label="Output HDF5")
        seq_dropdown = gr.Dropdown(label="Sequence", visible=False, interactive=True)
        output_plot = gr.Plot(label="Reactivity Profile")
        output_stats = gr.Markdown(label="Summary Statistics")
        output_log = gr.Textbox(label="Log", lines=15, max_lines=30)

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
            ],
            outputs=[output_file, output_plot, seq_dropdown, output_stats, output_log],
        )

        seq_dropdown.change(
            fn=select_profile,
            inputs=[seq_dropdown, output_file, group_name],
            outputs=[output_plot],
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

            All uploaded data is processed in ephemeral temporary directories
            and deleted after the pipeline completes. No data is stored
            persistently, and no user accounts or tracking are used.
            """
        )


if __name__ == "__main__":
    demo.launch(server_name="0.0.0.0", server_port=7860)
