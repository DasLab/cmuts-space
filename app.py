#!/usr/bin/env python3
"""Gradio app for cmuts: chemical mutation profiling for RNA structure analysis."""

from __future__ import annotations

import glob
import os
import shutil
import subprocess
import tempfile

import gradio as gr


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
) -> tuple[str | None, list[str], str]:
    """Run the full cmuts pipeline: align -> core -> normalize."""
    if fasta_file is None or mod_fastq is None:
        return None, [], "Please upload a FASTA file and at least one modified FASTQ file."

    workdir = tempfile.mkdtemp(prefix="cmuts_")
    outdir = os.path.join(workdir, "outputs")
    os.makedirs(outdir)

    group_name = group_name.strip() or "profile"
    log_lines = []

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
        # Copy inputs to workdir
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
        fastq_files = sorted(glob.glob(os.path.join(fastq_dir, "*")))
        align_cmd = [
            "cmuts", "align",
            "--fasta", fasta_path,
            "--output", os.path.join(outdir, "alignments"),
            *fastq_files,
        ]
        if not run(align_cmd, cwd=outdir):
            return None, [], "\n".join(log_lines)

        # Step 2: Count mutations
        # Use relative paths so HDF5 groups match what normalize expects
        # (e.g. "alignments/bicine-2A3" not "/tmp/.../alignments/bicine-2A3")
        log("\n=== Step 2: Counting mutations ===")
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
            return None, [], "\n".join(log_lines)

        # Step 3: Normalize
        log("\n=== Step 3: Normalizing reactivities ===")
        profiles_path = os.path.join(outdir, "profiles.h5")

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
            return None, [], "\n".join(log_lines)

        # Collect figures
        fig_dir = os.path.join(outdir, "figures")
        figures = []
        if os.path.isdir(fig_dir):
            figures = sorted(glob.glob(os.path.join(fig_dir, "*.png")))

        log(f"\nDone. Generated {len(figures)} figures.")
        return profiles_path, figures, "\n".join(log_lines)

    except subprocess.TimeoutExpired:
        log("Pipeline timed out (10 minute limit).")
        return None, [], "\n".join(log_lines)
    except Exception as e:
        log(f"Error: {e}")
        return None, [], "\n".join(log_lines)


with gr.Blocks(title="cmuts — RNA Chemical Probing Analysis") as demo:
    gr.Markdown(
        """
        # cmuts — RNA Chemical Probing Analysis

        Upload a FASTA reference and FASTQ file(s) from a MaP-seq experiment
        to compute normalized reactivity profiles.

        **Pipeline:** `cmuts align` → `cmuts core` → `cmuts normalize`
        """
    )

    with gr.Row():
        with gr.Column():
            fasta_input = gr.File(label="Reference FASTA", file_types=[".fasta", ".fa"])
            mod_input = gr.File(label="Modified FASTQ (required)", file_types=[".fastq", ".fq", ".fastq.gz", ".fq.gz"])
            nomod_input = gr.File(label="Control FASTQ (optional)", file_types=[".fastq", ".fq", ".fastq.gz", ".fq.gz"])
            group_name = gr.Textbox(label="Group name", value="profile", placeholder="e.g. DMS, 2A3")

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

    with gr.Row():
        output_file = gr.File(label="Output HDF5")
        output_gallery = gr.Gallery(label="Figures", columns=2, height=400)

    output_log = gr.Textbox(label="Log", lines=15, max_lines=30)

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
        outputs=[output_file, output_gallery, output_log],
    )


if __name__ == "__main__":
    demo.launch(server_name="0.0.0.0", server_port=7860)
