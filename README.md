---
title: cmuts
emoji: 🧬
colorFrom: blue
colorTo: purple
sdk: docker
pinned: false
license: mit
---

# cmuts web UI

FastAPI + Jinja + HTMX frontend for the
[cmuts](https://github.com/hmblair/cmuts) RNA chemical-probing pipeline.

## Local development

### Option A — Docker (full pipeline)

The Dockerfile builds `cmuts` from source along with its system
dependencies (bowtie2, samtools, libhts). This is the only way to run the
full pipeline locally on macOS.

```bash
docker build -t cmuts-space .
docker run --rm -p 7860:7860 cmuts-space
# Open http://localhost:7860
```

### Option B — Python venv (UI only; pipeline returns an error)

For working on the UI without rebuilding the cmuts toolchain:

```bash
uv sync
uv run python app.py
# Open http://localhost:7860
```

The form, group rows, results page, and error states all work.
Submitting a job will fail cleanly with a "cmuts is not installed"
message.

### Option C — venv with cmuts linked from a sibling checkout

If you have the `cmuts` repo built locally at `../cmuts` (sibling to
this directory), plus `bowtie2` and `samtools` on `PATH`, the pipeline
runs without Docker:

```bash
brew install bowtie2 samtools hdf5 htslib autoconf automake libtool libomp
# Build cmuts in ../cmuts (see github.com/hmblair/cmuts) — `./configure`
uv sync --extra cmuts
uv run python app.py
```

`uv sync` strips packages not declared in pyproject. The `cmuts` extra
declares the sibling path so re-running sync keeps it installed.

## Layout

- `app.py` — FastAPI routes, request handling
- `pipeline.py` — subprocess orchestration, HDF5/plot generation
- `templates/` — Jinja2 HTML
- `static/` — CSS + JS (Plotly via CDN)
- `examples/` — bundled FASTA/FASTQ for the "Run with example data" button
