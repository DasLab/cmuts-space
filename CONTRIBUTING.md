# Environment and Runs

## Setting up

The app needs `cmuts` (with `cmuts-align` and `cmuts-plot`) on `PATH`, plus `minimap2`, `samtools`, and `fastp` for alignment, and a Python environment holding the dependencies:

```sh
uv venv .venv
uv sync
```

The example FASTQ files go through Git LFS (see `.gitattributes`), so a clone needs `git lfs` installed.

## Running locally

With a sibling cmuts checkout built at `../cmuts`:

```sh
PATH=$PWD/../cmuts/build/release:$PATH uv run python app.py
```

Then open http://127.0.0.1:7860 and run the bundled examples. Results land under `/tmp/cmuts-space-results` by default. The tunable limits are the `CMUTS_*` environment variables read at the top of `pipeline.py` and `report.py`.

## Docker

```sh
docker build -t cmuts-space .
docker run --rm -p 7860:7860 cmuts-space
```

The Dockerfile pins cmuts by commit through the `CMUTS_SHA` build argument. The `Pin space to this commit` workflow in [DasLab/cmuts](https://github.com/DasLab/cmuts) updates the pin on every push to its `main` branch, so do not pin it by hand.

# Deployment

The GitHub repository ([DasLab/cmuts-space](https://github.com/DasLab/cmuts-space)) is the source of truth. Push only to GitHub.

The `Sync HF Space` workflow force-pushes every commit on `main` to the [Hugging Face space](https://huggingface.co/spaces/daslab-stanford/cmuts), which builds the Dockerfile and deploys it. Do not push to the space directly; a direct push is overwritten by the next sync. The workflow needs a repository secret named `HF_TOKEN` holding a Hugging Face token with write access to the space.

The pin workflow in DasLab/cmuts needs a `SPACE_TOKEN` secret there: a fine-grained token with contents write access to this repository only. The current token expires on September 12, 2027, and must be refreshed then.

The pin workflow pushes to this repository on every push to DasLab/cmuts, so always pull before editing; a stale checkout carries an out-of-date `CMUTS_SHA`.

The previous version of the space lives on the `v1` branch.
