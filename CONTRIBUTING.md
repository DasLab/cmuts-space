# Environment and Runs

## Setting up

The app needs `cmuts` (with `cmuts-plot`) on `PATH`. Alignment also needs `minimap2`, `samtools`, and `vsearch`. The Python dependencies go in a virtual environment:

```sh
uv venv .venv
uv sync
```

The form and the bundled examples submit a job the same way. They post a job description to `/run`, and `job_description.py` describes its format. A bundled example is a directory under `examples/` that holds its read files and a `job.json` description. That description names each file as `examples/<example>/<file>` and holds every option the example runs with. The page loads that description into the form and uploads the files it names, so a user can change anything before running it. `scripts/smoke.sh` posts the description unchanged.

The example FASTQ files go through Git LFS (see `.gitattributes`), so a clone needs `git lfs` installed.

## Running locally

```sh
uv run python app.py
```

This starts the server at http://127.0.0.1:7860. Results are written under `/tmp/cmuts-space-results` by default.

## Docker

```sh
docker build -t cmuts-space .
docker run --rm -p 7860:7860 cmuts-space
```

The Dockerfile pins cmuts by commit through the `CMUTS_SHA` build argument. The `Pin space to this commit` workflow in [DasLab/cmuts](https://github.com/DasLab/cmuts) updates the pin on every push to its `main` branch, so do not pin it by hand.

# Deployment

The GitHub repository ([DasLab/cmuts-space](https://github.com/DasLab/cmuts-space)) is the source of truth. Push only to GitHub.

The `Sync HF Space` workflow deploys every commit on `main` to the [Hugging Face space](https://huggingface.co/spaces/daslab-stanford/cmuts). It first builds the Docker image on the runner and runs the bundled examples against it with `scripts/smoke.sh`. Only if that passes does it push to the space, wait for the build there, and run the examples once more against the deployment.

The sync workflow needs a repository secret named `HF_TOKEN` holding a Hugging Face token with write access to the space. The pin workflow in DasLab/cmuts needs a `SPACE_TOKEN` secret there. That is a fine-grained token which may write the contents of this repository and nothing else. The current token expires on September 12, 2027, and must be refreshed then.

The pin workflow pushes to this repository on every push to DasLab/cmuts, so always pull before editing. A stale checkout carries an out-of-date `CMUTS_SHA`.

The previous version of the space lives on the `v1` branch.
