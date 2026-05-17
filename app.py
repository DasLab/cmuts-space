"""FastAPI app for the cmuts web UI.

The HTTP layer; pipeline logic lives in ``pipeline.py``.
"""

from __future__ import annotations

import asyncio
import io
import json
import os
import shutil
import time
import uuid
import zipfile
from typing import Annotated

from fastapi import BackgroundTasks, FastAPI, Form, HTTPException, Request, UploadFile
from fastapi.responses import (
    FileResponse,
    HTMLResponse,
    JSONResponse,
    PlainTextResponse,
    RedirectResponse,
    StreamingResponse,
)
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from pipeline import (
    AlignConfig,
    CoreConfig,
    GroupInput,
    JobState,
    MAX_FASTQ_MB,
    MAX_GROUPS,
    NormConfig,
    PLOT_KEYS,
    RESULTS_DIR,
    RESULTS_TTL_HOURS,
    build_diff_plot,
    build_perref_plot,
    build_profile_plot,
    cleanup_old_results,
    file_size_mb,
    job_dir_for,
    read_meta,
    run_pipeline,
    fastq_safe_name,
    sanitize_group_name,
)


# --- App setup ---

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
TEMPLATES_DIR = os.path.join(BASE_DIR, "templates")
STATIC_DIR = os.path.join(BASE_DIR, "static")

app = FastAPI(title="cmuts")
app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")
templates = Jinja2Templates(directory=TEMPLATES_DIR)


def _asset_version() -> str:
    """Derive a cache-busting tag from static file mtimes so any deploy
    invalidates browser-cached CSS/JS."""
    try:
        mtimes = [
            os.path.getmtime(os.path.join(STATIC_DIR, f))
            for f in ("app.css", "results.js")
        ]
        return str(int(max(mtimes)))
    except OSError:
        return "0"


templates.env.globals["asset_version"] = _asset_version()


# In-memory state for jobs currently running in this process. Once a job
# completes, its results live on disk under RESULTS_DIR; this dict is
# pruned. A job that's missing from this dict but present on disk is
# treated as complete (status read from meta.json).
JOB_STATES: dict[str, JobState] = {}
_STATES_LOCK = asyncio.Lock()


@app.on_event("startup")
def _startup() -> None:
    cleanup_old_results()


# --- Helpers ---


def _save_upload(upload: UploadFile, dest_dir: str, filename: str) -> str:
    """Save an uploaded file under ``dest_dir`` with the given filename and
    return its absolute path."""
    os.makedirs(dest_dir, exist_ok=True)
    out = os.path.join(dest_dir, filename)
    with open(out, "wb") as f:
        shutil.copyfileobj(upload.file, f)
    return out


def _is_real_upload(upload: UploadFile | None) -> bool:
    if upload is None:
        return False
    name = (upload.filename or "").strip()
    return bool(name)


def _job_status(job_id: str) -> tuple[str, list[str], str | None]:
    """Return (status, log_lines, error) for a job. Reads from in-memory
    state if the job is running here; otherwise from disk."""
    state = JOB_STATES.get(job_id)
    if state is not None:
        return state.status, list(state.log_lines), state.error
    meta = read_meta(job_dir_for(job_id))
    if meta is None:
        return "missing", [], None
    log_path = os.path.join(job_dir_for(job_id), "log.txt")
    log_lines: list[str] = []
    if os.path.isfile(log_path):
        with open(log_path) as f:
            log_lines = f.read().splitlines()
    if meta.get("status") == "error":
        return "error", log_lines, meta.get("error")
    return "done", log_lines, None


# --- Routes: form + dynamic rows ---


@app.get("/", response_class=HTMLResponse)
def index(request: Request) -> HTMLResponse:
    return templates.TemplateResponse(
        request,
        "index.html",
        {
            "max_groups": MAX_GROUPS,
            "ttl_hours": RESULTS_TTL_HOURS,
            "max_fastq_mb": MAX_FASTQ_MB,
        },
    )


@app.get("/group-row", response_class=HTMLResponse)
def group_row(request: Request) -> HTMLResponse:
    """HTMX partial: returns one new empty group row."""
    return templates.TemplateResponse(
        request, "_group_row.html", {"initial": None},
    )


EXAMPLES_DIR = os.environ.get(
    "CMUTS_EXAMPLES_DIR", os.path.join(BASE_DIR, "examples"),
)


@app.post("/run-example/{name}")
def run_example(name: str, background_tasks: BackgroundTasks):
    """Submit a bundled example dataset (subdirectory under EXAMPLES_DIR)."""
    safe = os.path.basename(name)
    src_dir = os.path.join(EXAMPLES_DIR, safe)
    if not os.path.isdir(src_dir):
        raise HTTPException(404, f"Unknown example dataset: {safe}")

    fasta = None
    treated = None
    untreated = None
    for f in sorted(os.listdir(src_dir)):
        path = os.path.join(src_dir, f)
        lower = f.lower()
        if lower.endswith((".fasta", ".fa")):
            fasta = path
        elif "untreated" in lower or "nomod" in lower or "control" in lower:
            untreated = path
        elif lower.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
            treated = path
    if fasta is None or treated is None:
        raise HTTPException(500, f"Example dataset {safe} is missing required files.")

    job_id = uuid.uuid4().hex[:12]
    job_dir = job_dir_for(job_id)
    uploads_dir = os.path.join(job_dir, "uploads")
    os.makedirs(uploads_dir, exist_ok=True)
    fasta_path = os.path.join(uploads_dir, "ref.fasta")
    shutil.copy(fasta, fasta_path)
    mod_path = os.path.join(uploads_dir, "example__" + os.path.basename(treated))
    shutil.copy(treated, mod_path)
    nomod_path: str | None = None
    if untreated is not None:
        nomod_path = os.path.join(uploads_dir, "example__" + os.path.basename(untreated))
        shutil.copy(untreated, nomod_path)

    state = JobState(job_id=job_id)
    state.log(f"Submitted example job {job_id}.")
    JOB_STATES[job_id] = state

    groups = [GroupInput(name="example", mod_fastq=mod_path, nomod_fastq=nomod_path)]

    def _runner() -> None:
        try:
            run_pipeline(
                job_id=job_id, job_dir=job_dir,
                fasta_path=fasta_path, groups=groups,
                align_cfg=AlignConfig(),
                core_cfg=CoreConfig(),
                norm_cfg=NormConfig(),
                cif_path=None, state=state,
            )
        finally:
            asyncio.get_event_loop().call_later(60, JOB_STATES.pop, job_id, None)

    background_tasks.add_task(asyncio.to_thread, _runner)
    return RedirectResponse(f"/results/{job_id}", status_code=303)


# --- Routes: submit + status ---


@app.post("/run")
async def run(
    request: Request,
    background_tasks: BackgroundTasks,
    fasta: UploadFile,
    cif: UploadFile | None = None,
    group_name: list[str] = Form(default=[]),
    mod_fastq: list[UploadFile] = Form(default=[]),
    nomod_fastq: list[UploadFile] = Form(default=[]),
    # Alignment
    trim_5: str = Form(default=""),
    trim_3: str = Form(default=""),
    local_align: bool = Form(default=False),
    # Core
    min_mapq: int = Form(default=10),
    min_phred: int = Form(default=10),
    min_length: int = Form(default=2),
    max_length: int = Form(default=1024),
    no_mismatches: bool = Form(default=False),
    strand: str = Form(default="both"),
    compute_pairwise: bool = Form(default=False),
    sig: float = Form(default=0.05),
    # Norm
    norm_method: str = Form(default="ubr"),
    no_insertions: bool = Form(default=True),
    no_deletions: bool = Form(default=False),
    clip_below: str = Form(default=""),
    clip_above: str = Form(default=""),
    blank_5p: int = Form(default=0),
    blank_3p: int = Form(default=0),
    blank_cutoff: int = Form(default=10),
    norm_cutoff: int = Form(default=500),
    norm_percentile: int = Form(default=90),
):
    if not _is_real_upload(fasta):
        raise HTTPException(400, "A reference FASTA is required.")

    # Validate group inputs and stage files to the job dir.
    job_id = uuid.uuid4().hex[:12]
    job_dir = job_dir_for(job_id)
    os.makedirs(job_dir, exist_ok=True)
    uploads_dir = os.path.join(job_dir, "uploads")
    os.makedirs(uploads_dir, exist_ok=True)

    fasta_path = _save_upload(fasta, uploads_dir, "ref.fasta")
    cif_path: str | None = None
    if _is_real_upload(cif):
        cif_path = _save_upload(cif, uploads_dir, cif.filename or "ref.cif")

    # Pair the form list values. Browsers always send equal-length lists,
    # but be defensive.
    n_rows = max(len(group_name), len(mod_fastq), len(nomod_fastq))
    while len(group_name) < n_rows:
        group_name.append("")
    while len(mod_fastq) < n_rows:
        mod_fastq.append(UploadFile(filename="", file=io.BytesIO()))  # type: ignore[call-arg]
    while len(nomod_fastq) < n_rows:
        nomod_fastq.append(UploadFile(filename="", file=io.BytesIO()))  # type: ignore[call-arg]

    groups: list[GroupInput] = []
    seen_names: set[str] = set()
    for i in range(n_rows):
        mod = mod_fastq[i]
        if not _is_real_upload(mod):
            continue
        gn = sanitize_group_name(group_name[i]) or f"group_{i + 1}"
        # Disambiguate clashing names.
        base = gn
        k = 2
        while gn in seen_names:
            gn = f"{base}_{k}"
            k += 1
        seen_names.add(gn)

        mod_path = _save_upload(mod, uploads_dir, f"{fastq_safe_name(gn)}__{mod.filename}")
        if file_size_mb(mod_path) > MAX_FASTQ_MB:
            raise HTTPException(
                400,
                f"Modified FASTQ for group '{gn}' is "
                f"{file_size_mb(mod_path):.0f} MB; the limit is "
                f"{MAX_FASTQ_MB} MB.",
            )
        nomod_path: str | None = None
        if _is_real_upload(nomod_fastq[i]):
            nomod_path = _save_upload(
                nomod_fastq[i], uploads_dir,
                f"{fastq_safe_name(gn)}__{nomod_fastq[i].filename}",
            )
            if file_size_mb(nomod_path) > MAX_FASTQ_MB:
                raise HTTPException(
                    400,
                    f"Control FASTQ for group '{gn}' is "
                    f"{file_size_mb(nomod_path):.0f} MB; the limit is "
                    f"{MAX_FASTQ_MB} MB.",
                )
        groups.append(GroupInput(name=gn, mod_fastq=mod_path, nomod_fastq=nomod_path))

    if not groups:
        raise HTTPException(400, "At least one group with a Modified FASTQ is required.")

    align_cfg = AlignConfig(trim_5=trim_5, trim_3=trim_3, local_align=local_align)
    core_cfg = CoreConfig(
        min_mapq=min_mapq, min_phred=min_phred,
        min_length=min_length, max_length=max_length,
        no_insertions=no_insertions, no_mismatches=no_mismatches,
        strand=strand, compute_pairwise=compute_pairwise,
    )
    norm_cfg = NormConfig(
        norm_method=norm_method,
        no_insertions=no_insertions, no_deletions=no_deletions,
        clip_below=float(clip_below) if clip_below.strip() else None,
        clip_above=float(clip_above) if clip_above.strip() else None,
        blank_5p=blank_5p, blank_3p=blank_3p,
        blank_cutoff=blank_cutoff, norm_cutoff=norm_cutoff,
        norm_percentile=norm_percentile, sig=sig,
    )

    state = JobState(job_id=job_id)
    state.log(f"Submitted job {job_id} with {len(groups)} group(s).")
    JOB_STATES[job_id] = state

    def _runner() -> None:
        try:
            run_pipeline(
                job_id=job_id,
                job_dir=job_dir,
                fasta_path=fasta_path,
                groups=groups,
                align_cfg=align_cfg,
                core_cfg=core_cfg,
                norm_cfg=norm_cfg,
                cif_path=cif_path,
                state=state,
            )
        finally:
            # Keep the state briefly so the UI can read the final transition,
            # then drop it. The on-disk meta + log are authoritative after.
            import threading
            threading.Timer(60.0, JOB_STATES.pop, args=(job_id, None)).start()

    background_tasks.add_task(asyncio.to_thread, _runner)
    return RedirectResponse(f"/results/{job_id}", status_code=303)


@app.get("/results/{job_id}/status", response_class=JSONResponse)
def status(job_id: str) -> JSONResponse:
    status_, log_lines, error = _job_status(job_id)
    return JSONResponse({
        "status": status_,
        "log": "\n".join(log_lines),
        "error": error,
    })


# --- Routes: results page + plot data ---


@app.get("/results/{job_id}", response_class=HTMLResponse)
def results(request: Request, job_id: str) -> HTMLResponse:
    status_, log_lines, error = _job_status(job_id)
    if status_ == "missing":
        return templates.TemplateResponse(
            request,
            "results.html",
            {
                "job_id": job_id,
                "status": "missing",
                "ttl_hours": RESULTS_TTL_HOURS,
                "log": "",
                "error": None,
                "meta": None,
                "first_group": None,
            },
            status_code=404,
        )
    meta = read_meta(job_dir_for(job_id))
    first_group = None
    if meta and meta.get("group_names"):
        first_group = meta["group_names"][0]
    return templates.TemplateResponse(
        request,
        "results.html",
        {
            "job_id": job_id,
            "status": status_,
            "ttl_hours": RESULTS_TTL_HOURS,
            "log": "\n".join(log_lines),
            "error": error,
            "meta": meta,
            "first_group": first_group,
            "plot_keys": PLOT_KEYS,
        },
    )


def _read_plot_json(path: str) -> str | None:
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        return f.read()


@app.get("/results/{job_id}/plot/combined", response_class=PlainTextResponse)
def plot_combined(job_id: str) -> PlainTextResponse:
    path = os.path.join(job_dir_for(job_id), "combined_profile.json")
    body = _read_plot_json(path)
    if body is None:
        raise HTTPException(404, "Combined plot not available.")
    return PlainTextResponse(body, media_type="application/json")


@app.get("/results/{job_id}/plot/diff", response_class=PlainTextResponse)
def plot_diff(job_id: str, a: str, b: str) -> PlainTextResponse:
    body = build_diff_plot(job_dir_for(job_id), a, b)
    if body is None:
        raise HTTPException(404, "Diff plot not available for these groups.")
    return PlainTextResponse(body, media_type="application/json")


@app.get("/results/{job_id}/plot/{group}/{key}", response_class=PlainTextResponse)
def plot_group_key(
    job_id: str, group: str, key: str, seq: str = "0",
) -> PlainTextResponse:
    if key not in PLOT_KEYS:
        raise HTTPException(404, "Unknown plot key.")
    job_dir = job_dir_for(job_id)
    saved_path = os.path.join(job_dir, "groups", group, f"{key}.json")

    # "all" is a special profile-only mode that returns the pre-saved
    # heatmap-across-references plot. For every other plot the pre-saved
    # JSON is the first-reference view.
    if seq == "all":
        if key == "profile":
            body = _read_plot_json(saved_path)
            if body is not None:
                return PlainTextResponse(body, media_type="application/json")
        # Other tiles fall back to first-reference view in "all" mode.
        seq = "0"

    try:
        seq_idx = int(seq)
    except ValueError:
        raise HTTPException(400, f"Invalid seq value: {seq}")

    if key == "profile":
        # Always build single-reference profiles on demand so seq=0 shows
        # the first reference, not the all-references heatmap.
        body = build_profile_plot(job_dir, group, seq_idx)
    elif key in {"mi", "correlation", "pairwise_coverage"}:
        if seq_idx == 0:
            body = _read_plot_json(saved_path)
            if body is None:
                body = build_perref_plot(job_dir, group, key, seq_idx)
        else:
            body = build_perref_plot(job_dir, group, key, seq_idx)
    else:
        body = _read_plot_json(saved_path)

    if body is None:
        raise HTTPException(404, "Plot not available for this group/sequence.")
    return PlainTextResponse(body, media_type="application/json")


# --- Routes: downloads ---


@app.get("/results/{job_id}/download/h5")
def download_h5(job_id: str) -> FileResponse:
    path = os.path.join(job_dir_for(job_id), "profiles.h5")
    if not os.path.isfile(path):
        raise HTTPException(404, "HDF5 not available.")
    meta = read_meta(job_dir_for(job_id)) or {}
    group_names = meta.get("group_names") or [job_id]
    filename = "-".join(fastq_safe_name(g) for g in group_names) + "-profiles.h5"
    return FileResponse(path, media_type="application/x-hdf5", filename=filename)


@app.get("/results/{job_id}/download/csv")
def download_csv(job_id: str) -> FileResponse:
    path = os.path.join(job_dir_for(job_id), "profiles.csv")
    if not os.path.isfile(path):
        raise HTTPException(404, "CSV not available.")
    meta = read_meta(job_dir_for(job_id)) or {}
    group_names = meta.get("group_names") or [job_id]
    filename = "-".join(fastq_safe_name(g) for g in group_names) + "-profiles.csv"
    return FileResponse(path, media_type="text/csv", filename=filename)


@app.get("/results/{job_id}/download/defattr/{name}")
def download_defattr(job_id: str, name: str) -> FileResponse:
    safe = os.path.basename(name)  # prevent traversal
    path = os.path.join(job_dir_for(job_id), "defattr", safe)
    if not os.path.isfile(path):
        raise HTTPException(404, "Defattr not available.")
    return FileResponse(path, media_type="text/plain", filename=safe)


def _plot_json_to_png(json_path: str) -> bytes | None:
    """Render a saved Plotly JSON to PNG bytes. Returns None on failure."""
    try:
        import plotly.graph_objects as go
        with open(json_path) as f:
            fig = go.Figure(json.load(f))
        return fig.to_image(format="png", width=1000, height=600, scale=2)
    except Exception:
        return None


@app.get("/results/{job_id}/download/all")
def download_all(job_id: str) -> StreamingResponse:
    job_dir = job_dir_for(job_id)
    if not os.path.isdir(job_dir):
        raise HTTPException(404, "Job not found.")

    def _gen():
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED) as zf:
            for root, _dirs, files in os.walk(job_dir):
                rel = os.path.relpath(root, job_dir)
                # Skip raw uploads to keep the bundle small.
                if rel.startswith("uploads"):
                    continue
                for name in files:
                    full = os.path.join(root, name)
                    arc_rel = os.path.relpath(full, job_dir)
                    # Convert plot JSONs to PNG; include the underlying h5,
                    # csv, log, defattr, meta as-is.
                    if name.endswith(".json") and (
                        rel.startswith("groups")
                        or name == "combined_profile.json"
                    ):
                        png = _plot_json_to_png(full)
                        if png is not None:
                            zf.writestr(
                                os.path.join(job_id, arc_rel[:-5] + ".png"),
                                png,
                            )
                        continue
                    if name == "meta.json":
                        # Skip — internal app state, not useful to the user.
                        continue
                    zf.write(full, os.path.join(job_id, arc_rel))
        buf.seek(0)
        yield buf.read()

    headers = {"Content-Disposition": f'attachment; filename="{job_id}-results.zip"'}
    return StreamingResponse(_gen(), media_type="application/zip", headers=headers)


# --- Local dev entry point ---


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(
        "app:app",
        host=os.environ.get("HOST", "127.0.0.1"),
        port=int(os.environ.get("PORT", "7860")),
        reload=bool(os.environ.get("RELOAD", "1") == "1"),
    )
