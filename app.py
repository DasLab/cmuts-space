"""FastAPI app for the cmuts web UI.

The HTTP layer only: the option form comes from options.py, the pipeline
runs in pipeline.py, and the report is a proxied `cmuts plot` server managed
by report.py.
"""

from __future__ import annotations

import asyncio
import io
import os
import re
import shutil
import signal
import threading
import uuid
import zipfile
from contextlib import asynccontextmanager

import httpx
from fastapi import BackgroundTasks, FastAPI, HTTPException, Request
from fastapi.responses import (
    FileResponse,
    HTMLResponse,
    JSONResponse,
    RedirectResponse,
    Response,
    StreamingResponse,
)
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from starlette.datastructures import FormData, UploadFile

import options
import pipeline
import report
from pipeline import (
    MAX_CONDITIONS,
    MAX_UPLOAD_MB,
    RESULTS_TTL_HOURS,
    ConditionInput,
    JobState,
    job_dir_for,
    read_meta,
)

if not options.cmuts_available():
    raise SystemExit(
        "cmuts is not on PATH. For local runs, prepend a build directory:\n"
        "    PATH=$PWD/../cmuts/build/release:$PATH uv run python app.py"
    )

# The option dumps, read once against the installed binary so the form can
# never disagree with it.
SPECS, REQUIRED_OPTIONS, FORM_SECTIONS = options.load_form_spec()

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
EXAMPLES_DIR = os.environ.get("CMUTS_EXAMPLES_DIR", os.path.join(BASE_DIR, "examples"))

def _install_sigterm_handler() -> None:
    """Marks every running job interrupted when the container is replaced,
    then hands the signal on to the server's own handler. Registered from
    the lifespan, which runs after the server installs its handler."""
    previous = signal.getsignal(signal.SIGTERM)

    def handle(signum, frame):
        for state in list(JOB_STATES.values()):
            pipeline.mark_interrupted(state)
        if callable(previous):
            previous(signum, frame)
        else:
            signal.signal(signal.SIGTERM, signal.SIG_DFL)
            os.kill(os.getpid(), signal.SIGTERM)

    signal.signal(signal.SIGTERM, handle)


@asynccontextmanager
async def _lifespan(app: FastAPI):
    pipeline.cleanup_old_results()
    report.start_reaper()
    _install_sigterm_handler()
    yield


app = FastAPI(title="cmuts", lifespan=_lifespan,
              docs_url=None, redoc_url=None, openapi_url=None)
app.mount("/static", StaticFiles(directory=os.path.join(BASE_DIR, "static")), name="static")
templates = Jinja2Templates(directory=os.path.join(BASE_DIR, "templates"))


def _asset_version() -> str:
    """A cache-busting tag from the static file mtimes, so a deploy
    invalidates browser-cached CSS and JS."""
    static_dir = os.path.join(BASE_DIR, "static")
    mtimes = [entry.stat().st_mtime for entry in os.scandir(static_dir)]
    return str(int(max(mtimes)))


templates.env.globals["asset_version"] = _asset_version()

# In-memory state for jobs running in this process. A job missing from this
# dict but present on disk is complete; its meta and log are authoritative.
JOB_STATES: dict[str, JobState] = {}

PROXY_CLIENT = httpx.AsyncClient(timeout=60)

CONDITION_ROLES = ("treated", "untreated", "denatured")

UPLOAD_CHUNK_BYTES = 1 << 20


# --- Helpers ---


def _real_uploads(values: list) -> list[UploadFile]:
    return [v for v in values
            if isinstance(v, UploadFile) and (v.filename or "").strip()]


def _save_upload(upload: UploadFile, dest_dir: str) -> str:
    """Saves one upload in chunks, refusing it once it passes the size
    limit rather than after it has fully landed on disk."""
    os.makedirs(dest_dir, exist_ok=True)
    out = os.path.join(dest_dir, os.path.basename(upload.filename))
    limit = MAX_UPLOAD_MB * 1024 * 1024
    written = 0
    with open(out, "wb") as f:
        while chunk := upload.file.read(UPLOAD_CHUNK_BYTES):
            written += len(chunk)
            if written > limit:
                raise HTTPException(
                    400,
                    f"{upload.filename} exceeds the {MAX_UPLOAD_MB} MB limit.",
                )
            f.write(chunk)
    return out


def _condition_indices(form: FormData) -> list[int]:
    indices = set()
    for key in form.keys():
        match = re.fullmatch(r"cond-(\d+)-treated", key)
        if match:
            indices.add(int(match.group(1)))
    return sorted(indices)


def _stage_condition(form: FormData, index: int, uploads_dir: str) -> ConditionInput | None:
    """Stages one row's uploads; None where the row has no treated reads."""
    files: dict[str, list[str]] = {}
    for role in CONDITION_ROLES:
        uploads = _real_uploads(form.getlist(f"cond-{index}-{role}"))
        if len(uploads) > 2:
            raise HTTPException(400, f"At most two {role} files per condition "
                                     "(a read file and its mate).")
        names = [os.path.basename(u.filename) for u in uploads]
        if len(set(names)) != len(names):
            raise HTTPException(400, f"The {role} files of one condition "
                                     "share a filename.")
        dest = os.path.join(uploads_dir, f"cond-{index}", role)
        files[role] = [_save_upload(u, dest) for u in uploads]
    if not files["treated"]:
        return None
    name = (form.get(f"cond-{index}-name") or "").strip() or f"condition {index + 1}"
    return ConditionInput(name=name, treated=files["treated"],
                         untreated=files["untreated"], denatured=files["denatured"])


def _submit_job(background_tasks: BackgroundTasks, job_id: str, job_dir: str,
                fasta_path: str, conditions: list[ConditionInput],
                extra: dict[str, list[str]]) -> RedirectResponse:
    state = JobState(job_id=job_id)
    state.log(f"Submitted job {job_id} with {len(conditions)} condition(s).")
    JOB_STATES[job_id] = state

    def _runner() -> None:
        try:
            pipeline.run_pipeline(job_id, job_dir, fasta_path, conditions, extra, state)
        finally:
            threading.Timer(60.0, JOB_STATES.pop, args=(job_id, None)).start()

    background_tasks.add_task(asyncio.to_thread, _runner)
    return RedirectResponse(f"/results/{job_id}", status_code=303)


def _job_status(job_id: str) -> tuple[str, str, str | None]:
    """(status, log, error) from memory for a running job, from disk after."""
    state = JOB_STATES.get(job_id)
    if state is not None:
        return state.status, "\n".join(state.log_lines), state.error
    meta = read_meta(job_dir_for(job_id))
    if meta is None:
        return "missing", "", None
    log_path = os.path.join(job_dir_for(job_id), "log.txt")
    log = ""
    if os.path.isfile(log_path):
        with open(log_path) as f:
            log = f.read()
    if meta.get("status") == "error":
        return "error", log, meta.get("error")
    return "done", log, None


# --- Routes: form ---


@app.get("/", response_class=HTMLResponse)
def index(request: Request) -> HTMLResponse:
    return templates.TemplateResponse(request, "index.html", {
        "required_options": REQUIRED_OPTIONS,
        "sections": FORM_SECTIONS,
        "ttl_hours": RESULTS_TTL_HOURS,
        "max_upload_mb": MAX_UPLOAD_MB,
        "max_conditions": MAX_CONDITIONS,
    })


@app.get("/condition-row", response_class=HTMLResponse)
def condition_row(request: Request, index: int) -> HTMLResponse:
    return templates.TemplateResponse(request, "_condition_row.html", {"index": index})


# --- Routes: submit ---


def _stage_run_inputs(form: FormData, job_dir: str) -> tuple[str, list[ConditionInput]]:
    """Stages the FASTA and every condition row into the job directory."""
    fasta = _real_uploads(form.getlist("fasta"))
    if not fasta:
        raise HTTPException(400, "A reference FASTA is required.")

    uploads_dir = os.path.join(job_dir, "uploads")
    fasta_path = _save_upload(fasta[0], uploads_dir)

    conditions = [
        staged
        for index in _condition_indices(form)
        if (staged := _stage_condition(form, index, uploads_dir)) is not None
    ]
    if not conditions:
        raise HTTPException(400, "At least one condition with treated reads is required.")
    if len(conditions) > MAX_CONDITIONS:
        raise HTTPException(400, f"At most {MAX_CONDITIONS} conditions per run.")
    return fasta_path, conditions


@app.post("/run")
async def run(request: Request, background_tasks: BackgroundTasks):
    form = await request.form()

    try:
        extra = options.all_option_args(SPECS, form)
    except ValueError as error:
        raise HTTPException(400, str(error))

    job_id = uuid.uuid4().hex[:12]
    job_dir = job_dir_for(job_id)
    try:
        fasta_path, conditions = _stage_run_inputs(form, job_dir)
    except HTTPException:
        shutil.rmtree(job_dir, ignore_errors=True)
        raise

    return _submit_job(background_tasks, job_id, job_dir, fasta_path, conditions, extra)


@app.post("/run-example/{name}")
def run_example(name: str, background_tasks: BackgroundTasks):
    """Submits a bundled example dataset (a subdirectory under examples/)."""
    src_dir = os.path.join(EXAMPLES_DIR, os.path.basename(name))
    if not os.path.isdir(src_dir):
        raise HTTPException(404, f"Unknown example dataset: {name}")

    fasta, treated, untreated = None, None, None
    for f in sorted(os.listdir(src_dir)):
        path = os.path.join(src_dir, f)
        if f.endswith((".fasta", ".fa")):
            fasta = path
        elif "untreated" in f:
            untreated = path
        elif f.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
            treated = path
    if fasta is None or treated is None:
        raise HTTPException(500, f"Example dataset {name} is missing required files.")

    job_id = uuid.uuid4().hex[:12]
    job_dir = job_dir_for(job_id)
    uploads_dir = os.path.join(job_dir, "uploads")
    os.makedirs(uploads_dir, exist_ok=True)

    def stage(src: str) -> str:
        dest = os.path.join(uploads_dir, os.path.basename(src))
        shutil.copy(src, dest)
        return dest

    conditions = [ConditionInput(
        name="example",
        treated=[stage(treated)],
        untreated=[stage(untreated)] if untreated else [],
    )]
    extra = options.all_option_args(SPECS, FormData())
    return _submit_job(background_tasks, job_id, job_dir, stage(fasta), conditions, extra)


# --- Routes: results ---


@app.get("/results/{job_id}/status", response_class=JSONResponse)
def status(job_id: str) -> JSONResponse:
    status_, log, error = _job_status(job_id)
    return JSONResponse({"status": status_, "log": log, "error": error})


@app.get("/results/{job_id}", response_class=HTMLResponse)
def results(request: Request, job_id: str) -> HTMLResponse:
    status_, log, error = _job_status(job_id)
    return templates.TemplateResponse(
        request, "results.html",
        {
            "job_id": job_id,
            "status": status_,
            "ttl_hours": RESULTS_TTL_HOURS,
            "log": log,
            "error": error,
            "meta": read_meta(job_dir_for(job_id)),
        },
        status_code=404 if status_ == "missing" else 200,
    )


# --- Routes: report proxy ---


@app.get("/results/{job_id}/report")
def report_root(job_id: str) -> RedirectResponse:
    return RedirectResponse(f"/results/{job_id}/report/")


@app.get("/results/{job_id}/report/{path:path}")
async def report_proxy(job_id: str, path: str, request: Request) -> Response:
    try:
        port = await asyncio.to_thread(report.port_for, job_id)
    except FileNotFoundError:
        raise HTTPException(404, "Report not available.")
    except RuntimeError as error:
        raise HTTPException(502, str(error))

    url = f"http://127.0.0.1:{port}/{path}"
    upstream = await PROXY_CLIENT.get(url, params=request.query_params)
    return Response(
        content=upstream.content,
        status_code=upstream.status_code,
        media_type=upstream.headers.get("content-type"),
    )


# --- Routes: downloads ---


def _output_files(job_id: str) -> list[str]:
    meta = read_meta(job_dir_for(job_id)) or {}
    files = []
    for condition in meta.get("conditions", []):
        files.extend([f"{condition['tag']}.h5", f"{condition['tag']}.csv"])
    return files


@app.get("/results/{job_id}/download/all")
def download_all(job_id: str) -> StreamingResponse:
    job_dir = job_dir_for(job_id)
    names = [n for n in _output_files(job_id) + ["log.txt"]
             if os.path.isfile(os.path.join(job_dir, n))]
    if not names:
        raise HTTPException(404, "Job not found.")

    def _gen():
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED) as zf:
            for name in names:
                zf.write(os.path.join(job_dir, name), os.path.join(job_id, name))
        buf.seek(0)
        yield buf.read()

    headers = {"Content-Disposition": f'attachment; filename="{job_id}-results.zip"'}
    return StreamingResponse(_gen(), media_type="application/zip", headers=headers)


@app.get("/results/{job_id}/download/{filename}")
def download(job_id: str, filename: str) -> FileResponse:
    if filename not in _output_files(job_id):
        raise HTTPException(404, "No such file.")
    path = os.path.join(job_dir_for(job_id), filename)
    if not os.path.isfile(path):
        raise HTTPException(404, "No such file.")
    media = "text/csv" if filename.endswith(".csv") else "application/x-hdf5"
    return FileResponse(path, media_type=media, filename=filename)


# --- Local dev entry point ---


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(
        "app:app",
        host=os.environ.get("HOST", "127.0.0.1"),
        port=int(os.environ.get("PORT", "7860")),
        reload=os.environ.get("RELOAD", "1") == "1",
    )
