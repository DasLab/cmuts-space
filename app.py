"""FastAPI app for the cmuts web UI.

The HTTP layer only: the option form comes from options.py, the pipeline
runs in pipeline.py, and the report is a proxied `cmuts plot` server managed
by report.py.
"""

from __future__ import annotations

import asyncio
import io
import json
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

import job_description
import options
import pipeline
import report
from pipeline import (
    CONDITION_ROLES,
    JOB_DESCRIPTION_FILE,
    MAX_CONDITIONS,
    MAX_UPLOAD_MB,
    RESULTS_TTL_HOURS,
    SETTINGS_FILE,
    UPLOADS_DIR,
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
    print(f"cmuts: results directory {pipeline.RESULTS_DIR}", flush=True)
    pipeline.start_cleaner()
    report.start_reaper()
    _install_sigterm_handler()
    yield


app = FastAPI(title="cmuts", lifespan=_lifespan,
              docs_url=None, redoc_url=None, openapi_url=None)
app.mount("/static", StaticFiles(directory=os.path.join(BASE_DIR, "static")), name="static")
app.mount("/examples", StaticFiles(directory=EXAMPLES_DIR), name="examples")
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

# The sources that a file reference in a job description can name.
UPLOADS_SOURCE = "uploads"
EXAMPLES_SOURCE = "examples"

# This route serves the job description and the staged inputs of a job. The
# form fetches the files of a job through it, as it fetches the files of an
# example.
JOB_FILES_ROUTE = "jobs"

# This query parameter of the form page names a job to load into the form.
EDIT_PARAMETER = "edit"

templates.env.globals["job_files_route"] = JOB_FILES_ROUTE
templates.env.globals["edit_parameter"] = EDIT_PARAMETER
templates.env.globals["job_description_file"] = JOB_DESCRIPTION_FILE

# A job id is twelve hexadecimal characters. A query string that holds
# anything else is refused before it reaches a path.
JOB_ID_PATTERN = re.compile(r"[0-9a-f]{12}")

UPLOAD_CHUNK_BYTES = 1 << 20

# A settings file holds one small JSON object; anything larger is refused
# before it is parsed.
SETTINGS_MAX_BYTES = 1 << 20


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
    # A job description can name one file part more than once.
    upload.file.seek(0)
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


def _uploaded_part(form: FormData, part: str) -> UploadFile:
    uploads = _real_uploads(form.getlist(part))
    if not uploads:
        raise HTTPException(400, f"The request has no file part named {part}.")
    return uploads[0]


def _example_path(relative: str) -> str:
    """Returns the path of one bundled example file. Refuses a path that
    leaves the examples directory."""
    root = os.path.realpath(EXAMPLES_DIR)
    path = os.path.realpath(os.path.join(root, relative))
    if not path.startswith(root + os.sep) or not os.path.isfile(path):
        raise HTTPException(400, f"There is no example file {relative}.")
    return path


def _copy_file(source: str, dest_dir: str) -> str:
    os.makedirs(dest_dir, exist_ok=True)
    return shutil.copy(source, os.path.join(dest_dir, os.path.basename(source)))


def _stage_file(form: FormData, reference: str, dest_dir: str) -> str:
    """Stages the file that one reference names into dest_dir. Returns the
    staged path."""
    source, _, name = reference.partition("/")
    if source == UPLOADS_SOURCE:
        return _save_upload(_uploaded_part(form, name), dest_dir)
    if source == EXAMPLES_SOURCE:
        return _copy_file(_example_path(name), dest_dir)
    raise HTTPException(400, f"The file reference {reference} has no known source.")


def _stage_role(form: FormData, references: list[str], dest_dir: str) -> list[str]:
    paths = [_stage_file(form, reference, dest_dir) for reference in references]
    if len(set(paths)) != len(paths):
        raise HTTPException(400, "Two files of one role in one condition "
                                 "share a filename.")
    return paths


def _stage_condition(form: FormData, index: int, condition: ConditionInput,
                     uploads_dir: str) -> ConditionInput:
    """Stages the files of one condition. Returns the condition with each
    reference replaced by its staged path."""
    files = {
        role: _stage_role(form, getattr(condition, role),
                          os.path.join(uploads_dir, f"cond-{index}", role))
        for role in CONDITION_ROLES
    }
    return ConditionInput(name=condition.name, **files)


def _job_file_reference(job_id: str, job_dir: str, path: str) -> str:
    """Returns the reference through which the form fetches one staged input
    of a job."""
    return f"{JOB_FILES_ROUTE}/{job_id}/{os.path.relpath(path, job_dir)}"


def _saved_condition(job_id: str, job_dir: str, condition: ConditionInput) -> dict:
    """Returns one staged condition in the format of a job description."""
    files = {role: [_job_file_reference(job_id, job_dir, path)
                    for path in getattr(condition, role)]
             for role in CONDITION_ROLES}
    return {"name": condition.name, **files}


def _saved_job_description(job_id: str, job_dir: str, fasta_path: str,
                           conditions: list[ConditionInput], provided: dict) -> dict:
    """Returns the job description that loads a job back into the form. It
    holds the staged inputs of the job and the options that the job gave."""
    return {
        "reference": _job_file_reference(job_id, job_dir, fasta_path),
        "conditions": [_saved_condition(job_id, job_dir, condition)
                       for condition in conditions],
        "options": provided,
    }


def _write_run_settings(job_dir: str, settings: dict) -> None:
    pipeline.write_settings(job_dir, options.settings_document(SPECS, settings))


def _submit_job(background_tasks: BackgroundTasks, job_id: str, job_dir: str,
                fasta_path: str, conditions: list[ConditionInput],
                extra: dict[str, list[str]]) -> JSONResponse:
    state = JobState(job_id=job_id)
    state.log(f"Submitted job {job_id} with {len(conditions)} condition(s).")
    JOB_STATES[job_id] = state

    def _runner() -> None:
        try:
            pipeline.run_pipeline(job_id, job_dir, fasta_path, conditions, extra, state)
        finally:
            threading.Timer(60.0, JOB_STATES.pop, args=(job_id, None)).start()

    background_tasks.add_task(asyncio.to_thread, _runner)
    return JSONResponse({"job_id": job_id, "url": f"/results/{job_id}"})


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


@app.get("/")
def index(request: Request, job: str | None = None) -> Response:
    """Serves the form, or redirects to the job that the query string names.
    The page embedding this app holds the current job in its own URL, so a
    reload of that page arrives here and returns to the job."""
    if job is not None and JOB_ID_PATTERN.fullmatch(job):
        return RedirectResponse(f"/results/{job}", status_code=303)
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


def _stage_run_inputs(form: FormData, description: job_description.JobDescription,
                      job_dir: str) -> tuple[str, list[ConditionInput]]:
    """Stages the reference and the files of every condition into the job
    directory."""
    uploads_dir = os.path.join(job_dir, UPLOADS_DIR)
    fasta_path = _stage_file(form, description.reference, uploads_dir)
    conditions = [_stage_condition(form, index, condition, uploads_dir)
                  for index, condition in enumerate(description.conditions)]
    return fasta_path, conditions


def _job_document(form: FormData) -> dict:
    """Returns the parsed job description that the "job" field holds."""
    raw = form.get("job")
    if not isinstance(raw, str):
        raise HTTPException(400, "The request has no job description.")
    try:
        return json.loads(raw)
    except json.JSONDecodeError:
        raise HTTPException(400, "The job description is not valid JSON.")


@app.post("/run", response_class=JSONResponse)
async def run(request: Request, background_tasks: BackgroundTasks) -> JSONResponse:
    """Submits one job. The request is a multipart form. Its "job" field holds
    the job description, and its file parts hold the uploads that the
    description names. The form and the bundled examples both submit here."""
    form = await request.form()
    document = _job_document(form)
    try:
        description = job_description.read_job_description(document)
        provided = options.document_settings(SPECS, document)
        settings = options.run_settings(SPECS, provided)
    except ValueError as error:
        raise HTTPException(400, str(error))
    extra = options.all_option_args(SPECS, settings)

    job_id = uuid.uuid4().hex[:12]
    job_dir = job_dir_for(job_id)
    try:
        fasta_path, conditions = _stage_run_inputs(form, description, job_dir)
    except HTTPException:
        shutil.rmtree(job_dir, ignore_errors=True)
        raise

    _write_run_settings(job_dir, settings)
    pipeline.write_job_description(job_dir, _saved_job_description(
        job_id, job_dir, fasta_path, conditions, provided))
    return _submit_job(background_tasks, job_id, job_dir, fasta_path, conditions, extra)


# --- Routes: settings ---


def _settings_upload(form: FormData) -> UploadFile:
    uploads = _real_uploads(form.getlist("settings"))
    if not uploads:
        raise HTTPException(400, "No settings file was chosen.")
    return uploads[0]


async def _uploaded_document(upload: UploadFile) -> dict:
    """The parsed settings file, refused where it is too large to be one or
    is not JSON at all."""
    raw = await upload.read(SETTINGS_MAX_BYTES + 1)
    if len(raw) > SETTINGS_MAX_BYTES:
        raise HTTPException(400, "The settings file is too large.")
    try:
        return json.loads(raw)
    except (json.JSONDecodeError, UnicodeDecodeError):
        raise HTTPException(400, "The settings file is not valid JSON.")


@app.post("/settings", response_class=JSONResponse)
async def check_settings(request: Request) -> JSONResponse:
    """Checks an uploaded settings file against the option dumps and answers
    with the form fields it sets, which the page writes into the form."""
    form = await request.form()
    document = await _uploaded_document(_settings_upload(form))
    try:
        fields = options.settings_fields(SPECS, document)
    except ValueError as error:
        raise HTTPException(400, str(error))
    return JSONResponse({"fields": fields})


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
            "settings_file": SETTINGS_FILE,
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
    return files + [SETTINGS_FILE]


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


# --- Routes: the inputs of a job ---


def _job_input_path(job_id: str, relative: str) -> str:
    """Returns the path of the job description or of one staged input of a
    job. Refuses every other path, so the outputs and the log stay behind
    their own routes."""
    if not JOB_ID_PATTERN.fullmatch(job_id):
        raise HTTPException(404, "Job not found.")
    root = os.path.realpath(job_dir_for(job_id))
    path = os.path.realpath(os.path.join(root, relative))
    staged = path.startswith(os.path.join(root, UPLOADS_DIR) + os.sep)
    described = path == os.path.join(root, JOB_DESCRIPTION_FILE)
    if not (staged or described) or not os.path.isfile(path):
        raise HTTPException(404, "No such file.")
    return path


@app.get(f"/{JOB_FILES_ROUTE}/{{job_id}}/{{relative:path}}")
def job_input(job_id: str, relative: str) -> FileResponse:
    """Serves the job description or one staged input of a job. The form
    loads them to run the job again."""
    return FileResponse(_job_input_path(job_id, relative))


# --- Local dev entry point ---


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(
        "app:app",
        host=os.environ.get("HOST", "127.0.0.1"),
        port=int(os.environ.get("PORT", "7860")),
        reload=os.environ.get("RELOAD", "1") == "1",
    )
