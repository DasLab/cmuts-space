"""Manages one live `cmuts plot` server per job.

The report page is the cmuts plot server itself: the app spawns one over a
job's final HDF5 files on the first request, proxies requests to it, and
kills it after an idle period. No report code lives in the space.
"""

from __future__ import annotations

import os
import re
import subprocess
import threading
import time
from dataclasses import dataclass

import pipeline

PLOT_IDLE_SEC = int(os.environ.get("CMUTS_PLOT_IDLE_SEC", "900"))
PLOT_START_TIMEOUT_SEC = 30

_URL_PATTERN = re.compile(r"http://127\.0\.0\.1:(\d+)/")


@dataclass
class PlotServer:
    process: subprocess.Popen
    port: int
    last_used: float


_SERVERS: dict[str, PlotServer] = {}
_LOCK = threading.Lock()


def plot_command(job_dir: str, meta: dict) -> list[str]:
    cmd = ["cmuts", "plot", "--host", "127.0.0.1", "--port", "0"]
    for condition in meta["conditions"]:
        cmd.append(os.path.join(job_dir, f"{condition['tag']}.h5"))
        cmd.extend(["--label", condition["name"]])
    return cmd


def spawn(job_dir: str, meta: dict) -> PlotServer:
    """Starts cmuts plot over the job's outputs and reads the port it bound
    from the URL it prints. A timer kills a process that prints no URL in
    time, so the blocking reads cannot hang past the deadline."""
    process = subprocess.Popen(
        plot_command(job_dir, meta),
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    killer = threading.Timer(PLOT_START_TIMEOUT_SEC, process.kill)
    killer.start()
    lines: list[str] = []
    try:
        for line in process.stdout:
            lines.append(line)
            match = _URL_PATTERN.search(line)
            if match:
                return PlotServer(process, int(match.group(1)), time.time())
    finally:
        killer.cancel()
    process.kill()
    detail = "".join(lines).strip() or "no output"
    raise RuntimeError(f"cmuts plot did not start: {detail}")


def port_for(job_id: str) -> int:
    """The proxy target for one job's report, spawning the plot server if
    none is running."""
    with _LOCK:
        server = _SERVERS.get(job_id)
        if server is not None and server.process.poll() is None:
            server.last_used = time.time()
            return server.port
        job_dir = pipeline.job_dir_for(job_id)
        meta = pipeline.read_meta(job_dir)
        if meta is None or "conditions" not in meta:
            raise FileNotFoundError(f"no results for job {job_id}")
        server = spawn(job_dir, meta)
        _SERVERS[job_id] = server
        return server.port


def reap_idle() -> None:
    """Kills plot servers idle past PLOT_IDLE_SEC. Called from a background
    thread."""
    cutoff = time.time() - PLOT_IDLE_SEC
    with _LOCK:
        idle = [job_id for job_id, s in _SERVERS.items()
                if s.last_used < cutoff or s.process.poll() is not None]
        for job_id in idle:
            _SERVERS.pop(job_id).process.kill()


def start_reaper() -> None:
    def loop() -> None:
        while True:
            time.sleep(60)
            reap_idle()

    threading.Thread(target=loop, daemon=True).start()
