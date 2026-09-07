"""Subprocess jobs for MCP; no workflow algorithms or protocol dependencies here."""

from __future__ import annotations

import heapq
import json
import os
import re
import signal
import subprocess
import sys
import threading
from pathlib import Path
from typing import Any
from uuid import uuid4

import psutil

from .provenance import atomic_write_json, utc_now

_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}\Z")
_RESERVED = {"CON", "PRN", "AUX", "NUL"} | {
    f"{prefix}{n}" for prefix in ("COM", "LPT") for n in range(1, 10)
}


def portable_name(value: str) -> str:
    """Accept one short path component, including on native Windows."""
    if not _NAME.fullmatch(value) or value.upper() in _RESERVED:
        raise ValueError("Use 1–64 letters, digits, underscores or hyphens; no device names")
    return value


class Workspace:
    """Bound tool-addressable paths; trusted YAML may reference external inputs."""

    def __init__(self, root: Path):
        self.root = root.expanduser().resolve(strict=True)
        if not self.root.is_dir():
            raise ValueError("workspace must be an existing directory")

    def path(self, value: str | Path, *, exists: bool = False) -> Path:
        path = Path(value).expanduser()
        if not path.is_absolute():
            path = self.root / path
        path = path.resolve(strict=exists)
        if not path.is_relative_to(self.root):
            raise ValueError("Tool paths must stay inside the configured workspace")
        relative = path.relative_to(self.root)
        if any(part in {".git", ".venv", ".codex", ".agents"} for part in relative.parts):
            raise ValueError("Repository metadata and environment paths are not tool targets")
        return path

    def relative(self, path: Path) -> str:
        return path.relative_to(self.root).as_posix()


def tail(path: Path, limit: int = 16_384) -> str:
    if not path.is_file():
        return ""
    with path.open("rb") as handle:
        handle.seek(max(0, path.stat().st_size - limit))
        return handle.read(limit).decode("utf-8", errors="replace")


class JobManager:
    """One bounded job pool per server, with durable metadata and separate logs.

    Completed jobs remain readable after a restart. Unowned running jobs are
    reported as unknown, never reattached or signalled using a stale PID.
    """

    def __init__(self, workspace: Workspace, max_jobs: int = 1):
        if not 1 <= max_jobs <= 8:
            raise ValueError("max_jobs must be between 1 and 8")
        self.workspace = workspace
        self.max_jobs = max_jobs
        self.root = workspace.path(".nanopore3-mcp/jobs")
        self.root.mkdir(parents=True, exist_ok=True)
        self._lock = threading.RLock()
        self._active: dict[str, subprocess.Popen] = {}
        self._records: dict[str, dict[str, Any]] = {}
        self._threads: dict[str, threading.Thread] = {}
        self._closed = False

    def _directory(self, job_id: str) -> Path:
        return self.workspace.path(self.root / portable_name(job_id))

    def _save(self, record: dict[str, Any]) -> None:
        atomic_write_json(self._directory(record["job_id"]) / "job.json", record)

    def start(self, arguments: list[str], *, run_dir: Path | None = None) -> dict[str, Any]:
        with self._lock:
            if self._closed:
                raise ValueError("Server is shutting down")
            if len(self._active) >= self.max_jobs:
                raise ValueError("Job capacity reached; poll or cancel an existing job first")
            if run_dir is not None and any(
                record["run_dir"] is not None
                and self.workspace.path(record["run_dir"]) == run_dir
                for record in self._records.values()
            ):
                raise ValueError("Another active job already owns this run directory")
            job_id = uuid4().hex
            directory = self._directory(job_id)
            directory.mkdir()
            record = {
                "job_id": job_id, "operation": arguments[0], "state": "starting",
                "created_utc": utc_now(), "completed_utc": None, "returncode": None,
                "arguments": arguments,
                "run_dir": self.workspace.relative(run_dir) if run_dir else None,
            }
            self._save(record)
            try:
                with (directory / "stdout.log").open("wb") as stdout, (
                    directory / "stderr.log"
                ).open("wb") as stderr:
                    process = subprocess.Popen(
                        [sys.executable, "-m", "nanopore3", *arguments],
                        cwd=self.workspace.root, stdin=subprocess.DEVNULL,
                        stdout=stdout, stderr=stderr, shell=False,
                        env={
                            **os.environ, "PYTHONUNBUFFERED": "1", "PYTHONUTF8": "1",
                            "MPLCONFIGDIR": os.environ.get(
                                "MPLCONFIGDIR",
                                str(self.workspace.path(".nanopore3-mcp/matplotlib")),
                            ),
                        },
                        start_new_session=os.name != "nt",
                    )
            except OSError as exc:
                record.update(state="failed", error=str(exc), completed_utc=utc_now())
                self._save(record)
                raise
            record.update(state="running", pid=process.pid)
            self._records[job_id] = record
            self._active[job_id] = process
            self._save(record)
            thread = threading.Thread(target=self._wait, args=(job_id, process), daemon=True)
            self._threads[job_id] = thread
            thread.start()
            return dict(record)

    def _wait(self, job_id: str, process: subprocess.Popen) -> None:
        code = process.wait()
        with self._lock:
            record = self._records[job_id]
            record.update(
                state="cancelled" if record["state"] == "cancelling" else (
                    "succeeded" if code == 0 else "failed"
                ),
                returncode=code, completed_utc=utc_now(),
            )
            self._save(record)
            self._active.pop(job_id, None)
            self._records.pop(job_id, None)
            self._threads.pop(job_id, None)

    def status(self, job_id: str) -> dict[str, Any]:
        directory = self._directory(job_id)
        with self._lock:
            record = json.loads((directory / "job.json").read_text(encoding="utf-8"))
            if record["state"] in {"starting", "running", "cancelling"}:
                if job_id not in self._active:
                    record.update(
                        state="unknown",
                        warning="Job belongs to another server or a previous session. "
                        "Check the original process before resuming; it is not managed here.",
                    )
        record["stdout_tail"] = tail(directory / "stdout.log")
        record["stderr_tail"] = tail(directory / "stderr.log")
        if record["operation"] == "validate" and record["state"] == "succeeded":
            output = directory / "stdout.log"
            if output.stat().st_size <= 262_144:
                record["result"] = json.loads(output.read_text(encoding="utf-8"))
        return record

    def list_jobs(self, limit: int) -> dict[str, Any]:
        if not 1 <= limit <= 100:
            raise ValueError("limit must be between 1 and 100")
        paths = heapq.nlargest(
            limit + 1,
            (p for p in self.root.iterdir() if p.is_dir() and not p.is_symlink()),
            key=lambda p: p.stat().st_mtime_ns,
        )
        entries = []
        for path in paths[:limit]:
            record = self.status(path.name)
            entries.append({key: record[key] for key in (
                "job_id", "operation", "state", "created_utc", "run_dir",
            )})
        return {"jobs": entries, "truncated": len(paths) > limit}

    def cancel(self, job_id: str) -> dict[str, Any]:
        portable_name(job_id)
        with self._lock:
            process = self._active.get(job_id)
            if process is None or process.poll() is not None:
                return self.status(job_id)
            self._records[job_id]["state"] = "cancelling"
            self._save(self._records[job_id])
            # Hold the lock until descendants are gone so a second job cannot
            # enter while cancelled process workers still consume resources.
            if os.name != "nt":
                # A dedicated session covers even workers orphaned as the
                # parent exits, and needs no system-wide process enumeration.
                for sig in (signal.SIGTERM, signal.SIGKILL):
                    try:
                        os.killpg(process.pid, sig)
                    except ProcessLookupError:
                        break
                    if sig == signal.SIGTERM:
                        try:
                            process.wait(timeout=2)
                        except subprocess.TimeoutExpired:
                            pass
                return self.status(job_id)
            try:
                parent = psutil.Process(process.pid)
                descendants = parent.children(recursive=True)
                # Stop the parent first to prevent it scheduling further work.
                for item in [parent, *descendants]:
                    try:
                        item.terminate()
                    except psutil.NoSuchProcess:
                        pass
                _, alive = psutil.wait_procs([parent, *descendants], timeout=2)
                for item in alive:
                    try:
                        item.kill()
                    except psutil.NoSuchProcess:
                        pass
            except psutil.NoSuchProcess:
                pass
        return self.status(job_id)

    def close(self) -> None:
        with self._lock:
            self._closed = True
            active = list(self._active)
            threads = list(self._threads.values())
        for job_id in active:
            self.cancel(job_id)
        for thread in threads:
            thread.join(timeout=5)
