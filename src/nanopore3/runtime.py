"""Cross-platform capability discovery and conservative resource planning."""

from __future__ import annotations

import importlib.metadata
import os
import platform
import shutil
import subprocess
import sys
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class ResourcePlan:
    jobs: int
    threads_per_job: int
    available_cpus: int


def plan_resources(jobs: int, threads_per_job: int) -> ResourcePlan:
    """Bound nested parallelism to the detected CPU budget."""

    available = max(1, os.cpu_count() or 1)
    if jobs < 1 or threads_per_job < 1:
        raise ValueError("jobs and threads_per_job must be positive")
    bounded_jobs = min(jobs, max(1, available // threads_per_job))
    return ResourcePlan(bounded_jobs, threads_per_job, available)


def _binary_version(name: str, args: Iterable[str]) -> dict[str, str | bool | None]:
    path = shutil.which(name)
    if path is None:
        return {"available": False, "path": None, "version": None}
    try:
        result = subprocess.run(
            [path, *args], capture_output=True, text=True, timeout=10, check=False, shell=False
        )
        text = (result.stdout or result.stderr).strip().splitlines()
        version = text[0] if text else f"exit {result.returncode}"
    except (OSError, subprocess.SubprocessError) as exc:
        version = f"probe failed: {exc}"
    return {"available": True, "path": str(Path(path).resolve()), "version": version}


def doctor_report() -> dict[str, object]:
    """Return machine-readable runtime and optional-backend capabilities."""

    packages: dict[str, str | None] = {}
    for package in ("nanopore3", "edlib", "PyYAML"):
        try:
            packages[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            packages[package] = None
    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cpu_count": max(1, os.cpu_count() or 1),
        "packages": packages,
        "portable_backend": True,
        "optional_binaries": {
            "mafft": _binary_version("mafft", ("--version",)),
            "spoa": _binary_version("spoa", ("--version",)),
            "minimap2": _binary_version("minimap2", ("--version",)),
        },
    }


def resource_plan_dict(plan: ResourcePlan) -> dict[str, int]:
    return asdict(plan)

