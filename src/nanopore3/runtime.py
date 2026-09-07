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


def _read(path: Path) -> str:
    try:
        return path.read_text(encoding="ascii").strip()
    except (OSError, UnicodeError):
        return ""


def _cgroup_directories(proc: Path = Path("/proc")) -> dict[str, list[Path]]:
    """Resolve this process's v1/v2 cgroups and their mounted ancestors."""
    memberships = []
    for line in _read(proc / "self/cgroup").splitlines():
        parts = line.split(":", 2)
        if len(parts) == 3:
            memberships.append((set(parts[1].split(",")) - {""}, Path(parts[2])))
    result: dict[str, list[Path]] = {}
    for line in _read(proc / "self/mountinfo").splitlines():
        fields = line.split()
        try:
            split = fields.index("-")
            kind, options = fields[split + 1], fields[split + 3]
            if kind not in {"cgroup", "cgroup2"}:
                continue

            # mountinfo uses octal escapes for spaces, tabs and backslashes.
            def unescape(value: str) -> str:
                import re

                return re.sub(r"\\([0-7]{3})", lambda m: chr(int(m[1], 8)), value)

            mount_root, mount = Path(unescape(fields[3])), Path(unescape(fields[4]))
            for controllers, member in memberships:
                names = (
                    {"cpu", "cpuset", "memory"}
                    if kind == "cgroup2"
                    else controllers & set(options.split(","))
                )
                if (kind == "cgroup2") != (not controllers) or not names:
                    continue
                try:
                    relative = member.relative_to(mount_root)
                except ValueError:
                    # A cgroup namespace may report membership relative to its
                    # root while mountinfo still names the host subtree.
                    relative = member.relative_to("/")
                if ".." in relative.parts:
                    continue
                current = mount / relative
                while True:
                    for name in names:
                        result.setdefault(name, []).append(current)
                    if current == mount:
                        break
                    current = current.parent
        except (ValueError, IndexError):
            continue
    return result


def _cpuset_count(value: str) -> int | None:
    try:
        ranges = []
        for part in value.split(","):
            endpoints = [int(n) for n in part.split("-")]
            if len(endpoints) not in (1, 2) or endpoints[0] < 0 or endpoints[-1] < endpoints[0]:
                return None
            ranges.append((endpoints[0], endpoints[-1]))
        total, end = 0, -1
        for start, stop in sorted(ranges):
            total += max(0, stop - max(start, end + 1) + 1)
            end = max(end, stop)
        return total or None
    except ValueError:
        return None


def effective_cpu_count(*, proc: Path = Path("/proc")) -> int:
    """Maximum whole-CPU budget allowed by host, affinity and cgroup limits."""
    limits = [max(1, os.cpu_count() or 1)]
    for probe in (
        getattr(os, "process_cpu_count", None),
        (lambda: len(os.sched_getaffinity(0))) if hasattr(os, "sched_getaffinity") else None,
    ):
        if probe:
            try:
                value = probe()
                if value:
                    limits.append(value)
            except (OSError, NotImplementedError):
                pass
    groups = _cgroup_directories(proc)
    for directory in groups.get("cpu", []):
        quota = _read(directory / "cpu.max").split()
        if not quota:
            quota = [_read(directory / "cpu.cfs_quota_us"), _read(directory / "cpu.cfs_period_us")]
        try:
            amount, period = map(int, quota)
            if amount > 0 and period > 0:
                limits.append(max(1, amount // period))
        except ValueError:
            pass
    for directory in groups.get("cpuset", []):
        value = _read(directory / "cpuset.cpus.effective") or _read(directory / "cpuset.cpus")
        count = _cpuset_count(value)
        if count:
            limits.append(count)
    return min(limits)


def group_memory_limit(configured_mb: int, *, proc: Path = Path("/proc")) -> int:
    """Retained-group estimate budget, reserving 3/4 of a cgroup for other work.

    This is admission control for Python groups, not an OS RSS hard limit.
    """
    limits = [configured_mb * 1024 * 1024]
    for directory in _cgroup_directories(proc).get("memory", []):
        value = _read(directory / "memory.max") or _read(directory / "memory.limit_in_bytes")
        try:
            amount = int(value)
            if amount > 0:
                limits.append(max(1, amount // 4))
        except ValueError:
            pass
    return min(limits)


@dataclass(frozen=True, slots=True)
class ResourcePlan:
    jobs: int
    threads_per_job: int
    available_cpus: int


def plan_resources(jobs: int, threads_per_job: int) -> ResourcePlan:
    """Bound nested parallelism to the detected CPU budget."""

    available = effective_cpu_count()
    if jobs < 1 or threads_per_job < 1:
        raise ValueError("jobs and threads_per_job must be positive")
    bounded_threads = min(threads_per_job, available)
    bounded_jobs = min(jobs, max(1, available // bounded_threads))
    return ResourcePlan(bounded_jobs, bounded_threads, available)


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
    from .config import ParallelSettings
    from .provenance import analysis_implementation

    defaults = ParallelSettings()
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
        "available_cpus": effective_cpu_count(),
        "parallel_defaults": asdict(defaults),
        "group_memory_limit_bytes": group_memory_limit(defaults.group_memory_mb),
        "analysis_implementation": analysis_implementation(),
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
