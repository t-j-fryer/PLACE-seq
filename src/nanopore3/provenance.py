"""Content hashing, provenance records, and immutable stage publication.

A completed stage is a small transaction: work happens in a unique sibling
directory, all outputs are checksummed, ``manifest.json`` and ``_SUCCESS`` are
written last, and the directory is atomically renamed into place.  Resume never
trusts file presence alone; it validates the fingerprint and artifact hashes.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import re
import subprocess
import sys
from collections.abc import Mapping
from dataclasses import asdict, dataclass, is_dataclass
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any
from uuid import uuid4

from .errors import Nanopore3Error

MANIFEST_NAME = "manifest.json"
SUCCESS_NAME = "_SUCCESS"
STAGE_SCHEMA_VERSION = 1
_SAFE_STAGE_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


class ProvenanceError(Nanopore3Error, RuntimeError):
    """Base error for an invalid or unsafe provenance operation."""


class StageExistsError(ProvenanceError):
    """A completed stage exists but cannot be safely reused."""


class StageValidationError(ProvenanceError):
    """A completed stage is missing, modified, or fingerprint-incompatible."""


def utc_now() -> str:
    """Return an RFC 3339 UTC timestamp without local-time ambiguity."""

    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def sha256_bytes(data: bytes) -> str:
    """Return the lowercase SHA-256 hexadecimal digest of bytes."""

    return hashlib.sha256(data).hexdigest()


def sha256_file(path: str | Path, *, chunk_size: int = 1024 * 1024) -> str:
    """Stream a regular file into SHA-256 without loading it into memory."""

    if chunk_size < 1:
        raise ValueError("chunk_size must be >= 1")
    source = Path(path).expanduser().resolve(strict=True)
    if not source.is_file():
        raise OSError(f"Not a regular file: {source}")
    digest = hashlib.sha256()
    with source.open("rb") as handle:
        for block in iter(lambda: handle.read(chunk_size), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_ready(value: Any) -> Any:
    if is_dataclass(value) and not isinstance(value, type):
        return _json_ready(asdict(value))
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Enum):
        return _json_ready(value.value)
    if isinstance(value, datetime):
        if value.tzinfo is None:
            raise ValueError("Naive datetimes are not canonical; attach a timezone")
        return value.astimezone(timezone.utc).isoformat().replace("+00:00", "Z")
    if isinstance(value, Mapping):
        result: dict[str, Any] = {}
        for key, item in value.items():
            if not isinstance(key, str):
                raise TypeError("Canonical JSON mappings must use string keys")
            result[key] = _json_ready(item)
        return result
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    if isinstance(value, (set, frozenset)):
        converted = [_json_ready(item) for item in value]
        return sorted(converted, key=lambda item: canonical_json(item))
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    raise TypeError(f"Cannot encode {type(value).__name__} as canonical JSON")


def canonical_json(value: Any) -> str:
    """Serialize supported Python values into stable, compact UTF-8 JSON text."""

    return json.dumps(
        _json_ready(value),
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    )


def canonical_digest(value: Any) -> str:
    """Return SHA-256 of :func:`canonical_json` for a structured value."""

    return sha256_bytes(canonical_json(value).encode("utf-8"))


def compute_stage_fingerprint(
    stage: str,
    *,
    pipeline_version: str,
    parameters: Mapping[str, Any] | None = None,
    input_digests: Mapping[str, str] | None = None,
    reference_digest: str | None = None,
    backend_versions: Mapping[str, str] | None = None,
    contract_version: int = 1,
) -> str:
    """Create a content fingerprint controlling whether a stage may be resumed."""

    return canonical_digest(
        {
            "stage": stage,
            "contract_version": contract_version,
            "pipeline_version": pipeline_version,
            "analysis_implementation": analysis_implementation(),
            "parameters": dict(parameters or {}),
            "input_digests": dict(input_digests or {}),
            "reference_digest": reference_digest,
            "backend_versions": dict(backend_versions or {}),
        }
    )


def analysis_implementation(package_root: Path | None = None) -> dict[str, str]:
    """Hash installed source and presets, including uncommitted/local edits.

    Paths are relative and source newlines normalized, so a wheel and checkout
    have the same identity across platforms. No Git executable is required.
    """
    root = package_root or Path(__file__).resolve().parent
    sources = sorted(set(root.rglob("*.py")) | set(root.glob("presets/*.yaml")))
    if not sources:
        raise StageValidationError(f"cannot establish analysis implementation: {root}")
    return {
        "scheme": "nanopore3-source-v1",
        "sha256": canonical_digest(
            {
                path.relative_to(root).as_posix(): sha256_bytes(
                    path.read_bytes().replace(b"\r\n", b"\n")
                )
                for path in sources
            }
        ),
    }


def require_current_implementation(metadata: Mapping[str, Any]) -> None:
    """Refuse ambiguous legacy or changed-code reuse before touching stages."""
    recorded = metadata.get("analysis_implementation")
    current = analysis_implementation()
    if current != _IMPORTED_IMPLEMENTATION:
        raise StageValidationError(
            "analysis source changed since this process started; restart Python/the MCP server "
            "before starting or reusing a run"
        )
    if recorded != current:
        reason = (
            "missing implementation identity"
            if recorded is None
            else "analysis implementation differs"
        )
        raise StageValidationError(
            f"resume/rerun refused: {reason}; start a new run with a new run ID. "
            "Existing results remain readable."
        )


_IMPORTED_IMPLEMENTATION = analysis_implementation()


def atomic_write_bytes(path: str | Path, data: bytes) -> Path:
    """Atomically replace a file using a flushed temporary sibling."""

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.tmp-{uuid4().hex}")
    try:
        with temporary.open("xb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    except BaseException:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
        raise
    return destination


def atomic_write_text(path: str | Path, text: str) -> Path:
    """Atomically write UTF-8 text with no implicit newline transformation."""

    return atomic_write_bytes(path, text.encode("utf-8"))


def atomic_write_json(path: str | Path, value: Any) -> Path:
    """Atomically write canonical JSON followed by one newline."""

    return atomic_write_text(path, canonical_json(value) + "\n")


@dataclass(frozen=True, slots=True)
class ArtifactManifest:
    """Integrity and interpretation metadata for one stage artifact."""

    path: str
    sha256: str
    size_bytes: int
    media_type: str
    record_count: int | None = None
    semantic_sha256: str | None = None

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> ArtifactManifest:
        allowed = {
            "path",
            "sha256",
            "size_bytes",
            "media_type",
            "record_count",
            "semantic_sha256",
        }
        unknown = set(value) - allowed
        if unknown:
            raise StageValidationError(
                "Unknown artifact manifest key(s): " + ", ".join(sorted(unknown))
            )
        try:
            return cls(
                path=str(value["path"]),
                sha256=str(value["sha256"]),
                size_bytes=int(value["size_bytes"]),
                media_type=str(value["media_type"]),
                record_count=(
                    None
                    if value.get("record_count") is None
                    else int(value["record_count"])
                ),
                semantic_sha256=(
                    None
                    if value.get("semantic_sha256") is None
                    else str(value["semantic_sha256"])
                ),
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise StageValidationError(f"Malformed artifact manifest: {exc}") from exc


@dataclass(frozen=True, slots=True)
class StageManifest:
    """Immutable declaration of a successfully published stage."""

    schema_version: int
    stage: str
    fingerprint: str
    pipeline_version: str
    created_utc: str
    completed_utc: str
    parameters: Mapping[str, Any]
    input_digests: Mapping[str, str]
    backend_versions: Mapping[str, str]
    runtime: Mapping[str, Any]
    artifacts: tuple[ArtifactManifest, ...]

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> StageManifest:
        allowed = {
            "schema_version",
            "stage",
            "fingerprint",
            "pipeline_version",
            "created_utc",
            "completed_utc",
            "parameters",
            "input_digests",
            "backend_versions",
            "runtime",
            "artifacts",
        }
        unknown = set(value) - allowed
        missing = allowed - set(value)
        if unknown or missing:
            details = []
            if unknown:
                details.append("unknown=" + ",".join(sorted(unknown)))
            if missing:
                details.append("missing=" + ",".join(sorted(missing)))
            raise StageValidationError("Invalid stage manifest keys: " + "; ".join(details))
        try:
            artifacts_value = value["artifacts"]
            if not isinstance(artifacts_value, list):
                raise TypeError("artifacts must be a list")
            manifest = cls(
                schema_version=int(value["schema_version"]),
                stage=str(value["stage"]),
                fingerprint=str(value["fingerprint"]),
                pipeline_version=str(value["pipeline_version"]),
                created_utc=str(value["created_utc"]),
                completed_utc=str(value["completed_utc"]),
                parameters=dict(value["parameters"]),
                input_digests={
                    str(key): str(item) for key, item in value["input_digests"].items()
                },
                backend_versions={
                    str(key): str(item) for key, item in value["backend_versions"].items()
                },
                runtime=dict(value["runtime"]),
                artifacts=tuple(
                    ArtifactManifest.from_mapping(item) for item in artifacts_value
                ),
            )
        except (KeyError, TypeError, ValueError, AttributeError) as exc:
            raise StageValidationError(f"Malformed stage manifest: {exc}") from exc
        if manifest.schema_version != STAGE_SCHEMA_VERSION:
            raise StageValidationError(
                f"Unsupported stage manifest schema {manifest.schema_version}"
            )
        return manifest


def runtime_provenance() -> dict[str, Any]:
    """Capture a compact, stdlib-only runtime snapshot for a stage manifest."""

    return {
        "python": platform.python_version(),
        "implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "executable": sys.executable,
        "cpu_count": os.cpu_count(),
        "argv": list(sys.argv),
        "analysis_implementation": analysis_implementation(),
    }


def git_provenance(path: str | Path) -> dict[str, Any]:
    """Return commit and dirty-state evidence without requiring GitPython."""

    directory = Path(path).expanduser().resolve(strict=False)
    try:
        top = subprocess.run(
            ["git", "-C", str(directory), "rev-parse", "--show-toplevel"],
            capture_output=True,
            text=True,
            timeout=10,
            check=True,
            shell=False,
        ).stdout.strip()
        commit = subprocess.run(
            ["git", "-C", top, "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=10,
            check=True,
            shell=False,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "-C", top, "status", "--porcelain"],
            capture_output=True,
            text=True,
            timeout=10,
            check=True,
            shell=False,
        ).stdout
    except (OSError, subprocess.SubprocessError):
        return {"available": False, "root": None, "commit": None, "dirty": None}
    return {
        "available": True,
        "root": top,
        "commit": commit,
        "dirty": bool(status),
    }


def _safe_artifact_path(stage_dir: Path, relative: str) -> Path:
    relative_path = Path(relative)
    if relative_path.is_absolute() or ".." in relative_path.parts or relative_path == Path("."):
        raise StageValidationError(f"Unsafe artifact path in manifest: {relative!r}")
    candidate = stage_dir / relative_path
    if candidate.is_symlink():
        raise StageValidationError(f"Stage artifacts may not be symlinks: {relative}")
    try:
        candidate.resolve(strict=False).relative_to(stage_dir.resolve(strict=True))
    except ValueError as exc:
        raise StageValidationError(f"Artifact escapes stage directory: {relative}") from exc
    return candidate


def _read_manifest(path: Path) -> StageManifest:
    try:
        with path.open("r", encoding="utf-8") as handle:
            value = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise StageValidationError(f"Cannot read {path}: {exc}") from exc
    if not isinstance(value, Mapping):
        raise StageValidationError(f"Stage manifest is not an object: {path}")
    return StageManifest.from_mapping(value)


def validate_stage_directory(
    stage_dir: str | Path,
    *,
    expected_fingerprint: str | None = None,
    expected_stage: str | None = None,
    verify_checksums: bool = True,
) -> StageManifest:
    """Validate a completed stage before resume and return its manifest."""

    directory = Path(stage_dir).expanduser().resolve(strict=True)
    if not directory.is_dir() or directory.is_symlink():
        raise StageValidationError(f"Stage path is not a regular directory: {directory}")
    manifest_path = directory / MANIFEST_NAME
    success_path = directory / SUCCESS_NAME
    if not manifest_path.is_file() or not success_path.is_file():
        raise StageValidationError(
            f"Stage lacks {MANIFEST_NAME} and/or {SUCCESS_NAME}: {directory}"
        )
    manifest = _read_manifest(manifest_path)
    try:
        success_fingerprint = success_path.read_text(encoding="ascii").strip()
    except (OSError, UnicodeDecodeError) as exc:
        raise StageValidationError(f"Cannot read {success_path}: {exc}") from exc
    if success_fingerprint != manifest.fingerprint:
        raise StageValidationError("_SUCCESS fingerprint differs from manifest")
    if expected_fingerprint is not None and manifest.fingerprint != expected_fingerprint:
        raise StageValidationError(
            "Completed stage fingerprint does not match the requested inputs/configuration"
        )
    if expected_stage is not None and manifest.stage != expected_stage:
        raise StageValidationError(
            f"Completed stage is {manifest.stage!r}, expected {expected_stage!r}"
        )
    declared_paths: set[str] = set()
    for artifact in manifest.artifacts:
        if artifact.path in declared_paths:
            raise StageValidationError(f"Duplicate artifact path: {artifact.path}")
        declared_paths.add(artifact.path)
        candidate = _safe_artifact_path(directory, artifact.path)
        if not candidate.is_file():
            raise StageValidationError(f"Missing stage artifact: {artifact.path}")
        if candidate.stat().st_size != artifact.size_bytes:
            raise StageValidationError(f"Artifact size changed: {artifact.path}")
        if verify_checksums and sha256_file(candidate) != artifact.sha256:
            raise StageValidationError(f"Artifact checksum changed: {artifact.path}")
    actual_paths = {
        item.relative_to(directory).as_posix()
        for item in directory.rglob("*")
        if item.is_file() and item.name not in {MANIFEST_NAME, SUCCESS_NAME}
    }
    if actual_paths != declared_paths:
        extra = sorted(actual_paths - declared_paths)
        missing = sorted(declared_paths - actual_paths)
        raise StageValidationError(
            f"Stage file set changed; untracked={extra!r}, missing={missing!r}"
        )
    return manifest


def _guess_media_type(path: Path) -> str:
    name = path.name.lower()
    if name.endswith((".fastq.gz", ".fq.gz")):
        return "application/gzip; profile=fastq"
    if name.endswith((".fasta.gz", ".fa.gz", ".fna.gz")):
        return "application/gzip; profile=fasta"
    if name.endswith((".fastq", ".fq")):
        return "text/x-fastq"
    if name.endswith((".fasta", ".fa", ".fna")):
        return "text/x-fasta"
    if name.endswith(".json"):
        return "application/json"
    if name.endswith(".csv"):
        return "text/csv"
    if name.endswith(".parquet"):
        return "application/vnd.apache.parquet"
    return "application/octet-stream"


class StageDirectory:
    """Context manager for immutable, atomic and safely resumable stage output.

    Example::

        with StageDirectory(run_dir, "02_plate_demux", fingerprint) as stage:
            if not stage.reused:
                stage.output_path("calls.csv").write_text(...)

    On normal exit, all files below the temporary directory are checksummed and
    published together.  If a valid stage already exists and ``resume`` is true,
    the context points at it with ``reused=True`` and performs no mutation.
    """

    def __init__(
        self,
        run_directory: str | Path,
        stage: str,
        fingerprint: str,
        *,
        pipeline_version: str = "unknown",
        parameters: Mapping[str, Any] | None = None,
        input_digests: Mapping[str, str] | None = None,
        backend_versions: Mapping[str, str] | None = None,
        resume: bool = True,
        verify_on_resume: bool = True,
        keep_failed: bool = True,
    ) -> None:
        if not _SAFE_STAGE_NAME.fullmatch(stage):
            raise ValueError(
                "stage must contain only letters, digits, '.', '_' and '-' and "
                "must start with a letter or digit"
            )
        if len(fingerprint) != 64 or any(
            character not in "0123456789abcdefABCDEF" for character in fingerprint
        ):
            raise ValueError("fingerprint must be a 64-character hexadecimal SHA-256")
        self.run_directory = Path(run_directory).expanduser().resolve(strict=False)
        self.stage = stage
        self.fingerprint = fingerprint.lower()
        self.pipeline_version = pipeline_version
        self.parameters = dict(parameters or {})
        self.input_digests = dict(input_digests or {})
        self.backend_versions = dict(backend_versions or {})
        self.resume = resume
        self.verify_on_resume = verify_on_resume
        self.keep_failed = keep_failed
        self.reused = False
        self.manifest: StageManifest | None = None
        self._created_utc: str | None = None
        self._work_path: Path | None = None
        self.final_path = self.run_directory / "stages" / stage

    @property
    def path(self) -> Path:
        """Current work directory, or the validated completed directory on resume."""

        if self._work_path is None:
            raise RuntimeError("StageDirectory has not been entered")
        return self._work_path

    def output_path(self, relative: str | Path) -> Path:
        """Create parent directories and return a safe path inside active work."""

        if self.reused:
            raise StageExistsError("A reused stage is immutable")
        relative_path = Path(relative)
        if relative_path.is_absolute() or ".." in relative_path.parts or relative_path == Path("."):
            raise ValueError(f"Unsafe stage-relative path: {relative!r}")
        destination = self.path / relative_path
        destination.parent.mkdir(parents=True, exist_ok=True)
        return destination

    def __enter__(self) -> StageDirectory:
        stages = self.run_directory / "stages"
        stages.mkdir(parents=True, exist_ok=True)
        if self.final_path.exists():
            if not self.resume:
                raise StageExistsError(f"Stage already exists: {self.final_path}")
            self.manifest = validate_stage_directory(
                self.final_path,
                expected_fingerprint=self.fingerprint,
                expected_stage=self.stage,
                verify_checksums=self.verify_on_resume,
            )
            self.reused = True
            self._work_path = self.final_path
            return self
        self._created_utc = utc_now()
        self._work_path = stages / f".{self.stage}.tmp-{uuid4().hex}"
        self._work_path.mkdir(exist_ok=False)
        return self

    def _artifact_manifests(self) -> tuple[ArtifactManifest, ...]:
        artifacts: list[ArtifactManifest] = []
        for item in sorted(self.path.rglob("*"), key=lambda path: path.as_posix()):
            if item.is_symlink():
                raise ProvenanceError(f"Stage output may not be a symlink: {item}")
            if not item.is_file() or item.name in {MANIFEST_NAME, SUCCESS_NAME}:
                continue
            artifacts.append(
                ArtifactManifest(
                    path=item.relative_to(self.path).as_posix(),
                    sha256=sha256_file(item),
                    size_bytes=item.stat().st_size,
                    media_type=_guess_media_type(item),
                )
            )
        return tuple(artifacts)

    def _publish(self) -> None:
        artifacts = self._artifact_manifests()
        self.manifest = StageManifest(
            schema_version=STAGE_SCHEMA_VERSION,
            stage=self.stage,
            fingerprint=self.fingerprint,
            pipeline_version=self.pipeline_version,
            created_utc=self._created_utc or utc_now(),
            completed_utc=utc_now(),
            parameters=self.parameters,
            input_digests=self.input_digests,
            backend_versions=self.backend_versions,
            runtime=runtime_provenance(),
            artifacts=artifacts,
        )
        atomic_write_json(self.path / MANIFEST_NAME, self.manifest.as_dict())
        atomic_write_text(self.path / SUCCESS_NAME, self.fingerprint + "\n")
        if self.final_path.exists():
            raise StageExistsError(
                f"Another process published the stage concurrently: {self.final_path}"
            )
        # Both directories are siblings on one filesystem.  The destination is
        # required not to exist, preserving an already-completed run on all OSes.
        os.rename(self.path, self.final_path)
        self._work_path = self.final_path

    def __exit__(self, exc_type: Any, exc: BaseException | None, traceback: Any) -> bool:
        if self.reused:
            return False
        if exc is None:
            self._publish()
            return False
        if self.keep_failed and self._work_path is not None and self._work_path.exists():
            failure = {
                "stage": self.stage,
                "fingerprint": self.fingerprint,
                "failed_utc": utc_now(),
                "exception_type": type(exc).__name__,
                "message": str(exc),
            }
            atomic_write_json(self._work_path / "_FAILED.json", failure)
        return False


__all__ = [
    "ArtifactManifest",
    "ProvenanceError",
    "StageDirectory",
    "StageExistsError",
    "StageManifest",
    "StageValidationError",
    "atomic_write_bytes",
    "atomic_write_json",
    "atomic_write_text",
    "canonical_digest",
    "canonical_json",
    "compute_stage_fingerprint",
    "git_provenance",
    "runtime_provenance",
    "sha256_bytes",
    "sha256_file",
    "validate_stage_directory",
]
