"""Validated, portable configuration loading for Nanopore3.

Configuration files are deliberately strict: unknown keys are rejected so that a
misspelled threshold cannot silently change a run.  Every path is expanded and
resolved relative to the YAML file, never relative to the caller's current
working directory.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

import yaml


DNA_IUPAC = frozenset("ACGTRYSWKMBDHVN")


class ConfigError(ValueError):
    """Raised when a configuration file is malformed or internally inconsistent."""


def _mapping(value: Any, location: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ConfigError(f"{location} must be a mapping")
    return value


def _reject_unknown(
    value: Mapping[str, Any], allowed: set[str], location: str
) -> None:
    unknown = sorted(set(value) - allowed)
    if unknown:
        raise ConfigError(f"Unknown key(s) in {location}: {', '.join(unknown)}")


def _required(value: Mapping[str, Any], key: str, location: str) -> Any:
    if key not in value:
        raise ConfigError(f"Missing required key {location}.{key}")
    return value[key]


def _positive_int(value: Any, location: str, *, minimum: int = 1) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ConfigError(f"{location} must be an integer >= {minimum}")
    return value


def _number(
    value: Any, location: str, *, minimum: float, maximum: float
) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigError(f"{location} must be a number")
    result = float(value)
    if not minimum <= result <= maximum:
        raise ConfigError(
            f"{location} must be between {minimum:g} and {maximum:g}, inclusive"
        )
    return result


def _nonempty_string(value: Any, location: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigError(f"{location} must be a non-empty string")
    return value.strip()


def _path(value: Any, base_dir: Path, location: str) -> Path:
    raw = Path(_nonempty_string(value, location)).expanduser()
    if not raw.is_absolute():
        raw = base_dir / raw
    return raw.resolve(strict=False)


def _dna(value: Any, location: str) -> str:
    sequence = _nonempty_string(value, location).upper()
    invalid = sorted(set(sequence) - DNA_IUPAC)
    if invalid:
        raise ConfigError(
            f"{location} contains unsupported DNA symbol(s): {''.join(invalid)}"
        )
    return sequence


@dataclass(frozen=True, slots=True)
class InputSettings:
    """One raw sequencing input and the stable sample name assigned to it."""

    path: Path
    sample_id: str


@dataclass(frozen=True, slots=True)
class BarcodeSettings:
    """A named barcode set and the edit-distance search policy used for it."""

    sequences: Mapping[str, str] = field(default_factory=dict)
    max_edits: int = 4
    trim_bases: int = 0
    search_window: int = 400
    search_ends: tuple[str, ...] = ("head", "tail")
    allow_reverse_complement: bool = True
    minimum_margin: int = 1

    def __post_init__(self) -> None:
        if self.max_edits < 0:
            raise ConfigError("barcode max_edits must be >= 0")
        if self.trim_bases < 0:
            raise ConfigError("barcode trim_bases must be >= 0")
        if self.search_window < 1:
            raise ConfigError("barcode search_window must be >= 1")
        if self.minimum_margin < 0:
            raise ConfigError("barcode minimum_margin must be >= 0")
        if not self.search_ends or set(self.search_ends) - {"head", "tail"}:
            raise ConfigError(
                "barcode search_ends must contain only 'head' and/or 'tail'"
            )
        if len(set(self.search_ends)) != len(self.search_ends):
            raise ConfigError("barcode search_ends must not contain duplicates")
        if len(set(self.sequences)) != len(self.sequences):
            raise ConfigError("barcode identifiers must be unique")
        for barcode_id, sequence in self.sequences.items():
            _nonempty_string(barcode_id, "barcode identifier")
            _dna(sequence, f"barcode {barcode_id!r}")
            if self.trim_bases >= len(sequence):
                raise ConfigError(
                    f"barcode {barcode_id!r} is not longer than trim_bases"
                )


@dataclass(frozen=True, slots=True)
class ReferenceSettings:
    """Reference FASTA and conservative shortlist/alignment thresholds."""

    fasta: tuple[Path, ...]
    kmer_sizes: tuple[int, ...] = (15, 11, 9)
    max_kmer_owners: int = 10
    candidate_count: int = 5
    minimum_kmer_score: float = 2.0
    minimum_identity: float = 0.80
    minimum_query_coverage: float = 0.70
    minimum_reference_coverage: float = 0.70
    minimum_identity_margin: float = 0.02

    def __post_init__(self) -> None:
        if not self.fasta:
            raise ConfigError("references.fasta must contain at least one path")
        if not self.kmer_sizes:
            raise ConfigError("references.kmer_sizes must not be empty")
        if any(k < 3 for k in self.kmer_sizes):
            raise ConfigError("all references.kmer_sizes values must be >= 3")
        if len(set(self.kmer_sizes)) != len(self.kmer_sizes):
            raise ConfigError("references.kmer_sizes must not contain duplicates")
        for name in (
            "max_kmer_owners",
            "candidate_count",
        ):
            if getattr(self, name) < 1:
                raise ConfigError(f"references.{name} must be >= 1")
        if self.minimum_kmer_score < 0:
            raise ConfigError("references.minimum_kmer_score must be non-negative")
        for name in (
            "minimum_identity",
            "minimum_query_coverage",
            "minimum_reference_coverage",
            "minimum_identity_margin",
        ):
            value = getattr(self, name)
            if not 0.0 <= value <= 1.0:
                raise ConfigError(f"references.{name} must be between 0 and 1")


@dataclass(frozen=True, slots=True)
class LibrarySettings:
    """Assay/library boundaries that are independent of individual references."""

    name: str = "library"
    forward_motif: str | None = None
    reverse_motif: str | None = None
    motif_max_edits: int = 0
    minimum_read_length: int | None = None
    maximum_read_length: int | None = None
    minimum_mean_quality: float | None = None
    metadata: Path | None = None

    def __post_init__(self) -> None:
        _nonempty_string(self.name, "library.name")
        if self.forward_motif is not None:
            _dna(self.forward_motif, "library.forward_motif")
        if self.reverse_motif is not None:
            _dna(self.reverse_motif, "library.reverse_motif")
        if (self.forward_motif is None) != (self.reverse_motif is None):
            raise ConfigError(
                "library.forward_motif and library.reverse_motif must be set together"
            )
        if self.motif_max_edits < 0:
            raise ConfigError("library.motif_max_edits must be >= 0")
        if self.minimum_read_length is not None and self.minimum_read_length < 1:
            raise ConfigError("library.minimum_read_length must be >= 1")
        if self.maximum_read_length is not None and self.maximum_read_length < 1:
            raise ConfigError("library.maximum_read_length must be >= 1")
        if (
            self.minimum_read_length is not None
            and self.maximum_read_length is not None
            and self.minimum_read_length > self.maximum_read_length
        ):
            raise ConfigError(
                "library.minimum_read_length must not exceed maximum_read_length"
            )
        if self.minimum_mean_quality is not None and not 0 <= self.minimum_mean_quality <= 93:
            raise ConfigError("library.minimum_mean_quality must be between 0 and 93")


@dataclass(frozen=True, slots=True)
class ParallelSettings:
    """Portable process-level resource limits.

    ``jobs`` is the maximum number of independent Python worker processes.
    ``threads_per_job`` limits native/external-tool threads within each worker so
    callers can prevent nested oversubscription.
    """

    jobs: int = 1
    threads_per_job: int = 1
    chunk_reads: int = 25_000
    backend: str = "auto"

    def __post_init__(self) -> None:
        for name in ("jobs", "threads_per_job", "chunk_reads"):
            if getattr(self, name) < 1:
                raise ConfigError(f"parallel.{name} must be >= 1")
        if self.backend not in {"auto", "serial", "thread"}:
            raise ConfigError("parallel.backend must be auto, serial, or thread")


@dataclass(frozen=True, slots=True)
class ConsensusSettings:
    """Portable consensus thresholds and deterministic depth cap."""

    minimum_depth: int = 3
    maximum_reads: int = 100
    minimum_support: float = 0.60
    backend: str = "portable"

    def __post_init__(self) -> None:
        if self.minimum_depth < 1 or self.maximum_reads < 1:
            raise ConfigError("consensus depth settings must be positive")
        if self.minimum_depth > self.maximum_reads:
            raise ConfigError("consensus.minimum_depth must not exceed maximum_reads")
        if not 0.5 <= self.minimum_support <= 1.0:
            raise ConfigError("consensus.minimum_support must be between 0.5 and 1")
        if self.backend != "portable":
            raise ConfigError("v0.1 currently supports consensus.backend: portable")


@dataclass(frozen=True, slots=True)
class QcSettings:
    """Whole-consensus QC thresholds; region-aware QC can extend this contract."""

    minimum_identity: float = 0.98
    minimum_query_coverage: float = 0.95
    minimum_reference_coverage: float = 0.95
    length_tolerance: int = 10

    def __post_init__(self) -> None:
        for name in (
            "minimum_identity",
            "minimum_query_coverage",
            "minimum_reference_coverage",
        ):
            if not 0 <= getattr(self, name) <= 1:
                raise ConfigError(f"qc.{name} must be between 0 and 1")
        if self.length_tolerance < 0:
            raise ConfigError("qc.length_tolerance must be non-negative")


@dataclass(frozen=True, slots=True)
class PipelineConfig:
    """Fully resolved and validated Nanopore3 run configuration."""

    schema_version: int
    run_name: str
    output_root: Path
    inputs: tuple[InputSettings, ...]
    references: ReferenceSettings
    library: LibrarySettings = field(default_factory=LibrarySettings)
    plate_barcodes: BarcodeSettings = field(default_factory=BarcodeSettings)
    well_barcodes: BarcodeSettings = field(default_factory=BarcodeSettings)
    parallel: ParallelSettings = field(default_factory=ParallelSettings)
    consensus: ConsensusSettings = field(default_factory=ConsensusSettings)
    qc: QcSettings = field(default_factory=QcSettings)
    random_seed: int = 0
    source_path: Path | None = field(default=None, compare=False)

    def __post_init__(self) -> None:
        if self.schema_version != 1:
            raise ConfigError(
                f"Unsupported schema_version {self.schema_version!r}; expected 1"
            )
        _nonempty_string(self.run_name, "run_name")
        if not self.inputs:
            raise ConfigError("inputs must contain at least one FASTQ")
        sample_ids = [item.sample_id for item in self.inputs]
        if len(sample_ids) != len(set(sample_ids)):
            raise ConfigError("inputs sample_id values must be unique")
        input_paths = [item.path for item in self.inputs]
        if len(input_paths) != len(set(input_paths)):
            raise ConfigError("the same input FASTQ path is listed more than once")
        if isinstance(self.random_seed, bool) or not isinstance(self.random_seed, int):
            raise ConfigError("random_seed must be an integer")

    def as_dict(self) -> dict[str, Any]:
        """Return a serialization-friendly resolved representation."""

        return asdict(self)


# Public name used by the CLI and workflow layers.  ``PipelineConfig`` remains
# available as a descriptive compatibility alias for callers that prefer it.
RunConfig = PipelineConfig


def _parse_inputs(value: Any, base_dir: Path) -> tuple[InputSettings, ...]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ConfigError("inputs must be a list")
    result: list[InputSettings] = []
    for index, item in enumerate(value):
        location = f"inputs[{index}]"
        if isinstance(item, str):
            path = _path(item, base_dir, location)
            sample_id = path.name
            for suffix in (".gz", ".fastq", ".fq"):
                if sample_id.lower().endswith(suffix):
                    sample_id = sample_id[: -len(suffix)]
        else:
            mapping = _mapping(item, location)
            _reject_unknown(mapping, {"path", "sample_id"}, location)
            path = _path(_required(mapping, "path", location), base_dir, f"{location}.path")
            sample_id = _nonempty_string(
                mapping.get("sample_id", path.stem), f"{location}.sample_id"
            )
        result.append(InputSettings(path=path, sample_id=sample_id))
    return tuple(result)


def _parse_barcodes(value: Any, location: str) -> BarcodeSettings:
    if value is None:
        return BarcodeSettings()
    mapping = _mapping(value, location)
    allowed = {
        "sequences",
        "max_edits",
        "trim_bases",
        "search_window",
        "search_ends",
        "allow_reverse_complement",
        "minimum_margin",
    }
    _reject_unknown(mapping, allowed, location)
    sequences_value = mapping.get("sequences", {})
    sequences_mapping = _mapping(sequences_value, f"{location}.sequences")
    sequences: dict[str, str] = {}
    for barcode_id, sequence in sequences_mapping.items():
        clean_id = _nonempty_string(barcode_id, f"{location}.sequences key")
        sequences[clean_id] = _dna(
            sequence, f"{location}.sequences.{clean_id}"
        )
    ends_value = mapping.get("search_ends", ("head", "tail"))
    if not isinstance(ends_value, Sequence) or isinstance(ends_value, (str, bytes)):
        raise ConfigError(f"{location}.search_ends must be a list")
    ends = tuple(_nonempty_string(item, f"{location}.search_ends") for item in ends_value)
    allow_rc = mapping.get("allow_reverse_complement", True)
    if not isinstance(allow_rc, bool):
        raise ConfigError(f"{location}.allow_reverse_complement must be boolean")
    return BarcodeSettings(
        sequences=sequences,
        max_edits=_positive_int(
            mapping.get("max_edits", 4), f"{location}.max_edits", minimum=0
        ),
        trim_bases=_positive_int(
            mapping.get("trim_bases", 0), f"{location}.trim_bases", minimum=0
        ),
        search_window=_positive_int(
            mapping.get("search_window", 400), f"{location}.search_window"
        ),
        search_ends=ends,
        allow_reverse_complement=allow_rc,
        minimum_margin=_positive_int(
            mapping.get("minimum_margin", 1),
            f"{location}.minimum_margin",
            minimum=0,
        ),
    )


def _parse_references(value: Any, base_dir: Path) -> ReferenceSettings:
    location = "references"
    mapping = _mapping(value, location)
    allowed = {
        "fasta",
        "kmer_sizes",
        "max_kmer_owners",
        "candidate_count",
        "minimum_kmer_score",
        "minimum_identity",
        "minimum_query_coverage",
        "minimum_reference_coverage",
        "minimum_identity_margin",
    }
    _reject_unknown(mapping, allowed, location)
    fasta_value = _required(mapping, "fasta", location)
    if isinstance(fasta_value, str):
        fasta_items = [fasta_value]
    elif isinstance(fasta_value, Sequence):
        fasta_items = list(fasta_value)
    else:
        raise ConfigError("references.fasta must be a path or list of paths")
    fasta = tuple(
        _path(item, base_dir, f"references.fasta[{index}]")
        for index, item in enumerate(fasta_items)
    )
    kmer_value = mapping.get("kmer_sizes", (15, 11, 9))
    if not isinstance(kmer_value, Sequence) or isinstance(kmer_value, (str, bytes)):
        raise ConfigError("references.kmer_sizes must be a list")
    kmer_sizes = tuple(
        _positive_int(item, f"references.kmer_sizes[{index}]", minimum=3)
        for index, item in enumerate(kmer_value)
    )
    return ReferenceSettings(
        fasta=fasta,
        kmer_sizes=kmer_sizes,
        max_kmer_owners=_positive_int(
            mapping.get("max_kmer_owners", 10), "references.max_kmer_owners"
        ),
        candidate_count=_positive_int(
            mapping.get("candidate_count", 5), "references.candidate_count"
        ),
        minimum_kmer_score=_number(
            mapping.get("minimum_kmer_score", 2.0),
            "references.minimum_kmer_score",
            minimum=0,
            maximum=1_000_000,
        ),
        minimum_identity=_number(
            mapping.get("minimum_identity", 0.80),
            "references.minimum_identity",
            minimum=0,
            maximum=1,
        ),
        minimum_query_coverage=_number(
            mapping.get("minimum_query_coverage", 0.70),
            "references.minimum_query_coverage",
            minimum=0,
            maximum=1,
        ),
        minimum_reference_coverage=_number(
            mapping.get("minimum_reference_coverage", 0.70),
            "references.minimum_reference_coverage",
            minimum=0,
            maximum=1,
        ),
        minimum_identity_margin=_number(
            mapping.get("minimum_identity_margin", 0.02),
            "references.minimum_identity_margin",
            minimum=0,
            maximum=1,
        ),
    )


def _parse_library(value: Any, base_dir: Path) -> LibrarySettings:
    if value is None:
        return LibrarySettings()
    location = "library"
    mapping = _mapping(value, location)
    allowed = {
        "name",
        "forward_motif",
        "reverse_motif",
        "motif_max_edits",
        "minimum_read_length",
        "maximum_read_length",
        "minimum_mean_quality",
        "metadata",
    }
    _reject_unknown(mapping, allowed, location)
    forward = mapping.get("forward_motif")
    reverse = mapping.get("reverse_motif")
    return LibrarySettings(
        name=_nonempty_string(mapping.get("name", "library"), "library.name"),
        forward_motif=None if forward is None else _dna(forward, "library.forward_motif"),
        reverse_motif=None if reverse is None else _dna(reverse, "library.reverse_motif"),
        motif_max_edits=_positive_int(
            mapping.get("motif_max_edits", 0),
            "library.motif_max_edits",
            minimum=0,
        ),
        minimum_read_length=(
            None
            if mapping.get("minimum_read_length") is None
            else _positive_int(
                mapping["minimum_read_length"], "library.minimum_read_length"
            )
        ),
        maximum_read_length=(
            None
            if mapping.get("maximum_read_length") is None
            else _positive_int(
                mapping["maximum_read_length"], "library.maximum_read_length"
            )
        ),
        minimum_mean_quality=(
            None
            if mapping.get("minimum_mean_quality") is None
            else _number(
                mapping["minimum_mean_quality"],
                "library.minimum_mean_quality",
                minimum=0,
                maximum=93,
            )
        ),
        metadata=(
            None
            if mapping.get("metadata") is None
            else _path(mapping["metadata"], base_dir, "library.metadata")
        ),
    )


def _parse_parallel(value: Any) -> ParallelSettings:
    if value is None:
        return ParallelSettings()
    location = "parallel"
    mapping = _mapping(value, location)
    _reject_unknown(mapping, {"jobs", "threads_per_job", "chunk_reads", "backend"}, location)
    return ParallelSettings(
        jobs=_positive_int(mapping.get("jobs", 1), "parallel.jobs"),
        threads_per_job=_positive_int(
            mapping.get("threads_per_job", 1), "parallel.threads_per_job"
        ),
        chunk_reads=_positive_int(
            mapping.get("chunk_reads", 25_000), "parallel.chunk_reads"
        ),
        backend=_nonempty_string(mapping.get("backend", "auto"), "parallel.backend"),
    )


def _parse_consensus(value: Any) -> ConsensusSettings:
    if value is None:
        return ConsensusSettings()
    location = "consensus"
    mapping = _mapping(value, location)
    _reject_unknown(
        mapping,
        {"minimum_depth", "maximum_reads", "minimum_support", "backend"},
        location,
    )
    return ConsensusSettings(
        minimum_depth=_positive_int(mapping.get("minimum_depth", 3), "consensus.minimum_depth"),
        maximum_reads=_positive_int(mapping.get("maximum_reads", 100), "consensus.maximum_reads"),
        minimum_support=_number(
            mapping.get("minimum_support", 0.60),
            "consensus.minimum_support",
            minimum=0.5,
            maximum=1,
        ),
        backend=_nonempty_string(mapping.get("backend", "portable"), "consensus.backend"),
    )


def _parse_qc(value: Any) -> QcSettings:
    if value is None:
        return QcSettings()
    location = "qc"
    mapping = _mapping(value, location)
    _reject_unknown(
        mapping,
        {
            "minimum_identity",
            "minimum_query_coverage",
            "minimum_reference_coverage",
            "length_tolerance",
        },
        location,
    )
    return QcSettings(
        minimum_identity=_number(mapping.get("minimum_identity", 0.98), "qc.minimum_identity", minimum=0, maximum=1),
        minimum_query_coverage=_number(mapping.get("minimum_query_coverage", 0.95), "qc.minimum_query_coverage", minimum=0, maximum=1),
        minimum_reference_coverage=_number(mapping.get("minimum_reference_coverage", 0.95), "qc.minimum_reference_coverage", minimum=0, maximum=1),
        length_tolerance=_positive_int(mapping.get("length_tolerance", 10), "qc.length_tolerance", minimum=0),
    )


def load_config(path: str | Path) -> PipelineConfig:
    """Load and validate a Nanopore3 YAML configuration.

    The loader validates structure and value ranges but intentionally does not
    require inputs to exist.  A separate ingest validation step can therefore
    produce one complete, user-friendly report for missing and malformed inputs.
    """

    source_path = Path(path).expanduser().resolve(strict=True)
    if not source_path.is_file():
        raise ConfigError(f"Configuration is not a file: {source_path}")
    try:
        with source_path.open("r", encoding="utf-8") as handle:
            raw = yaml.safe_load(handle)
    except yaml.YAMLError as exc:
        raise ConfigError(f"Invalid YAML in {source_path}: {exc}") from exc
    root = _mapping(raw, "configuration")
    allowed = {
        "schema_version",
        "run_name",
        "output_root",
        "inputs",
        "references",
        "library",
        "barcodes",
        "parallel",
        "consensus",
        "qc",
        "random_seed",
    }
    _reject_unknown(root, allowed, "configuration")
    base_dir = source_path.parent
    barcodes_value = root.get("barcodes", {})
    barcodes = _mapping(barcodes_value, "barcodes")
    _reject_unknown(barcodes, {"plate", "well"}, "barcodes")
    schema_version = _positive_int(
        root.get("schema_version", 1), "schema_version"
    )
    random_seed = root.get("random_seed", 0)
    if isinstance(random_seed, bool) or not isinstance(random_seed, int):
        raise ConfigError("random_seed must be an integer")
    return PipelineConfig(
        schema_version=schema_version,
        run_name=_nonempty_string(root.get("run_name", source_path.stem), "run_name"),
        output_root=_path(
            _required(root, "output_root", "configuration"),
            base_dir,
            "output_root",
        ),
        inputs=_parse_inputs(_required(root, "inputs", "configuration"), base_dir),
        references=_parse_references(
            _required(root, "references", "configuration"), base_dir
        ),
        library=_parse_library(root.get("library"), base_dir),
        plate_barcodes=_parse_barcodes(barcodes.get("plate"), "barcodes.plate"),
        well_barcodes=_parse_barcodes(barcodes.get("well"), "barcodes.well"),
        parallel=_parse_parallel(root.get("parallel")),
        consensus=_parse_consensus(root.get("consensus")),
        qc=_parse_qc(root.get("qc")),
        random_seed=random_seed,
        source_path=source_path,
    )


__all__ = [
    "BarcodeSettings",
    "ConfigError",
    "ConsensusSettings",
    "InputSettings",
    "LibrarySettings",
    "ParallelSettings",
    "PipelineConfig",
    "QcSettings",
    "ReferenceSettings",
    "RunConfig",
    "load_config",
]
