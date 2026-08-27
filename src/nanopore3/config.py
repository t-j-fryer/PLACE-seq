"""Validated, portable configuration loading for Nanopore3.

Configuration files are deliberately strict: unknown keys are rejected so that a
misspelled threshold cannot silently change a run.  Every path is expanded and
resolved relative to the YAML file, never relative to the caller's current
working directory.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
import os
import re
from pathlib import Path
from typing import Any, Mapping, Sequence

import yaml

from .barcodes import BarcodeRegistryError, read_barcode_panel


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


def _boolean(value: Any, location: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigError(f"{location} must be true or false")
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


_VARIABLE = re.compile(r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}")


def expand_variables(value: str, location: str) -> str:
    """Substitute ``${NAME}`` from the environment.

    A configuration describes an experiment; where that experiment's data sits is
    a property of the machine reading it.  Writing the second into the first is
    what makes a config unshareable - and, once committed, what puts somebody's
    home directory in a public repository.  A named variable keeps the two apart:
    the config says *which* dataset, the environment says *where*.

    An unset variable is an error naming it, never an empty string: silently
    expanding to "" would produce a path like "/AI_DBTL.fastq" and a
    file-not-found three lines later that says nothing about the cause.
    """

    missing: list[str] = []

    def replace(match: re.Match[str]) -> str:
        name = match.group(1)
        resolved = os.environ.get(name)
        if resolved is None:
            missing.append(name)
            return ""
        return resolved

    expanded = _VARIABLE.sub(replace, value)
    if missing:
        names = ", ".join(sorted(set(missing)))
        raise ConfigError(
            f"{location} refers to unset environment variable(s): {names}. "
            f"Set them for this machine, for example:\n"
            + "\n".join(f"  export {name}=/path/to/data" for name in sorted(set(missing)))
        )
    return expanded


def _path(value: Any, base_dir: Path, location: str) -> Path:
    text = expand_variables(_nonempty_string(value, location), location)
    raw = Path(text).expanduser()
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
    decision_policy: str = "best_margin"
    registry_csv: Path | None = None
    family_id: str | None = None
    registry_sha256: str | None = None

    def __post_init__(self) -> None:
        if self.max_edits < 0:
            raise ConfigError("barcode max_edits must be >= 0")
        if self.trim_bases < 0:
            raise ConfigError("barcode trim_bases must be >= 0")
        if self.search_window < 1:
            raise ConfigError("barcode search_window must be >= 1")
        if self.minimum_margin < 0:
            raise ConfigError("barcode minimum_margin must be >= 0")
        if self.decision_policy not in {
            "best_margin",
            "legacy_unique_threshold",
            "legacy_unique_best",
        }:
            raise ConfigError(
                "barcode decision_policy must be best_margin, "
                "legacy_unique_threshold, or legacy_unique_best"
            )
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
            if 2 * self.trim_bases >= len(sequence):
                raise ConfigError(
                    f"barcode {barcode_id!r} is not longer than twice trim_bases"
                )


@dataclass(frozen=True, slots=True)
class ReferenceSettings:
    """Reference FASTA and conservative shortlist/alignment thresholds."""

    fasta: tuple[Path, ...]
    # Optional oPool *_FULL_INFO.csv design table. When supplied it is the
    # authoritative source of each gene's assembly block, joined by sequence
    # because the FASTA and the design table use different identifiers.
    fragments_csv: Path | None = None
    kmer_sizes: tuple[int, ...] = (15, 11, 9)
    max_kmer_owners: int = 10
    candidate_count: int = 5
    minimum_kmer_score: float = 2.0
    minimum_identity: float = 0.80
    minimum_query_coverage: float = 0.70
    minimum_reference_coverage: float = 0.70
    minimum_identity_margin: float = 0.02
    # How hard to look when no k-mer size proposes an acceptable reference.
    # "kmer" widens the shortlist using any k-mer evidence; "all" restores the
    # exhaustive whole-library sweep; "none" accepts the shortlist verdict.
    rescue_policy: str = "kmer"
    rescue_candidates: int = 25
    # Optional per-library overrides.  One run can carry libraries built on
    # different constructs, whose amplicons end at different constant regions;
    # a single global motif pair would extract the wrong span for some of them.
    # Unset values fall back to the run-level library/qc settings.
    forward_motif: str | None = None
    reverse_motif: str | None = None
    motif_max_edits: int | None = None
    qc_upstream_constant: str | None = None
    qc_downstream_constant: str | None = None
    # Full-length mode.  An oligo-pool tool emits only the variable insert, but
    # the sequenced molecule is the whole amplicon, so the constant regions are
    # joined on at load and the reconstructed region is bounded by the primer
    # sites at their outer ends.  Give the two constant regions directly, or give
    # one example construct - a binder already in the vector - and let them be
    # derived from it.
    flanks_upstream: str | None = None
    flanks_downstream: str | None = None
    flanks_template: Path | None = None
    # Third route: the FASTA already holds assembled constructs, and the shared
    # backbone is found as their longest common prefix and suffix. Costs the user
    # nothing to declare and unlocks per-region accuracy, insert-scoped thresholds
    # and insert-scoped chimera detection for a library that was never expressed
    # as inserts-plus-backbone.
    flanks_derive: bool = False
    flanks_anchor_length: int = 20
    # Apply minimum_identity and minimum_query_coverage to the insert region
    # rather than the whole amplicon. Only meaningful for a full-length library,
    # and on by default there, because the constant flanks are ~70% of every
    # reference: a read carrying 250 nt of foreign insert still scored 0.93
    # identity and 0.86 query coverage over the amplicon, clearing floors of 0.80
    # and 0.70 that were set when the reference was the insert. Over the insert
    # alone the same read covers 0.53 and fails.
    insert_thresholds: bool = True

    @property
    def full_length(self) -> bool:
        """Whether references are flanked into whole amplicons."""

        return bool(
            (self.flanks_upstream and self.flanks_downstream)
            or self.flanks_template
            or self.flanks_derive
        )

    def __post_init__(self) -> None:
        if not self.fasta:
            raise ConfigError("references.fasta must contain at least one path")
        if len(set(self.fasta)) != len(self.fasta):
            raise ConfigError("references.fasta must not contain duplicate paths")
        if not self.kmer_sizes:
            raise ConfigError("references.kmer_sizes must not be empty")
        if any(k < 3 for k in self.kmer_sizes):
            raise ConfigError("all references.kmer_sizes values must be >= 3")
        if len(set(self.kmer_sizes)) != len(self.kmer_sizes):
            raise ConfigError("references.kmer_sizes must not contain duplicates")
        # The cascade consults the most specific k first and rescues with the
        # most sensitive one, so the declared order is load-bearing.
        if list(self.kmer_sizes) != sorted(self.kmer_sizes, reverse=True):
            raise ConfigError(
                "references.kmer_sizes must be listed in descending order, "
                "most specific first"
            )
        for name in (
            "max_kmer_owners",
            "candidate_count",
            "rescue_candidates",
        ):
            if getattr(self, name) < 1:
                raise ConfigError(f"references.{name} must be >= 1")
        for name in (
            "forward_motif",
            "reverse_motif",
            "qc_upstream_constant",
            "qc_downstream_constant",
        ):
            value = getattr(self, name)
            if value is not None:
                _dna(value, f"references.{name}")
        # These are overrides layered on the run-level library motifs, so either
        # may be given alone; the unset one falls back. Requiring both together
        # would force a library that only ends differently to restate the shared
        # 5' motif, inviting the two copies to drift apart.
        if (self.qc_upstream_constant is None) != (self.qc_downstream_constant is None):
            raise ConfigError(
                "references.qc_upstream_constant and references.qc_downstream_constant "
                "must be supplied together"
            )
        if self.motif_max_edits is not None and self.motif_max_edits < 0:
            raise ConfigError("references.motif_max_edits must be >= 0")
        for name in ("flanks_upstream", "flanks_downstream"):
            value = getattr(self, name)
            if value is not None:
                _dna(value, f"references.{name}")
        explicit = (self.flanks_upstream is not None, self.flanks_downstream is not None)
        if any(explicit) and not all(explicit):
            raise ConfigError(
                "references.flanks.upstream and references.flanks.downstream must be "
                "supplied together, or use references.flanks.template instead"
            )
        routes = [
            name
            for name, chosen in (
                ("upstream/downstream", all(explicit)),
                ("template", self.flanks_template is not None),
                ("derive", self.flanks_derive),
            )
            if chosen
        ]
        if len(routes) > 1:
            raise ConfigError(
                "references.flanks accepts exactly one of upstream/downstream, "
                f"template or derive; got {', '.join(routes)}. Two definitions of "
                "the same constant regions can disagree"
            )
        if self.flanks_anchor_length < 8:
            raise ConfigError("references.flanks.anchor_length must be >= 8")
        if self.rescue_policy not in ("none", "kmer", "all"):
            raise ConfigError(
                "references.rescue_policy must be 'none', 'kmer', or 'all'"
            )
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
    # Basecallers occasionally emit zero-length reads. They are rejected as
    # malformed input by default, because a run should not silently depend on
    # how a parser treats them. Enabling this keeps them as a counted
    # `empty_read` demultiplexing state instead of failing the whole file.
    allow_empty_reads: bool = False

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
    """Portable worker and nested-thread resource limits.

    ``jobs`` is the maximum number of independent Python workers; **0 means use
    whatever this machine has**, so one configuration runs sensibly on a workstation
    and on a two-core hosted notebook without being edited.  Any value is bounded by
    the detected CPU budget either way.
    ``threads_per_job`` limits native/external-tool threads within each worker so
    callers can prevent nested oversubscription.
    """

    jobs: int = 1
    threads_per_job: int = 1
    chunk_reads: int = 25_000
    backend: str = "auto"

    def __post_init__(self) -> None:
        if self.jobs < 0:
            raise ConfigError("parallel.jobs must be >= 0 (0 means detect)")
        for name in ("threads_per_job", "chunk_reads"):
            if getattr(self, name) < 1:
                raise ConfigError(f"parallel.{name} must be >= 1")
        if self.backend not in {"auto", "serial", "thread", "process"}:
            raise ConfigError(
                "parallel.backend must be auto, serial, thread, or process"
            )


@dataclass(frozen=True, slots=True)
class ConsensusSettings:
    """Portable consensus thresholds and deterministic depth cap."""

    minimum_depth: int = 3
    # A cap discards depth exactly where a well is deepest, which is where the
    # evidence for a second allele lives. 300 sits above the 90th percentile of
    # reads per clone in the 260608 run while keeping memory bounded: the cap also
    # limits how many reads are held while streaming, and 0 (unbounded) costs
    # ~1.5 GB there. Beyond a few hundred reads a binomial test learns nothing
    # more, so this is depth where it counts rather than depth for its own sake.
    maximum_reads: int = 300
    # Retained for the "majority" caller used when polishing a draft. The
    # statistical caller does not use it: a flat fraction cannot tell 61% against
    # 32% (two alleles) from 61% against scattered error (one allele plus noise).
    minimum_support: float = 0.60
    # A base observation below this Phred is not counted. Quality gates whether an
    # observation is admitted; it does not scale a vote, which would make the
    # support fraction a fraction of weighted votes rather than of reads.
    minimum_base_quality: int = 10
    # Adjusted p below which an allele is held to exceed the group's own
    # background error rate.
    significance: float = 0.05
    # And the share of reads an allele needs before it counts as an allele at all.
    # The binomial test assumes independent errors at a uniform rate; nanopore
    # errors cluster in homopolymers, so the test alone flags systematic error as
    # a second allele. Both bars are required.
    minimum_minor_fraction: float = 0.20
    backend: str = "portable"

    def __post_init__(self) -> None:
        if self.minimum_depth < 1 or self.maximum_reads < 0:
            raise ConfigError("consensus depth settings must be positive")
        if self.maximum_reads and self.minimum_depth > self.maximum_reads:
            raise ConfigError("consensus.minimum_depth must not exceed maximum_reads")
        if not 0.5 <= self.minimum_support <= 1.0:
            raise ConfigError("consensus.minimum_support must be between 0.5 and 1")
        if self.minimum_base_quality < 0:
            raise ConfigError("consensus.minimum_base_quality must be >= 0")
        if not 0.0 < self.significance < 1.0:
            raise ConfigError("consensus.significance must be between 0 and 1")
        if not 0.0 < self.minimum_minor_fraction <= 0.5:
            raise ConfigError(
                "consensus.minimum_minor_fraction must be between 0 and 0.5"
            )
        if self.backend not in {"portable", "mafft_spoa"}:
            raise ConfigError(
                "consensus.backend must be portable or mafft_spoa"
            )


@dataclass(frozen=True, slots=True)
class QcSettings:
    """Whole-consensus QC thresholds; region-aware QC can extend this contract."""

    minimum_identity: float = 0.98
    minimum_query_coverage: float = 0.95
    minimum_reference_coverage: float = 0.95
    length_tolerance: int = 10
    # Constant sequence flanking the consensus in the full open reading frame.
    # `upstream_constant` must begin at the start codon. Supplying both enables
    # the reading-frame and internal-stop checks, which are otherwise reported
    # as not_evaluable because their answer would depend on an assumed frame.
    upstream_constant: str | None = None
    downstream_constant: str | None = None

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
        if (self.upstream_constant is None) != (self.downstream_constant is None):
            raise ConfigError(
                "qc.upstream_constant and qc.downstream_constant must be supplied "
                "together; the reading frame is undefined without both"
            )
        if self.upstream_constant is not None:
            _dna(self.upstream_constant, "qc.upstream_constant")
            _dna(self.downstream_constant, "qc.downstream_constant")
            if not self.upstream_constant.upper().startswith(("ATG", "GTG", "TTG")):
                raise ConfigError(
                    "qc.upstream_constant must begin at the start codon so the "
                    "reading frame is unambiguous"
                )


@dataclass(frozen=True, slots=True)
class ChimeraSettings:
    """Positional detection of chimeric clones and their consensus sequences.

    A chimeric molecule is a real clone in a real well, so its reads are worth a
    consensus.  ``write_pcr_origin`` is false because a chimera whose parents sit
    in different assembly blocks can only have formed during PCR: it is an
    artefact, and writing it as a FASTA beside genuine clones invites misreading.
    Such groups are still counted.
    """

    enabled: bool = False
    window: int = 90
    step: int = 45
    minimum_depth: int | None = None
    write_pcr_origin: bool = False
    # Which chimeric reads stop feeding the reference-guided consensuses.
    # "written" excludes only reads belonging to a chimeric clone that was
    # written, so every excluded read still contributes to a consensus - itself.
    # "all" also excludes PCR-origin reads, which are artefacts and are not
    # written, so those reads then contribute to nothing. "none" lets one molecule
    # appear twice: once as itself and once as a mutated version of one parent.
    #
    # Default "written". Four clones removed by it looked genuine - one reported
    # 100% identity from 146 reads - and were checked one at a time: every one is
    # a chimera of two designs only ~50% identical to each other, in reads 500-590
    # nt long against designs of 285-405 nt. They reported perfect identity only
    # because a reference-guided consensus represents the part of a molecule that
    # aligns and silently drops the rest. Removing them is right.
    exclude_reads: str = "written"

    def __post_init__(self) -> None:
        if self.exclude_reads not in ("written", "all", "none"):
            raise ConfigError(
                "chimera.exclude_reads must be 'written', 'all', or 'none'"
            )
        if self.window < 1 or self.step < 1:
            raise ConfigError("chimera.window and chimera.step must be >= 1")
        if self.step > self.window:
            raise ConfigError("chimera.step must not exceed chimera.window")
        if self.minimum_depth is not None and self.minimum_depth < 1:
            raise ConfigError("chimera.minimum_depth must be >= 1")


@dataclass(frozen=True, slots=True)
class CompressedPcrSettings:
    """Recover the source culture plate when several were pooled for colony PCR.

    ``pcr_plates`` lists the culture plates loaded into each colony PCR plate,
    keyed by plate barcode.  ``blocks`` lists the culture plate(s) each assembly
    block was picked into.  Together they let an assigned gene identify which
    pooled culture plate a read came from.
    """

    enabled: bool = False
    pcr_plates: Mapping[str, tuple[str, ...]] = field(default_factory=dict)
    # Nested by reference library, because block numbering restarts in each one.
    blocks: Mapping[str, Mapping[str, tuple[str, ...]]] = field(default_factory=dict)
    # Per colony PCR plate: "per_block" when the design put one clone per block,
    # so a well's expected diversity follows from the layout; "unspecified" for
    # scraped or otherwise unpredictable picking.
    clonality: Mapping[str, str] = field(default_factory=dict)
    # Applied to reference identifiers when no design table is configured.
    block_pattern: str = r"^Block_(\d+)_"

    def __post_init__(self) -> None:
        if not self.enabled:
            return
        if not self.pcr_plates:
            raise ConfigError("compressed_pcr.pcr_plates is required when enabled")
        if not self.blocks:
            raise ConfigError("compressed_pcr.blocks is required when enabled")
        try:
            re.compile(self.block_pattern)
        except re.error as exc:
            raise ConfigError(f"compressed_pcr.block_pattern is not a regex: {exc}") from exc
        for plate, sources in self.pcr_plates.items():
            _nonempty_string(plate, "compressed_pcr.pcr_plates key")
            if not sources:
                raise ConfigError(
                    f"compressed_pcr.pcr_plates[{plate!r}] must list at least one culture plate"
                )
        for library, by_block in self.blocks.items():
            _nonempty_string(library, "compressed_pcr.blocks library key")
            if not by_block:
                raise ConfigError(
                    f"compressed_pcr.blocks[{library!r}] must map at least one block"
                )
            for block, plates in by_block.items():
                _nonempty_string(block, "compressed_pcr.blocks block key")
                if not plates:
                    raise ConfigError(
                        f"compressed_pcr.blocks[{library!r}][{block!r}] must list "
                        "at least one culture plate"
                    )
        # Validate the layout itself: an unresolvable pooling design is an error
        # that can be caught now instead of appearing as ambiguous reads later.
        from .deconvolution import CompressedPcrPlan, DeconvolutionError

        try:
            CompressedPcrPlan(self.pcr_plates, self.blocks, self.clonality)
        except DeconvolutionError as exc:
            raise ConfigError(str(exc)) from exc


@dataclass(frozen=True, slots=True)
class PipelineConfig:
    """Fully resolved and validated Nanopore3 run configuration."""

    schema_version: int
    run_name: str
    output_root: Path
    inputs: tuple[InputSettings, ...]
    references: ReferenceSettings | None
    reference_libraries: Mapping[str, ReferenceSettings] = field(default_factory=dict)
    plate_reference_map: Mapping[str, str] = field(default_factory=dict)
    library: LibrarySettings = field(default_factory=LibrarySettings)
    plate_barcodes: BarcodeSettings = field(default_factory=BarcodeSettings)
    well_barcodes: BarcodeSettings = field(default_factory=BarcodeSettings)
    parallel: ParallelSettings = field(default_factory=ParallelSettings)
    consensus: ConsensusSettings = field(default_factory=ConsensusSettings)
    qc: QcSettings = field(default_factory=QcSettings)
    compressed_pcr: CompressedPcrSettings = field(default_factory=CompressedPcrSettings)
    chimera: ChimeraSettings = field(default_factory=ChimeraSettings)
    random_seed: int = 0
    # The graded per-plate tree is the output most runs are actually read from -
    # one FASTA per clone, named with its design and a single-word grade, so a
    # screening decision can be made from the file listing alone. It used to
    # require knowing about a separate script, which meant most runs never
    # produced it. Set false for a very large run where the file count matters.
    consensus_tree: bool = True
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
        if self.references is not None and self.reference_libraries:
            raise ConfigError(
                "configuration must use either references or reference_libraries, not both"
            )
        if self.references is None and not self.reference_libraries:
            raise ConfigError(
                "configuration must define references or reference_libraries"
            )
        for library_id in self.reference_libraries:
            _nonempty_string(library_id, "reference_libraries identifier")
        for plate_id, library_id in self.plate_reference_map.items():
            _nonempty_string(plate_id, "plate_reference_map identifier")
            clean_library_id = _nonempty_string(
                library_id, f"plate_reference_map.{plate_id}"
            )
            if clean_library_id not in self.reference_libraries:
                raise ConfigError(
                    f"plate_reference_map.{plate_id} names unknown reference "
                    f"library {clean_library_id!r}"
                )
        if self.references is not None and self.plate_reference_map:
            raise ConfigError(
                "plate_reference_map is only valid with reference_libraries"
            )

    @property
    def reference_sets(self) -> Mapping[str, ReferenceSettings]:
        """Return all reference settings under stable library identifiers.

        Legacy configurations with one ``references`` section are exposed as a
        single library named ``default``.  Multi-library mappings are sorted so
        manifests and stage fingerprints do not depend on YAML key order.
        """

        if self.references is not None:
            return {"default": self.references}
        return {
            library_id: self.reference_libraries[library_id]
            for library_id in sorted(self.reference_libraries)
        }

    def with_rescue_policy(self, rescue_policy: str) -> "PipelineConfig":
        """Return a copy whose reference libraries all use one rescue policy.

        This exists so a benchmark or audit can compare rescue policies against
        the very configuration a run used, without editing the run file.
        """

        if self.references is not None:
            return replace(
                self, references=replace(self.references, rescue_policy=rescue_policy)
            )
        return replace(
            self,
            reference_libraries={
                library_id: replace(settings, rescue_policy=rescue_policy)
                for library_id, settings in self.reference_libraries.items()
            },
        )

    def reference_library_id_for_plate(self, plate_barcode_id: str) -> str:
        """Resolve a plate barcode to its configured reference-library ID.

        A single reference set is an unambiguous fallback.  Configurations with
        multiple sets must map every plate encountered by the pipeline.
        """

        plate_id = _nonempty_string(plate_barcode_id, "plate barcode identifier")
        mapped = self.plate_reference_map.get(plate_id)
        if mapped is not None:
            return mapped
        reference_sets = self.reference_sets
        if len(reference_sets) == 1:
            return next(iter(reference_sets))
        raise KeyError(
            f"No reference library is configured for plate barcode {plate_id!r}"
        )

    def reference_library_for_plate(
        self, plate_barcode_id: str
    ) -> ReferenceSettings:
        """Return the reference settings applicable to one plate barcode."""

        return self.reference_sets[
            self.reference_library_id_for_plate(plate_barcode_id)
        ]

    def as_dict(self) -> dict[str, Any]:
        """Return a serialization-friendly resolved representation."""

        result = asdict(self)
        # The YAML's own location is loader context, not a scientific parameter.
        # All data paths above are already resolved, so retaining source_path would
        # make otherwise identical copied configs hash differently.
        result.pop("source_path", None)
        result["reference_libraries"] = {
            library_id: asdict(settings)
            for library_id, settings in self.reference_sets.items()
            if self.references is None
        }
        result["plate_reference_map"] = {
            plate_id: self.plate_reference_map[plate_id]
            for plate_id in sorted(self.plate_reference_map)
        }
        return result


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


def _parse_barcodes(value: Any, location: str, base_dir: Path | None = None) -> BarcodeSettings:
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
        "decision_policy",
        "registry_csv",
        "family_id",
    }
    _reject_unknown(mapping, allowed, location)
    has_inline = "sequences" in mapping
    has_registry = "registry_csv" in mapping or "family_id" in mapping
    if has_inline and has_registry:
        raise ConfigError(
            f"{location} must use either sequences or registry_csv/family_id, not both"
        )
    registry_path: Path | None = None
    family_id: str | None = None
    registry_sha256: str | None = None
    if has_registry:
        if "registry_csv" not in mapping or "family_id" not in mapping:
            raise ConfigError(
                f"{location}.registry_csv and {location}.family_id must be provided together"
            )
        raw_path = Path(_nonempty_string(mapping["registry_csv"], f"{location}.registry_csv")).expanduser()
        registry_path = (
            raw_path if raw_path.is_absolute() else (base_dir or Path.cwd()) / raw_path
        ).resolve(strict=False)
        family_id = _nonempty_string(mapping["family_id"], f"{location}.family_id")
        try:
            panel = read_barcode_panel(registry_path, family_id)
        except (BarcodeRegistryError, OSError) as exc:
            raise ConfigError(f"invalid {location} registry: {exc}") from exc
        sequences = panel.sequences
        registry_path = panel.source_path
        registry_sha256 = panel.source_sha256
    else:
        sequences_value = mapping.get("sequences", {})
        sequences_mapping = _mapping(sequences_value, f"{location}.sequences")
        sequences = {}
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
        decision_policy=_nonempty_string(
            mapping.get("decision_policy", "best_margin"),
            f"{location}.decision_policy",
        ),
        registry_csv=registry_path,
        family_id=family_id,
        registry_sha256=registry_sha256,
    )


def _parse_references(
    value: Any, base_dir: Path, *, location: str = "references"
) -> ReferenceSettings:
    mapping = _mapping(value, location)
    allowed = {
        "fasta",
        "fragments_csv",
        "kmer_sizes",
        "max_kmer_owners",
        "candidate_count",
        "minimum_kmer_score",
        "minimum_identity",
        "minimum_query_coverage",
        "minimum_reference_coverage",
        "minimum_identity_margin",
        "rescue_policy",
        "rescue_candidates",
        "forward_motif",
        "reverse_motif",
        "motif_max_edits",
        "qc_upstream_constant",
        "qc_downstream_constant",
        "flanks",
        "insert_thresholds",
    }
    _reject_unknown(mapping, allowed, location)
    flanks_value = mapping.get("flanks")
    flanks_fields: dict[str, Any] = {}
    if flanks_value is not None:
        flanks_mapping = _mapping(flanks_value, f"{location}.flanks")
        _reject_unknown(
            flanks_mapping,
            {"upstream", "downstream", "template", "derive", "anchor_length"},
            f"{location}.flanks",
        )
        for key in ("upstream", "downstream"):
            raw = flanks_mapping.get(key)
            flanks_fields[f"flanks_{key}"] = (
                None if raw is None else _dna(raw, f"{location}.flanks.{key}")
            )
        if flanks_mapping.get("derive") is not None:
            flanks_fields["flanks_derive"] = _boolean(
                flanks_mapping["derive"], f"{location}.flanks.derive"
            )
        template = flanks_mapping.get("template")
        flanks_fields["flanks_template"] = (
            None if template is None else _path(template, base_dir, f"{location}.flanks.template")
        )
        if flanks_mapping.get("anchor_length") is not None:
            flanks_fields["flanks_anchor_length"] = _positive_int(
                flanks_mapping["anchor_length"], f"{location}.flanks.anchor_length", minimum=8
            )
    fasta_value = _required(mapping, "fasta", location)
    if isinstance(fasta_value, str):
        fasta_items = [fasta_value]
    elif isinstance(fasta_value, Sequence):
        fasta_items = list(fasta_value)
    else:
        raise ConfigError(f"{location}.fasta must be a path or list of paths")
    fasta = tuple(
        _path(item, base_dir, f"{location}.fasta[{index}]")
        for index, item in enumerate(fasta_items)
    )
    kmer_value = mapping.get("kmer_sizes", (15, 11, 9))
    if not isinstance(kmer_value, Sequence) or isinstance(kmer_value, (str, bytes)):
        raise ConfigError(f"{location}.kmer_sizes must be a list")
    kmer_sizes = tuple(
        _positive_int(item, f"{location}.kmer_sizes[{index}]", minimum=3)
        for index, item in enumerate(kmer_value)
    )
    fragments_value = mapping.get("fragments_csv")
    optional_dna = {
        name: (
            None if mapping.get(name) is None else _dna(mapping[name], f"{location}.{name}")
        )
        for name in (
            "forward_motif",
            "reverse_motif",
            "qc_upstream_constant",
            "qc_downstream_constant",
        )
    }
    return ReferenceSettings(
        **flanks_fields,
        **optional_dna,
        insert_thresholds=_boolean(
            mapping.get("insert_thresholds", True), f"{location}.insert_thresholds"
        ),
        motif_max_edits=(
            None
            if mapping.get("motif_max_edits") is None
            else _positive_int(mapping["motif_max_edits"], f"{location}.motif_max_edits", minimum=0)
        ),
        fasta=fasta,
        fragments_csv=(
            None
            if fragments_value is None
            else _path(fragments_value, base_dir, f"{location}.fragments_csv")
        ),
        rescue_policy=_nonempty_string(
            mapping.get("rescue_policy", "kmer"), f"{location}.rescue_policy"
        ),
        rescue_candidates=_positive_int(
            mapping.get("rescue_candidates", 25), f"{location}.rescue_candidates"
        ),
        kmer_sizes=kmer_sizes,
        max_kmer_owners=_positive_int(
            mapping.get("max_kmer_owners", 10), f"{location}.max_kmer_owners"
        ),
        candidate_count=_positive_int(
            mapping.get("candidate_count", 5), f"{location}.candidate_count"
        ),
        minimum_kmer_score=_number(
            mapping.get("minimum_kmer_score", 2.0),
            f"{location}.minimum_kmer_score",
            minimum=0,
            maximum=1_000_000,
        ),
        minimum_identity=_number(
            mapping.get("minimum_identity", 0.80),
            f"{location}.minimum_identity",
            minimum=0,
            maximum=1,
        ),
        minimum_query_coverage=_number(
            mapping.get("minimum_query_coverage", 0.70),
            f"{location}.minimum_query_coverage",
            minimum=0,
            maximum=1,
        ),
        minimum_reference_coverage=_number(
            mapping.get("minimum_reference_coverage", 0.70),
            f"{location}.minimum_reference_coverage",
            minimum=0,
            maximum=1,
        ),
        minimum_identity_margin=_number(
            mapping.get("minimum_identity_margin", 0.02),
            f"{location}.minimum_identity_margin",
            minimum=0,
            maximum=1,
        ),
    )


def _parse_reference_libraries(
    value: Any, base_dir: Path
) -> Mapping[str, ReferenceSettings]:
    location = "reference_libraries"
    mapping = _mapping(value, location)
    if not mapping:
        raise ConfigError("reference_libraries must contain at least one library")
    result: dict[str, ReferenceSettings] = {}
    for raw_library_id in sorted(mapping, key=str):
        library_id = _nonempty_string(raw_library_id, f"{location} identifier")
        result[library_id] = _parse_references(
            mapping[raw_library_id],
            base_dir,
            location=f"{location}.{library_id}",
        )
    return result


def _parse_plate_reference_map(value: Any) -> Mapping[str, str]:
    location = "plate_reference_map"
    mapping = _mapping(value, location)
    result: dict[str, str] = {}
    for raw_plate_id in sorted(mapping, key=str):
        plate_id = _nonempty_string(raw_plate_id, f"{location} identifier")
        result[plate_id] = _nonempty_string(
            mapping[raw_plate_id], f"{location}.{plate_id}"
        )
    return result


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
        "allow_empty_reads",
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
        allow_empty_reads=_boolean(
            mapping.get("allow_empty_reads", False), "library.allow_empty_reads"
        ),
    )


def _parse_parallel(value: Any) -> ParallelSettings:
    if value is None:
        return ParallelSettings()
    location = "parallel"
    mapping = _mapping(value, location)
    _reject_unknown(mapping, {"jobs", "threads_per_job", "chunk_reads", "backend"}, location)
    return ParallelSettings(
        jobs=_positive_int(mapping.get("jobs", 1), "parallel.jobs", minimum=0),
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
        {
            "minimum_depth",
            "maximum_reads",
            "minimum_support",
            "minimum_base_quality",
            "significance",
            "minimum_minor_fraction",
            "backend",
        },
        location,
    )
    return ConsensusSettings(
        minimum_depth=_positive_int(mapping.get("minimum_depth", 3), "consensus.minimum_depth"),
        maximum_reads=_positive_int(
            mapping.get("maximum_reads", 300), "consensus.maximum_reads", minimum=0
        ),
        minimum_support=_number(
            mapping.get("minimum_support", 0.60),
            "consensus.minimum_support",
            minimum=0.5,
            maximum=1,
        ),
        minimum_base_quality=_positive_int(
            mapping.get("minimum_base_quality", 10),
            "consensus.minimum_base_quality",
            minimum=0,
        ),
        significance=_number(
            mapping.get("significance", 0.05), "consensus.significance", minimum=0, maximum=1
        ),
        minimum_minor_fraction=_number(
            mapping.get("minimum_minor_fraction", 0.20),
            "consensus.minimum_minor_fraction",
            minimum=0,
            maximum=0.5,
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
            "upstream_constant",
            "downstream_constant",
        },
        location,
    )
    upstream = mapping.get("upstream_constant")
    downstream = mapping.get("downstream_constant")
    return QcSettings(
        upstream_constant=(
            None if upstream is None else _dna(upstream, "qc.upstream_constant")
        ),
        downstream_constant=(
            None if downstream is None else _dna(downstream, "qc.downstream_constant")
        ),
        minimum_identity=_number(mapping.get("minimum_identity", 0.98), "qc.minimum_identity", minimum=0, maximum=1),
        minimum_query_coverage=_number(mapping.get("minimum_query_coverage", 0.95), "qc.minimum_query_coverage", minimum=0, maximum=1),
        minimum_reference_coverage=_number(mapping.get("minimum_reference_coverage", 0.95), "qc.minimum_reference_coverage", minimum=0, maximum=1),
        length_tolerance=_positive_int(mapping.get("length_tolerance", 10), "qc.length_tolerance", minimum=0),
    )


def _parse_chimera(value: Any) -> ChimeraSettings:
    if value is None:
        return ChimeraSettings()
    location = "chimera"
    mapping = _mapping(value, location)
    _reject_unknown(
        mapping,
        {
            "enabled",
            "window",
            "step",
            "minimum_depth",
            "write_pcr_origin",
            "exclude_reads",
        },
        location,
    )
    return ChimeraSettings(
        enabled=_boolean(mapping.get("enabled", True), f"{location}.enabled"),
        window=_positive_int(mapping.get("window", 90), f"{location}.window"),
        step=_positive_int(mapping.get("step", 45), f"{location}.step"),
        minimum_depth=(
            None
            if mapping.get("minimum_depth") is None
            else _positive_int(mapping["minimum_depth"], f"{location}.minimum_depth")
        ),
        write_pcr_origin=_boolean(
            mapping.get("write_pcr_origin", False), f"{location}.write_pcr_origin"
        ),
        exclude_reads=_nonempty_string(
            mapping.get("exclude_reads", "written"), f"{location}.exclude_reads"
        ),
    )


def _plate_lists(value: Any, location: str) -> dict[str, tuple[str, ...]]:
    """Accept either a single plate name or a list of them, per key."""

    mapping = _mapping(value, location)
    result: dict[str, tuple[str, ...]] = {}
    for key, raw in mapping.items():
        name = _nonempty_string(key, f"{location} key")
        if isinstance(raw, str):
            items = [raw]
        elif isinstance(raw, Sequence):
            items = list(raw)
        else:
            raise ConfigError(f"{location}[{name!r}] must be a name or list of names")
        plates = tuple(
            _nonempty_string(item, f"{location}[{name!r}][{index}]")
            for index, item in enumerate(items)
        )
        if len(set(plates)) != len(plates):
            raise ConfigError(f"{location}[{name!r}] lists a plate more than once")
        result[name] = plates
    return result


def _parse_compressed_pcr(value: Any, base_dir: Path) -> CompressedPcrSettings:
    if value is None:
        return CompressedPcrSettings()
    location = "compressed_pcr"
    mapping = _mapping(value, location)
    _reject_unknown(
        mapping,
        {"enabled", "pcr_plates", "blocks", "block_pattern", "clonality", "layout_csv"},
        location,
    )
    enabled = _boolean(mapping.get("enabled", True), f"{location}.enabled")
    # A layout table and inline YAML are two sources of truth for the same design.
    # Refuse both rather than pick one: silently preferring either would make a
    # stale block map invisible.
    if "layout_csv" in mapping:
        clashing = sorted(
            key for key in ("pcr_plates", "blocks", "clonality") if key in mapping
        )
        if clashing:
            raise ConfigError(
                f"{location}: layout_csv replaces {', '.join(clashing)}; "
                "give one or the other, not both"
            )
    blocks_value = _mapping(mapping.get("blocks", {}), f"{location}.blocks")
    blocks = {
        _nonempty_string(library, f"{location}.blocks library key"): _plate_lists(
            by_block, f"{location}.blocks[{library!r}]"
        )
        for library, by_block in blocks_value.items()
    }
    clonality_value = _mapping(mapping.get("clonality", {}), f"{location}.clonality")
    clonality = {
        _nonempty_string(plate, f"{location}.clonality key"): _nonempty_string(
            mode, f"{location}.clonality[{plate!r}]"
        )
        for plate, mode in clonality_value.items()
    }
    pcr_plates = _plate_lists(mapping.get("pcr_plates", {}), f"{location}.pcr_plates")
    layout_csv = mapping.get("layout_csv")
    if layout_csv is not None:
        from .layout import LayoutError, read_layout_csv

        path = _path(layout_csv, base_dir, f"{location}.layout_csv")
        try:
            layout = read_layout_csv(path)
        except LayoutError as exc:
            raise ConfigError(str(exc)) from exc
        pcr_plates = layout["pcr_plates"]
        blocks = layout["blocks"]
        clonality = layout["clonality"]
    return CompressedPcrSettings(
        enabled=enabled,
        pcr_plates=pcr_plates,
        blocks=blocks,
        clonality=clonality,
        block_pattern=_nonempty_string(
            mapping.get("block_pattern", r"^Block_(\d+)_"), f"{location}.block_pattern"
        ),
    )


def available_presets() -> list[str]:
    """Names of the tuning presets shipped with the package."""

    from importlib.resources import files

    directory = files("nanopore3").joinpath("presets")
    return sorted(
        item.name[: -len(".yaml")]
        for item in directory.iterdir()
        if item.name.endswith(".yaml")
    )


def _deep_merge(base: Mapping[str, Any], override: Mapping[str, Any]) -> dict[str, Any]:
    """Merge ``override`` onto ``base``; a mapping merges, anything else replaces.

    A list replaces rather than extends.  Appending to an inherited list would make
    a preset's contents depend on the order two files were written in, and there is
    no way to remove an inherited item.
    """

    merged = dict(base)
    for key, value in override.items():
        if (
            key in merged
            and isinstance(merged[key], Mapping)
            and isinstance(value, Mapping)
        ):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def _apply_preset(root: Mapping[str, Any]) -> dict[str, Any]:
    """Resolve ``preset:`` into the configuration it stands for.

    The preset supplies tuning; the configuration supplies the experiment, and wins
    wherever the two mention the same key.  The result is what is validated,
    digested and recorded, so a preset is a way of writing a configuration, never a
    hidden layer under one - two runs naming the same preset are as reproducible as
    two runs spelling it out.
    """

    name = root.get("preset")
    if name is None:
        return dict(root)
    if not isinstance(name, str) or not name.strip():
        raise ConfigError("preset must be the name of a shipped preset")

    from importlib.resources import files

    source = files("nanopore3").joinpath("presets").joinpath(f"{name}.yaml")
    if not source.is_file():
        raise ConfigError(
            f"unknown preset {name!r}; available: {', '.join(available_presets())}"
        )
    loaded = yaml.safe_load(source.read_text(encoding="utf-8")) or {}
    merged = _deep_merge(_mapping(loaded, f"preset {name}"), root)
    del merged["preset"]
    # Kept so a reader of the run's recorded configuration can see where the
    # values came from without having to diff them against the preset.
    merged["run_name"] = merged.get("run_name")
    return merged


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
    root = _apply_preset(root)
    allowed = {
        "schema_version",
        "preset",
        "run_name",
        "output_root",
        "inputs",
        "references",
        "reference_libraries",
        "plate_reference_map",
        "library",
        "barcodes",
        "parallel",
        "consensus",
        "qc",
        "compressed_pcr",
        "chimera",
        "random_seed",
        "consensus_tree",
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
    has_references = "references" in root
    has_reference_libraries = "reference_libraries" in root
    if has_references == has_reference_libraries:
        raise ConfigError(
            "configuration must define exactly one of references or "
            "reference_libraries"
        )
    reference_libraries = (
        _parse_reference_libraries(root["reference_libraries"], base_dir)
        if has_reference_libraries
        else {}
    )
    plate_reference_map = (
        _parse_plate_reference_map(root.get("plate_reference_map", {}))
        if has_reference_libraries
        else {}
    )
    if not has_reference_libraries and "plate_reference_map" in root:
        raise ConfigError(
            "plate_reference_map is only valid with reference_libraries"
        )
    return PipelineConfig(
        schema_version=schema_version,
        run_name=_nonempty_string(root.get("run_name", source_path.stem), "run_name"),
        output_root=_path(
            _required(root, "output_root", "configuration"),
            base_dir,
            "output_root",
        ),
        inputs=_parse_inputs(_required(root, "inputs", "configuration"), base_dir),
        references=(
            _parse_references(root["references"], base_dir)
            if has_references
            else None
        ),
        reference_libraries=reference_libraries,
        plate_reference_map=plate_reference_map,
        library=_parse_library(root.get("library"), base_dir),
        plate_barcodes=_parse_barcodes(
            barcodes.get("plate"), "barcodes.plate", base_dir
        ),
        well_barcodes=_parse_barcodes(
            barcodes.get("well"), "barcodes.well", base_dir
        ),
        parallel=_parse_parallel(root.get("parallel")),
        consensus=_parse_consensus(root.get("consensus")),
        qc=_parse_qc(root.get("qc")),
        compressed_pcr=_parse_compressed_pcr(root.get("compressed_pcr"), base_dir),
        consensus_tree=_boolean(root.get("consensus_tree", True), "consensus_tree"),
        chimera=_parse_chimera(root.get("chimera")),
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
