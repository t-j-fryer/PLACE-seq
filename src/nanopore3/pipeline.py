"""Safe, portable orchestration for the Nanopore3 workflow."""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from contextlib import contextmanager
import csv
import logging
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
import gzip
import hashlib
import heapq
import io
import json
from multiprocessing import get_context
import re
from pathlib import Path
from typing import Any, Callable, Iterable, Iterator, Mapping, Sequence, TextIO, TypeVar

from . import __version__
from .assignment import (
    AssignmentCall,
    ReferenceIndex,
    _align_geometry,
    assign_read,
    extract_insert,
)
from .config import BarcodeSettings, PipelineConfig
from .consensus import (
    ConsensusRead,
    build_reference_consensus,
    require_mafft_spoa,
)
from .deconvolution import (
    NOT_CONFIGURED,
    CompressedPcrPlan,
    Deconvolution,
    block_from_reference_id,
)
from .chimera import (
    ASSEMBLY_ORIGIN,
    group_chimeras,
    signature_for,
    synthesise_reference,
)
from .export import normalized_well, safe_name, write_consensus_tree
from .demux import (
    BarcodeCall,
    PreparedBarcodePanel,
    call_barcode,
    prepare_barcode_panel,
    validate_barcodes,
)
from .io import FastqRecord, iter_fastq
from .provenance import (
    StageDirectory,
    atomic_write_json,
    canonical_digest,
    compute_stage_fingerprint,
    git_provenance,
    sha256_file,
)
from .fragments import FragmentLibrary
from .qc import (
    classify_mixed_positions,
    designed_allele_fraction,
    evaluate_consensus,
    mixed_signature,
    worst_protein_effect,
)
from .qc import _locate as _locate_motif
from .flanks import Flanks, from_assembled, from_sequences, from_template
from .references import read_fasta, read_reference_libraries
from .report import write_html_report
from .runtime import doctor_report, plan_resources
from .sequence import reverse_complement


LOGGER = logging.getLogger(__name__)

T = TypeVar("T")
U = TypeVar("U")

# Assignment costs milliseconds per read, so a small batch already dwarfs the
# per-future overhead while keeping the tail of the run evenly balanced across
# workers.  Demultiplexing is far cheaper per read and uses config.chunk_reads.
_ASSIGNMENT_BATCH_READS = 32

# Tail of the 3' constant region used to locate the terminal stop in a
# full-length consensus. Long enough to be unique, short enough that a couple of
# consensus errors cannot hide it.
CODING_STOP_ANCHOR = 24
# Same tolerance the coding QC allows when locating an anchor.
CODING_ANCHOR_EDITS = 3


class PipelineError(RuntimeError):
    """A user-actionable workflow failure."""


def _trimmed_barcodes(settings: BarcodeSettings) -> dict[str, str]:
    if settings.trim_bases == 0:
        return dict(settings.sequences)
    return {
        name: sequence[settings.trim_bases : -settings.trim_bases]
        for name, sequence in settings.sequences.items()
    }


def resolve_flanks(config: PipelineConfig) -> dict[str, Flanks]:
    """Resolve the constant flanks declared by each reference library.

    Explicit sequences are used as given.  A template is read and the insert it
    contains is found, which needs the inserts themselves, so the insert FASTA is
    parsed once plainly before it is parsed again flanked.  Parsing a thousand
    short records twice costs nothing next to being unable to check the template.
    """

    resolved: dict[str, Flanks] = {}
    for library_id, settings in config.reference_sets.items():
        if not settings.full_length:
            continue
        if settings.flanks_derive:
            records = read_fasta(settings.fasta).records
            flanks = from_assembled(
                [record.sequence for record in records],
                anchor_length=settings.flanks_anchor_length,
            )
            LOGGER.info(
                "library %s: flanks derived from %d assembled reference(s), "
                "%d constant bases per reference",
                library_id,
                len(records),
                flanks.constant_bases,
            )
            resolved[library_id] = flanks
            continue
        if settings.flanks_template is not None:
            inserts = [record.sequence for record in read_fasta(settings.fasta).records]
            template = read_fasta(settings.flanks_template).records
            if len(template) != 1:
                raise PipelineError(
                    f"reference library {library_id!r}: flanks.template must hold exactly "
                    f"one example construct, found {len(template)}"
                )
            flanks, matched = from_template(
                template[0].sequence,
                inserts,
                anchor_length=settings.flanks_anchor_length,
            )
            LOGGER.info(
                "library %s: flanks derived from %s (matched a %d nt insert), "
                "%d constant bases per reference",
                library_id,
                template[0].id,
                len(matched),
                flanks.constant_bases,
            )
        else:
            flanks = from_sequences(
                settings.flanks_upstream or "",
                settings.flanks_downstream or "",
                anchor_length=settings.flanks_anchor_length,
            )
        resolved[library_id] = flanks
    return resolved


def apply_flanks(config: PipelineConfig, flanks: Mapping[str, Flanks]) -> PipelineConfig:
    """Return a config whose full-length libraries bound the whole amplicon.

    The region a read is trimmed to is defined by the motif pair, so in
    full-length mode the motifs become the primer anchors at the outer ends of
    the constant regions.  An explicitly configured motif still wins, so an
    unusual construct can override.
    """

    if not flanks:
        return config
    libraries = dict(config.reference_sets)
    for library_id, flank in flanks.items():
        settings = libraries[library_id]
        libraries[library_id] = replace(
            settings,
            forward_motif=settings.forward_motif or flank.left_anchor,
            reverse_motif=settings.reverse_motif or flank.right_anchor,
        )
    if config.references is not None:
        return replace(config, references=libraries["default"])
    return replace(config, reference_libraries=libraries)


def flank_transforms(flanks: Mapping[str, Flanks]):
    """Sequence transforms that join the constant regions onto each insert."""

    return {
        library_id: (lambda _identifier, sequence, flank=flank: flank.transform(sequence))
        for library_id, flank in flanks.items()
    }


@dataclass(frozen=True)
class InsertGate:
    """Insert-scoped acceptance floors for one full-length library.

    A full-length reference is ~70% constant flank, so identity and coverage
    measured over the whole amplicon are roughly 3.5x less sensitive to anything
    wrong with the insert than the same floors were when the reference *was* the
    insert.  A read carrying 250 nt of a different design cleared 0.80 identity
    and 0.70 coverage comfortably, and went on to build a designed consensus that
    reported a perfect match - because a reference-guided consensus represents
    only the part of a molecule that aligns.

    This re-applies the library's own floors where they discriminate: to the
    insert.  Plain data, so it survives pickling to a process worker.
    """

    left_anchor: str
    right_anchor: str
    inner_left_anchor: str
    inner_right_anchor: str
    head: int
    tail: int
    motif_max_edits: int
    minimum_identity: float
    minimum_query_coverage: float
    inserts: Mapping[str, str]

    def insert_of(self, sequence: str) -> str | None:
        """The read's insert: by its own boundaries, else by offset."""

        exact = extract_insert(
            sequence, self.inner_left_anchor, self.inner_right_anchor,
            max_edits=self.motif_max_edits,
        )
        if exact.status == "found" and exact.sequence:
            return exact.sequence
        found = extract_insert(
            sequence, self.left_anchor, self.right_anchor, max_edits=self.motif_max_edits
        )
        if found.status != "found" or not found.sequence:
            return None
        region = found.sequence
        if len(region) <= self.head + self.tail:
            return None
        return region[self.head : len(region) - self.tail]

    def verdict(self, sequence: str, reference_id: str) -> tuple[float, float, bool] | None:
        """Insert identity, insert query coverage, and whether it passes."""

        reference = self.inserts.get(reference_id)
        insert = self.insert_of(sequence)
        if not reference or not insert:
            return None
        geometry = _align_geometry(insert, reference)
        passed = (
            geometry.identity >= self.minimum_identity
            and geometry.query_coverage >= self.minimum_query_coverage
        )
        return geometry.identity, geometry.query_coverage, passed


def build_insert_gates(
    config: PipelineConfig,
    references_by_library: Mapping[str, Mapping[str, str]],
) -> dict[str, InsertGate]:
    """One gate per full-length library that asked for insert-scoped floors."""

    gates: dict[str, InsertGate] = {}
    flanks = resolve_flanks(config)
    if not flanks:
        return gates
    inserts, _motifs = insert_view(references_by_library, flanks)
    for library_id, flank in flanks.items():
        settings = config.reference_sets[library_id]
        if not settings.insert_thresholds:
            continue
        inner_left, inner_right = flank.insert_anchors()
        gates[library_id] = InsertGate(
            left_anchor=flank.left_anchor,
            right_anchor=flank.right_anchor,
            inner_left_anchor=inner_left,
            inner_right_anchor=inner_right,
            head=len(flank.inner_upstream),
            tail=len(flank.inner_downstream),
            motif_max_edits=(
                config.library.motif_max_edits
                if settings.motif_max_edits is None
                else settings.motif_max_edits
            ),
            minimum_identity=settings.minimum_identity,
            minimum_query_coverage=settings.minimum_query_coverage,
            inserts=inserts.get(library_id, {}),
        )
    return gates


def insert_extractors(
    config: PipelineConfig,
    flanks: Mapping[str, Flanks],
) -> dict[str, Callable[[str], str | None]]:
    """Per-library "give me the insert of this read" functions.

    A full-length library looks for the insert boundaries first, and only if that
    search fails slices the insert positionally out of the primer-anchored region
    the assignment stage already matched.  The order matters and was measured:
    slicing positionally for *every* read gives fuzzier insert ends, which
    fragments signatures and cost 11 of 90 chimeric clones to the depth
    threshold.  Falling back only when needed is strictly additive - it recovers
    reads that were assignable but not profilable, like the six in RP07 A12 that
    escaped chimera detection and went on to build a spurious designed consensus,
    without moving any read that already worked.
    """

    extractors: dict[str, Callable[[str], str | None]] = {}
    for library_id, settings in config.reference_sets.items():
        max_edits = (
            config.library.motif_max_edits
            if settings.motif_max_edits is None
            else settings.motif_max_edits
        )
        flank = flanks.get(library_id)
        if flank is None:
            left = settings.forward_motif or config.library.forward_motif
            right = settings.reverse_motif or config.library.reverse_motif

            def extract(sequence: str, left=left, right=right, edits=max_edits):
                found = extract_insert(sequence, left, right, max_edits=edits)
                return found.sequence if found.status == "found" else None

        else:
            inner_left, inner_right = flank.insert_anchors()
            outer_left, outer_right = flank.left_anchor, flank.right_anchor
            head, tail = len(flank.inner_upstream), len(flank.inner_downstream)

            def extract(
                sequence: str,
                inner_left=inner_left,
                inner_right=inner_right,
                outer_left=outer_left,
                outer_right=outer_right,
                edits=max_edits,
                head=head,
                tail=tail,
            ):
                exact = extract_insert(sequence, inner_left, inner_right, max_edits=edits)
                if exact.status == "found" and exact.sequence:
                    return exact.sequence
                # Fallback: the read reached both primer sites even though its
                # insert boundary did not match, so take the insert by offset.
                found = extract_insert(sequence, outer_left, outer_right, max_edits=edits)
                if found.status != "found" or not found.sequence:
                    return None
                region = found.sequence
                if len(region) <= head + tail:
                    return None
                return region[head : len(region) - tail]

        extractors[library_id] = extract
    return extractors


def insert_view(
    references_by_library: Mapping[str, Mapping[str, str]],
    flanks: Mapping[str, Flanks],
) -> tuple[dict[str, dict[str, str]], dict[str, tuple[str, str]]]:
    """Insert-only references and motifs, for analyses that need discrimination.

    Positional profiling compares windows of a read against the reference set, so
    it is only informative where the references differ.  With full-length
    references, 849 of ~1,200 bases are shared by all of them and every window
    there matches everything.  Slicing the insert back out - by the same span the
    flanks define - restores the insert-only view without re-reading any file.
    """

    references: dict[str, dict[str, str]] = {}
    motifs: dict[str, tuple[str, str]] = {}
    for library_id, sequences in references_by_library.items():
        flank = flanks.get(library_id)
        if flank is None:
            references[library_id] = dict(sequences)
            continue
        sliced: dict[str, str] = {}
        for name, sequence in sequences.items():
            start, end = flank.insert_span(len(sequence))
            sliced[name] = sequence[start:end]
        references[library_id] = sliced
        motifs[library_id] = flank.insert_anchors()
    return references, motifs


def qc_regions(flank: Flanks, reference_length: int) -> dict[str, tuple[int, int]]:
    """Named spans for per-region accuracy in a full-length reference."""

    start, end = flank.insert_span(reference_length)
    spans: dict[str, tuple[int, int]] = {}
    if start:
        spans["flank_5p"] = (0, start)
    spans["insert"] = (start, end)
    if end < reference_length:
        spans["flank_3p"] = (end, reference_length)
    return spans


def check_pooling_layout(config: PipelineConfig, reference_collection) -> list[str]:
    """Cross-check the supplied pooling layout against the references and barcodes.

    The layout is experimental design and is never inferred from the reads - but it
    is transcribed by hand, and a block missing from it produces reads that simply
    cannot be attributed to a culture plate, three stages later and without saying
    why.  Every mismatch that can be named now is named now.
    """

    settings = config.compressed_pcr
    if not settings.enabled:
        return []
    pattern = re.compile(settings.block_pattern)
    problems: list[str] = []

    for library_id in sorted(config.reference_sets):
        try:
            records = reference_collection.get(library_id).records
        except Exception:  # a library the collection could not load reports elsewhere
            continue
        seen = set()
        for record in records:
            match = pattern.search(record.id)
            if match:
                seen.add(match.group(1) if match.groups() else match.group(0))
        if not seen:
            continue
        mapped = set(settings.blocks.get(library_id, {}))
        missing = sorted(seen - mapped)
        if missing:
            problems.append(
                f"library {library_id!r}: {len(missing)} block(s) in the references have "
                f"no culture plate: {', '.join(missing[:8])}"
                + (" ..." if len(missing) > 8 else "")
            )
        unused = sorted(mapped - seen)
        if unused:
            problems.append(
                f"library {library_id!r}: {len(unused)} block(s) in the layout match no "
                f"reference: {', '.join(unused[:8])}"
                + (" ..." if len(unused) > 8 else "")
            )

    pooled = {plate for plates in settings.pcr_plates.values() for plate in plates}
    placed = {
        plate
        for by_block in settings.blocks.values()
        for plates in by_block.values()
        for plate in plates
    }
    orphan = sorted(placed - pooled)
    if orphan:
        problems.append(
            f"{len(orphan)} culture plate(s) hold blocks but were not pooled into any "
            f"barcode: {', '.join(orphan[:8])}" + (" ..." if len(orphan) > 8 else "")
        )
    empty = sorted(pooled - placed)
    if empty:
        problems.append(
            f"{len(empty)} culture plate(s) were pooled but hold no block: "
            f"{', '.join(empty[:8])}" + (" ..." if len(empty) > 8 else "")
        )
    unknown = sorted(set(settings.pcr_plates) - set(config.plate_reference_map))
    if unknown:
        problems.append(
            f"{len(unknown)} barcode(s) in the layout are not in plate_reference_map: "
            f"{', '.join(unknown[:8])}"
        )
    return problems


def parse_mixed_alleles(
    value: str,
) -> tuple[tuple[int, str, str, tuple[float, ...]], ...]:
    """Read back the consensus stage's ``12:G>GT:0.61,0.39`` column.

    The fractions are optional, so a run written before they existed still parses.
    """

    out: list[tuple[int, str, str, tuple[float, ...]]] = []
    for item in filter(None, (value or "").split("|")):
        index, _, rest = item.partition(":")
        spec, _, shares = rest.partition(":")
        base, _, alleles = spec.partition(">")
        if not (index.isdigit() and base and alleles):
            continue
        fractions: tuple[float, ...] = ()
        try:
            fractions = tuple(float(f) for f in shares.split(",") if f)
        except ValueError:
            fractions = ()
        out.append((int(index), base, alleles, fractions))
    return tuple(out)


def reading_frame_span(
    reference: str, upstream: str, downstream: str
) -> tuple[int, int] | None:
    """Locate the ORF inside a full-length reference, as a half-open span.

    The same anchors the coding QC uses: the upstream constant begins at the start
    codon and the downstream one ends at the terminal stop, so the frame is
    measured rather than assumed.
    """

    start = _locate_motif(upstream, reference, CODING_ANCHOR_EDITS)
    stop = _locate_motif(downstream[-CODING_STOP_ANCHOR:], reference, CODING_ANCHOR_EDITS)
    if start is None or stop is None:
        return None
    return start[0], stop[1]


def validate_inputs(
    config: PipelineConfig,
    *,
    scan_fastq: bool = True,
    require_inputs: bool = True,
) -> dict[str, Any]:
    """Perform complete preflight without creating a run directory.

    ``require_inputs`` is False only for a rerun that inherits every stage which
    reads the sequencing data.  References are still required either way: the
    stages being recomputed compare against them.
    """

    missing = [] if not require_inputs else [
        str(item.path) for item in config.inputs if not item.path.is_file()
    ]
    missing.extend(
        str(path)
        for settings in config.reference_sets.values()
        for path in settings.fasta
        if not path.is_file()
    )
    if missing:
        raise PipelineError("missing input file(s): " + ", ".join(missing))
    flanks = resolve_flanks(config)
    config = apply_flanks(config, flanks)
    reference_collection = read_reference_libraries(
        {
            library_id: settings.fasta
            for library_id, settings in config.reference_sets.items()
        },
        transforms=flank_transforms(flanks),
    )
    if config.consensus.backend == "mafft_spoa":
        require_mafft_spoa()
    layout_problems = check_pooling_layout(config, reference_collection)
    if layout_problems:
        raise PipelineError(
            "compressed_pcr layout does not match the references:\n  "
            + "\n  ".join(layout_problems)
        )
    for label, settings in (
        ("plate", config.plate_barcodes),
        ("well", config.well_barcodes),
    ):
        if settings.sequences and settings.decision_policy == "best_margin":
            try:
                validate_barcodes(
                    _trimmed_barcodes(settings),
                    max_edits=settings.max_edits,
                    min_margin=settings.minimum_margin,
                    orientation_aware=settings.allow_reverse_complement,
                )
            except ValueError as exc:
                raise PipelineError(f"invalid {label} barcode panel: {exc}") from exc
    inputs: list[dict[str, Any]] = []
    for item in config.inputs:
        if not require_inputs and not item.path.is_file():
            continue
        digest = sha256_file(item.path)
        count = None
        if scan_fastq:
            count = sum(
                1
                for _ in iter_fastq(
                    item.path,
                    source_sha256=digest,
                    allow_empty_sequence=config.library.allow_empty_reads,
                )
            )
        inputs.append(
            {
                "sample_id": item.sample_id,
                "path": str(item.path),
                "sha256": digest,
                "size_bytes": item.path.stat().st_size,
                "records": count,
            }
        )
    return {
        "inputs": inputs,
        "references": sum(
            len(bundle.records) for _, bundle in reference_collection.libraries
        ),
        "reference_digest": reference_collection.digest,
        "reference_libraries": {
            library_id: {
                "references": len(bundle.records),
                "digest": bundle.digest,
                "alias_groups": [list(group) for group in bundle.alias_groups],
            }
            for library_id, bundle in reference_collection.libraries
        },
        "resources": asdict(
            plan_resources(config.parallel.jobs, config.parallel.threads_per_job)
        ),
    }


@contextmanager
def _open_gzip_text(path: Path) -> Iterator[TextIO]:
    """Open a deterministic gzip text writer and close every wrapper explicitly."""

    with path.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                yield text


def _iter_gzip_json(path: Path) -> Iterator[dict[str, Any]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            try:
                value = json.loads(line)
            except json.JSONDecodeError as exc:
                raise PipelineError(f"invalid JSONL at {path}:{line_number}: {exc}") from exc
            if not isinstance(value, dict):
                raise PipelineError(f"JSONL value at {path}:{line_number} is not an object")
            yield value


def _ordered_map(
    function: Callable[[T], U],
    values: Iterable[T],
    jobs: int,
    *,
    backend: str = "thread",
    initializer: Callable[..., None] | None = None,
    initargs: tuple[Any, ...] = (),
) -> Iterator[U]:
    """Map ``function`` over ``values`` in input order, bounded by ``jobs``.

    ``initializer`` runs once in each process worker, which lets a stage build
    expensive per-worker state such as reference indexes without pickling it
    with every task.  It is ignored by the serial path, which already has that
    state in the parent.
    """

    if jobs == 1:
        yield from map(function, values)
        return
    if backend not in {"thread", "process"}:
        raise ValueError("parallel map backend must be thread or process")
    iterator = iter(values)
    if backend == "process":
        pool_context = ProcessPoolExecutor(
            max_workers=jobs,
            mp_context=get_context("spawn"),
            initializer=initializer,
            initargs=initargs,
        )
    else:
        pool_context = ThreadPoolExecutor(
            max_workers=jobs,
            thread_name_prefix="nanopore3",
        )
    with pool_context as pool:
        pending = deque()
        # One pending unit per worker bounds memory when units are coarse batches.
        for _ in range(jobs):
            try:
                pending.append(pool.submit(function, next(iterator)))
            except StopIteration:
                break
        while pending:
            yield pending.popleft().result()
            try:
                pending.append(pool.submit(function, next(iterator)))
            except StopIteration:
                pass


def _batched(values: Iterable[T], size: int) -> Iterator[tuple[T, ...]]:
    """Yield bounded immutable work batches without reading the input eagerly."""

    if size < 1:
        raise ValueError("batch size must be positive")
    batch: list[T] = []
    for value in values:
        batch.append(value)
        if len(batch) == size:
            yield tuple(batch)
            batch = []
    if batch:
        yield tuple(batch)


def _call_dict(call: BarcodeCall, prefix: str) -> dict[str, Any]:
    return {
        f"{prefix}_status": call.status,
        f"{prefix}_id": call.barcode_id or "",
        f"{prefix}_orientation": call.orientation,
        f"{prefix}_best_edits": "" if call.best_distance is None else call.best_distance,
        f"{prefix}_second_edits": "" if call.second_distance is None else call.second_distance,
        f"{prefix}_margin": "" if call.margin is None else call.margin,
        f"{prefix}_matched_end": call.matched_end or "",
        f"{prefix}_start": "" if call.start is None else call.start,
        f"{prefix}_end": "" if call.end is None else call.end,
        f"{prefix}_reason": call.reason,
    }


def _disabled_barcode(identifier: str) -> BarcodeCall:
    return BarcodeCall(
        "assigned", identifier, "unknown", 0, None, None, None, None, None,
        "barcode stage disabled", (),
    )


def _barcode_call(
    sequence: str,
    settings: BarcodeSettings,
    disabled_id: str,
    prepared_panel: PreparedBarcodePanel | None = None,
) -> BarcodeCall:
    if not settings.sequences:
        return _disabled_barcode(disabled_id)
    return call_barcode(
        sequence,
        _trimmed_barcodes(settings),
        window_size=settings.search_window,
        max_edits=settings.max_edits,
        min_margin=settings.minimum_margin,
        search_ends=settings.search_ends,
        allow_reverse_complement=settings.allow_reverse_complement,
        decision_policy=settings.decision_policy,
        prepared_panel=prepared_panel,
    )


def _demux_one(
    job: tuple[
        FastqRecord,
        str,
        PipelineConfig,
        PreparedBarcodePanel | None,
        PreparedBarcodePanel | None,
    ]
) -> tuple[dict[str, Any], dict[str, Any] | None]:
    record, sample_id, config, plate_panel, well_panel = job
    if not record.sequence:
        # A zero-length basecall carries no evidence. It is reported as its own
        # state so it is counted rather than silently dropped or mistaken for a
        # barcode failure.
        row: dict[str, Any] = {
            "read_uid": record.read_uid,
            "original_read_id": record.name,
            "sample_id": sample_id,
            "record_index": record.record_index,
            "length": 0,
            "mean_q": "0.0000",
            "length_status": "empty_read",
            "quality_status": "not_evaluable",
            "call_status": "empty_read",
            "reason_code": "the basecaller emitted a zero-length read",
        }
        row.update(_call_dict(_disabled_barcode("not_attempted"), "plate"))
        row.update(_call_dict(_disabled_barcode("not_attempted"), "well"))
        return row, None
    length_status = "pass"
    if config.library.minimum_read_length is not None and len(record.sequence) < config.library.minimum_read_length:
        length_status = "out_of_length"
    if config.library.maximum_read_length is not None and len(record.sequence) > config.library.maximum_read_length:
        length_status = "out_of_length"
    quality_status = "pass"
    if (
        config.library.minimum_mean_quality is not None
        and record.mean_quality < config.library.minimum_mean_quality
    ):
        quality_status = "low_quality"

    if length_status != "pass" or quality_status != "pass":
        # These read-level gates already override any barcode verdict below, so
        # matching barcodes against a fragment, a dimer, or an unreadable read
        # only costs time. Rejecting here is what makes a length filter cheap.
        row = {
            "read_uid": record.read_uid,
            "original_read_id": record.name,
            "sample_id": sample_id,
            "record_index": record.record_index,
            "length": len(record.sequence),
            "mean_q": f"{record.mean_quality:.4f}",
            "length_status": length_status,
            "quality_status": quality_status,
            "call_status": length_status if length_status != "pass" else quality_status,
            "reason_code": (
                "read length outside configured range"
                if length_status != "pass"
                else "mean read quality below configured minimum"
            ),
        }
        row.update(_call_dict(_disabled_barcode("not_attempted"), "plate"))
        row.update(_call_dict(_disabled_barcode("not_attempted"), "well"))
        return row, None

    plate = _barcode_call(
        record.sequence, config.plate_barcodes, sample_id, plate_panel
    )
    sequence, quality = record.sequence, record.quality
    if plate.status == "assigned" and plate.orientation == "reverse":
        sequence, quality = reverse_complement(sequence), quality[::-1]
    well = (
        _barcode_call(sequence, config.well_barcodes, "well", well_panel)
        if plate.status == "assigned"
        else _disabled_barcode("not_attempted")
    )
    final_status = "assigned"
    reason = "plate and well calls passed"
    if length_status != "pass":
        final_status, reason = length_status, "read length outside configured range"
    elif quality_status != "pass":
        final_status, reason = quality_status, "mean read quality below configured minimum"
    elif plate.status != "assigned":
        final_status, reason = f"plate_{plate.status}", plate.reason
    elif well.status != "assigned":
        final_status, reason = f"well_{well.status}", well.reason
    elif (
        config.plate_barcodes.sequences
        and config.well_barcodes.sequences
        and well.orientation == "reverse"
    ):
        final_status = "orientation_conflict"
        reason = "well barcode orientation conflicts after plate-based orientation"
    elif not config.plate_barcodes.sequences and well.orientation == "reverse":
        sequence, quality = reverse_complement(sequence), quality[::-1]
    row: dict[str, Any] = {
        "read_uid": record.read_uid,
        "original_read_id": record.name,
        "sample_id": sample_id,
        "record_index": record.record_index,
        "length": len(record.sequence),
        "mean_q": f"{record.mean_quality:.4f}",
        "length_status": length_status,
        "quality_status": quality_status,
        "call_status": final_status,
        "reason_code": reason,
    }
    row.update(_call_dict(plate, "plate"))
    row.update(_call_dict(well, "well"))
    if final_status != "assigned":
        return row, None
    accepted = {
        "read_uid": record.read_uid,
        "original_read_id": record.name,
        "sample_id": sample_id,
        "plate_id": plate.barcode_id,
        "well_id": well.barcode_id,
        "sequence": sequence,
        "quality": quality,
    }
    return row, accepted


def _demux_batch(
    jobs: Sequence[
        tuple[
            FastqRecord,
            str,
            PipelineConfig,
            PreparedBarcodePanel | None,
            PreparedBarcodePanel | None,
        ]
    ],
) -> tuple[tuple[dict[str, Any], dict[str, Any] | None], ...]:
    """Process one coarse deterministic batch, amortizing executor overhead."""

    return tuple(_demux_one(job) for job in jobs)


def _assignment_row(
    read: Mapping[str, Any],
    call: AssignmentCall,
    k: int,
    reference_library_id: str,
    deconvolution: Deconvolution = NOT_CONFIGURED,
) -> dict[str, Any]:
    best, second = call.best, call.second
    return {
        "read_uid": read["read_uid"],
        "original_read_id": read["original_read_id"],
        "sample_id": read["sample_id"],
        "plate_id": read["plate_id"],
        "well_id": read["well_id"],
        "reference_library_id": reference_library_id,
        "assignment_status": call.status,
        "reference_ids": "|".join(call.reference_ids),
        "orientation": call.orientation,
        "motif_status": call.extraction.status if call.extraction else "not_configured",
        "kmer_size": k,
        "best_identity": "" if best is None else f"{best.identity:.6f}",
        "second_identity": "" if second is None else f"{second.identity:.6f}",
        "identity_margin": "" if call.identity_margin is None else f"{call.identity_margin:.6f}",
        "query_coverage": "" if best is None else f"{best.query_coverage:.6f}",
        "reference_coverage": "" if best is None else f"{best.reference_coverage:.6f}",
        "edit_distance": "" if best is None else best.edit_distance,
        "reason_code": call.reason,
        "assembly_block": deconvolution.block or "",
        "culture_plate": deconvolution.culture_plate or "",
        "culture_plate_status": deconvolution.status,
    }


def _assign_one(
    read: Mapping[str, Any],
    indexes_by_library: Mapping[str, Sequence[ReferenceIndex]],
    config: PipelineConfig,
    plan: CompressedPcrPlan | None = None,
    blocks_by_library: Mapping[str, Mapping[str, str]] | None = None,
    gates: Mapping[str, InsertGate] | None = None,
) -> tuple[dict[str, Any], dict[str, Any] | None]:
    try:
        library_id = config.reference_library_id_for_plate(str(read["plate_id"]))
    except KeyError:
        return (
            {
                "read_uid": read["read_uid"],
                "original_read_id": read["original_read_id"],
                "sample_id": read["sample_id"],
                "plate_id": read["plate_id"],
                "well_id": read["well_id"],
                "reference_library_id": "",
                "assignment_status": "unmapped_reference_library",
                "reference_ids": "",
                "orientation": "unknown",
                "motif_status": "not_attempted",
                "kmer_size": "",
                "best_identity": "",
                "second_identity": "",
                "identity_margin": "",
                "query_coverage": "",
                "reference_coverage": "",
                "edit_distance": "",
                "reason_code": (
                    f"plate barcode {read['plate_id']!r} has no configured "
                    "reference library"
                ),
                "assembly_block": "",
                "culture_plate": "",
                "culture_plate_status": "not_configured",
                "insert_identity": "",
                "insert_query_coverage": "",
            },
            None,
        )
    indexes = indexes_by_library[library_id]
    settings = config.reference_sets[library_id]
    final, used_k = assign_read(
        str(read["sequence"]),
        indexes,
        left_motif=settings.forward_motif or config.library.forward_motif,
        right_motif=settings.reverse_motif or config.library.reverse_motif,
        motif_max_edits=(
            config.library.motif_max_edits
            if settings.motif_max_edits is None
            else settings.motif_max_edits
        ),
        top_n=settings.candidate_count,
        min_kmer_score=settings.minimum_kmer_score,
        min_identity=settings.minimum_identity,
        min_query_coverage=settings.minimum_query_coverage,
        min_reference_coverage=settings.minimum_reference_coverage,
        min_identity_margin=settings.minimum_identity_margin,
        rescue=settings.rescue_policy,
        rescue_candidates=settings.rescue_candidates,
    )
    deconvolution = _deconvolve(
        read, final.reference_ids, library_id, plan, blocks_by_library or {}
    )
    row = _assignment_row(read, final, used_k, library_id, deconvolution)
    row.setdefault("insert_identity", "")
    row.setdefault("insert_query_coverage", "")
    # Insert-scoped floors. The amplicon-wide numbers can pass while the insert
    # is largely foreign, so a read that does not cover the design it was matched
    # to is held back from that design's consensus rather than diluting it.
    gate = (gates or {}).get(library_id)
    row["insert_identity"] = ""
    row["insert_query_coverage"] = ""
    if gate is not None and final.status in {"assigned_unique", "assigned_alias_set"}:
        verdict = gate.verdict(str(read["sequence"]), final.reference_ids[0])
        if verdict is not None:
            identity, coverage, passed = verdict
            row["insert_identity"] = f"{identity:.6f}"
            row["insert_query_coverage"] = f"{coverage:.6f}"
            if not passed:
                row["assignment_status"] = "insert_mismatch"
                row["reason_code"] = (
                    f"insert identity {identity:.3f} / coverage {coverage:.3f} "
                    f"below the library floor "
                    f"({gate.minimum_identity:.2f} / {gate.minimum_query_coverage:.2f}); "
                    "the read carries insert sequence this design does not explain"
                )
                return row, None
    eligible = None
    if final.status in {"assigned_unique", "assigned_alias_set"} and final.query_sequence:
        quality: str | None = str(read["quality"])
        if final.extraction is not None and final.extraction.left and final.extraction.right:
            if final.orientation == "reverse":
                quality = quality[::-1]
            quality = quality[
                final.extraction.left.end + 1 : final.extraction.right.start
            ]
        eligible = {
            **{
                key: read[key]
                for key in (
                    "read_uid",
                    "original_read_id",
                    "sample_id",
                    "plate_id",
                    "well_id",
                )
            },
            "quality": quality,
            "sequence": final.query_sequence,
            "reference_ids": list(final.reference_ids),
            "reference_library_id": library_id,
            "culture_plate": deconvolution.culture_plate or "",
        }
        if len(eligible["quality"]) != len(eligible["sequence"]):
            eligible["quality"] = None
    return row, eligible


def _build_block_map(
    config: PipelineConfig,
    references_by_library: Mapping[str, Mapping[str, str]],
) -> dict[str, dict[str, str]]:
    """Map each reference identifier to its assembly block, per library.

    A configured design table is authoritative and is joined by sequence, since
    the FASTA and the design table use different identifiers.  Otherwise the
    block is read from the reference identifier itself.
    """

    if not config.compressed_pcr.enabled:
        return {}
    pattern = re.compile(config.compressed_pcr.block_pattern)
    blocks: dict[str, dict[str, str]] = {}
    for library_id, references in references_by_library.items():
        settings = config.reference_sets[library_id]
        mapping: dict[str, str] = {}
        design = (
            FragmentLibrary.from_full_info_csv(settings.fragments_csv)
            if settings.fragments_csv is not None
            else None
        )
        for reference_id, sequence in references.items():
            block: str | None = None
            if design is not None:
                gene = design.gene_for_sequence(sequence)
                if gene is not None and gene.block:
                    block = gene.block
            if block is None:
                block = block_from_reference_id(reference_id, pattern)
            if block is not None:
                mapping[reference_id] = block
        blocks[library_id] = mapping
    return blocks


def _deconvolve(
    read: Mapping[str, Any],
    call_reference_ids: Sequence[str],
    library_id: str,
    plan: CompressedPcrPlan | None,
    blocks_by_library: Mapping[str, Mapping[str, str]],
) -> Deconvolution:
    """Resolve the source culture plate of one assigned read."""

    if plan is None:
        return NOT_CONFIGURED
    if not call_reference_ids:
        return Deconvolution(
            None, None, "unknown_block", "the read was not assigned to a reference"
        )
    blocks = blocks_by_library.get(library_id, {})
    observed = {blocks.get(reference_id) for reference_id in call_reference_ids}
    observed.discard(None)
    if len(observed) != 1:
        # An alias set spanning several blocks cannot name one culture plate.
        return Deconvolution(
            None, None, "unknown_block",
            "the assigned reference(s) do not identify exactly one block"
            if observed
            else "no assembly block is known for the assigned reference(s)",
        )
    return plan.resolve(str(read["plate_id"]), library_id, observed.pop())


def _build_reference_indexes(
    config: PipelineConfig,
    references_by_library: Mapping[str, Mapping[str, str]],
) -> dict[str, tuple[ReferenceIndex, ...]]:
    """Build one k-mer index per declared k-mer size, per reference library."""

    return {
        library_id: tuple(
            ReferenceIndex(
                references_by_library[library_id],
                k=k,
                max_kmer_owners=config.reference_sets[library_id].max_kmer_owners,
            )
            for k in config.reference_sets[library_id].kmer_sizes
        )
        for library_id in references_by_library
    }


# Assignment is dominated by Python-level k-mer work rather than by the
# GIL-releasing edlib calls, so process workers are the only way it scales.
# Each worker builds its own indexes once instead of receiving them per task.
_ASSIGNMENT_WORKER: dict[str, Any] = {}


def _init_assignment_worker(
    config: PipelineConfig,
    references_by_library: Mapping[str, Mapping[str, str]],
) -> None:
    """Prepare one process worker. Safe under the spawn start method."""

    _ASSIGNMENT_WORKER["config"] = config
    _ASSIGNMENT_WORKER["indexes"] = _build_reference_indexes(config, references_by_library)
    _ASSIGNMENT_WORKER["gates"] = build_insert_gates(config, references_by_library)
    _ASSIGNMENT_WORKER["blocks"] = _build_block_map(config, references_by_library)
    _ASSIGNMENT_WORKER["plan"] = (
        CompressedPcrPlan(
            config.compressed_pcr.pcr_plates,
            config.compressed_pcr.blocks,
            config.compressed_pcr.clonality,
        )
        if config.compressed_pcr.enabled
        else None
    )


def _assign_batch_worker(
    batch: Sequence[Mapping[str, Any]],
) -> tuple[tuple[dict[str, Any], dict[str, Any] | None], ...]:
    """Assign one coarse batch inside a process worker."""

    config = _ASSIGNMENT_WORKER["config"]
    indexes = _ASSIGNMENT_WORKER["indexes"]
    plan = _ASSIGNMENT_WORKER["plan"]
    blocks = _ASSIGNMENT_WORKER["blocks"]
    gates = _ASSIGNMENT_WORKER["gates"]
    return tuple(
        _assign_one(read, indexes, config, plan, blocks, gates) for read in batch
    )


def _detect_chimeras(
    config: PipelineConfig,
    demux_reads: Path,
    references_by_library: Mapping[str, Mapping[str, str]],
    root: Path,
    plan: CompressedPcrPlan | None = None,
    region_motifs: Mapping[str, tuple[str, str]] | None = None,
    extractors: Mapping[str, Callable[[str], str | None]] | None = None,
) -> dict[str, Any]:
    """Group reads by positional signature and write a consensus per clone.

    ``references_by_library`` and ``region_motifs`` scope the work: a full-length
    run passes the insert-only view, because a window of shared constant sequence
    matches every reference and carries no signal about a junction.
    """

    settings = config.chimera
    minimum_depth = settings.minimum_depth or config.consensus.minimum_depth
    block_pattern = (
        re.compile(config.compressed_pcr.block_pattern)
        if config.compressed_pcr.block_pattern
        else None
    )
    indexes = {
        library: ReferenceIndex(
            references,
            k=config.reference_sets[library].kmer_sizes[0],
            max_kmer_owners=config.reference_sets[library].max_kmer_owners,
        )
        for library, references in references_by_library.items()
    }
    wells: dict[tuple[str, str, str], list[tuple[str, str]]] = defaultdict(list)
    for read in _iter_gzip_json(demux_reads):
        plate = str(read["plate_id"])
        try:
            library = config.reference_library_id_for_plate(plate)
        except KeyError:
            continue
        wells[(plate, str(read["well_id"]), library)].append(
            (str(read["read_uid"]), str(read["sequence"]))
        )

    root.mkdir(parents=True, exist_ok=True)
    totals: Counter[str] = Counter()
    rows: list[dict[str, Any]] = []
    scaffolds: list[str] = []
    # read_uid -> the chimeric clone that accounts for it, so the consensus stage
    # can stop the same molecule appearing a second time under a design's name.
    claimed: list[tuple[str, str]] = []
    for (plate, well, library), reads in sorted(wells.items()):
        reference_settings = config.reference_sets[library]
        index = indexes[library]
        signatures = []
        inserts: dict[str, str] = {}
        scoped = (region_motifs or {}).get(library)
        extract = (extractors or {}).get(library)
        if extract is None:
            left_motif, right_motif = scoped or (
                reference_settings.forward_motif or config.library.forward_motif,
                reference_settings.reverse_motif or config.library.reverse_motif,
            )
            max_edits = (
                config.library.motif_max_edits
                if reference_settings.motif_max_edits is None
                else reference_settings.motif_max_edits
            )

            def extract(sequence: str, left=left_motif, right=right_motif, edits=max_edits):
                found = extract_insert(sequence, left, right, max_edits=edits)
                return found.sequence if found.status == "found" else None

        for read_id, sequence in reads:
            insert = extract(sequence)
            if not insert:
                continue
            inserts[read_id] = insert
            signatures.append(
                signature_for(read_id, insert, index,
                              window=settings.window, step=settings.step)
            )
        if not signatures:
            continue
        groups, outcome = group_chimeras(
            signatures, minimum_depth=minimum_depth, block_pattern=block_pattern
        )
        totals.update(outcome)
        for group in groups:
            keep = group.origin == ASSEMBLY_ORIGIN or settings.write_pcr_origin
            if not keep and settings.exclude_reads == "all":
                # Not a clone, so nothing is written; but its reads are not
                # evidence for any design either.
                for read_id in group.read_ids:
                    claimed.append((read_id, f"pcr:{group.label}"))
            if not keep:
                # Counted above, but not written: a PCR template-switch product
                # is not a clone, and a FASTA beside real ones invites misreading.
                totals["pcr_origin_not_written"] += 1
                continue
            members = [
                ConsensusRead(rid, inserts[rid])
                for rid in group.read_ids[: config.consensus.maximum_reads]
                if rid in inserts
            ]
            scaffold = synthesise_reference(
                group.signature, references_by_library[library],
                group.junction_window, window=settings.window, step=settings.step,
            )
            if not scaffold or not members:
                totals["no_scaffold"] += 1
                continue
            result = build_reference_consensus(
                scaffold, members, group_id=group.label, min_depth=1,
                max_reads=len(members), min_support=config.consensus.minimum_support,
                seed=config.random_seed,
            )
            if not result.sequence:
                totals["no_consensus"] += 1
                continue
            directory = root / safe_name(plate) / normalized_well(well)
            directory.mkdir(parents=True, exist_ok=True)
            label = "__".join(safe_name(p, limit=40) for p in group.signature)
            name = f"{safe_name(plate)}_{normalized_well(well)}__{label}__chimera.fasta"
            header = (
                f">{safe_name(plate)}_{normalized_well(well)}_chimera "
                f"parents={'|'.join(group.signature)} origin={group.origin} "
                f"reads={group.size} junction_window={group.junction_window} "
                f"library={library} region={'insert' if scoped else 'amplicon'}"
            )
            wrapped = "\n".join(
                result.sequence[i : i + 60] for i in range(0, len(result.sequence), 60)
            )
            (directory / name).write_text(f"{header}\n{wrapped}\n", encoding="ascii")
            totals["written"] += 1
            # Both parents share a block for an assembly-origin chimera, so the
            # culture plate is well defined; resolve it so the clone files under
            # the same provenance as every other consensus from this well.
            culture_plate = ""
            if plan is not None and block_pattern is not None:
                match = block_pattern.search(group.signature[0])
                if match is not None:
                    resolved = plan.resolve(
                        plate, library,
                        match.group(1) if match.groups() else match.group(0),
                    )
                    culture_plate = resolved.culture_plate or ""
            chimera_id = f"chim-{canonical_digest({'g': group.label, 'p': plate, 'w': well})[:16]}"
            scaffolds.append(f">{chimera_id}\n{scaffold}\n")
            if settings.exclude_reads in ("written", "all"):
                claimed.extend((read_id, chimera_id) for read_id in group.read_ids)
            rows.append({
                "chimera_id": chimera_id, "culture_plate": culture_plate,
                "plate_id": plate, "well_id": well, "reference_library_id": library,
                "parents": " >> ".join(group.signature), "n_parents": len(group.signature),
                "origin": group.origin, "reads": group.size,
                # A full-length run detects chimeras on the insert, so the clone
                # written here is insert-scoped while its designed neighbours span
                # the whole amplicon. Say which, rather than leave it to be
                # inferred from a length.
                "region": "insert" if scoped else "amplicon",
                "junction_window": group.junction_window,
                "ambiguous_bases": result.ambiguous_bases,
                "length": len(result.sequence),
                "file": str((directory / name).relative_to(root)),
            })
    if rows:
        with (root / "clones.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    # The spliced scaffold is what QC must grade a chimera against; storing it
    # here keeps QC from having to re-derive a junction it did not compute.
    (root / "scaffolds.fasta").write_text("".join(scaffolds), encoding="ascii")
    with (root / "claimed_reads.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["read_uid", "chimera_id"])
        writer.writerows(sorted(claimed))
    totals["reads_claimed"] = len(claimed)
    return dict(sorted(totals.items()))


def _write_csv(path: Path, rows: Iterable[Mapping[str, Any]], fields: Sequence[str]) -> Counter[str]:
    counts: Counter[str] = Counter()
    with _open_gzip_text(path) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
            status = row.get("call_status", row.get("assignment_status", "rows"))
            counts[str(status)] += 1
    return counts


def _stage_fingerprint(
    stage: str,
    parameters: Mapping[str, Any],
    inputs: Mapping[str, str],
    *,
    backend_versions: Mapping[str, str] | None = None,
) -> str:
    return compute_stage_fingerprint(
        stage,
        pipeline_version=__version__,
        parameters=parameters,
        input_digests=inputs,
        backend_versions=backend_versions,
    )


def run_pipeline(
    config: PipelineConfig,
    *,
    output_root: Path | None = None,
    run_id: str | None = None,
    resume: bool = False,
    inherit: Mapping[str, Any] | None = None,
) -> Path:
    """Execute the portable workflow into a new immutable run directory.

    ``inherit`` is the record returned by :func:`nanopore3.rerun.prepare_rerun`:
    earlier stages already placed in the run directory, carried from a finished
    run.  When it covers the stages that read the FASTQ, the FASTQ itself is not
    required - those stages will not run, and their manifests already carry its
    digest.
    """

    from .rerun import INPUT_CONSUMING_STAGES

    inherited_stages = set((inherit or {}).get("inherited_stages", ()))
    skip_inputs = bool(inherit) and all(
        stage in inherited_stages for stage in INPUT_CONSUMING_STAGES
    )
    if skip_inputs:
        # Provenance still names the input and its checksum: they come from the
        # inherited run rather than from re-reading a file nobody will open.
        preflight = validate_inputs(config, scan_fastq=False, require_inputs=False)
        preflight["inputs"] = list(inherit["source_preflight"].get("inputs", []))
        preflight["inherited_inputs"] = True
    else:
        preflight = validate_inputs(config, scan_fastq=True)
    # The digest covers the configuration as written, which is what a reader
    # declared; the derived anchors and flank digests are recorded in the
    # assignment stage manifest instead.
    config_digest = canonical_digest(config.as_dict())
    flanks = resolve_flanks(config)
    config = apply_flanks(config, flanks)
    if run_id is None:
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        run_id = f"{stamp}-{config_digest[:12]}"
    root = (output_root or config.output_root).expanduser().resolve(strict=False)
    run_dir = root / run_id
    if run_dir.exists() and not resume:
        raise PipelineError(f"run already exists; choose another run ID or use --resume: {run_dir}")
    run_dir.mkdir(parents=True, exist_ok=True)
    if not resume and any(run_dir.iterdir()):
        raise PipelineError(f"new run directory is not empty: {run_dir}")
    runtime_report = doctor_report()
    edlib_version = str(runtime_report["packages"].get("edlib") or "unknown")
    edlib_backend_versions = {"edlib": edlib_version}
    consensus_backend_versions = dict(edlib_backend_versions)
    if config.consensus.backend == "mafft_spoa":
        optional = runtime_report["optional_binaries"]
        for name in ("mafft", "spoa"):
            details = optional[name]
            consensus_backend_versions[name] = str(details["version"] or "unknown")
    run_metadata = {
        "schema_version": 1,
        "pipeline_version": __version__,
        "run_id": run_id,
        "config_digest": config_digest,
        "config": config.as_dict(),
        "preflight": preflight,
        "runtime": runtime_report,
        "source_control": git_provenance(Path(__file__).resolve()),
    }
    if inherit:
        run_metadata["inherited"] = dict(inherit)
    metadata_path = run_dir / "run.json"
    if metadata_path.exists():
        existing = json.loads(metadata_path.read_text(encoding="utf-8"))
        if existing.get("config_digest") != config_digest and not inherit:
            raise PipelineError(
                "resume refused: run configuration differs. To reuse this run's "
                "intermediates under a changed configuration, use `nanopore3 rerun`, "
                "which writes a new run and re-derives every stage the change reaches"
            )
    else:
        atomic_write_json(metadata_path, run_metadata)

    input_digests = {item["sample_id"]: item["sha256"] for item in preflight["inputs"]}
    requested_jobs = 1 if config.parallel.backend == "serial" else config.parallel.jobs
    resources = plan_resources(requested_jobs, config.parallel.threads_per_job)

    ingest_parameters = {"schema": 1}
    ingest_fp = _stage_fingerprint("01_ingest", ingest_parameters, input_digests)
    with StageDirectory(
        run_dir,
        "01_ingest",
        ingest_fp,
        pipeline_version=__version__,
        parameters=ingest_parameters,
        input_digests=input_digests,
        resume=resume,
    ) as stage:
        if not stage.reused:
            atomic_write_json(stage.output_path("inputs.json"), preflight)

    demux_parameters = {
        "plate": asdict(config.plate_barcodes),
        "well": asdict(config.well_barcodes),
        "library": asdict(config.library),
    }
    demux_fp = _stage_fingerprint(
        "02_demux",
        demux_parameters,
        input_digests,
        backend_versions=edlib_backend_versions,
    )
    with StageDirectory(
        run_dir,
        "02_demux",
        demux_fp,
        pipeline_version=__version__,
        parameters=demux_parameters,
        input_digests=input_digests,
        backend_versions=edlib_backend_versions,
        resume=resume,
    ) as stage:
        if not stage.reused:
            calls_path = stage.output_path("demux_calls.csv.gz")
            reads_path = stage.output_path("demuxed_reads.jsonl.gz")
            fields: list[str] | None = None
            counts: Counter[str] = Counter()
            plate_panel = (
                prepare_barcode_panel(_trimmed_barcodes(config.plate_barcodes))
                if config.plate_barcodes.sequences
                else None
            )
            well_panel = (
                prepare_barcode_panel(_trimmed_barcodes(config.well_barcodes))
                if config.well_barcodes.sequences
                else None
            )
            with _open_gzip_text(calls_path) as call_handle, _open_gzip_text(reads_path) as read_handle:
                writer: csv.DictWriter[str] | None = None
                for item in config.inputs:
                    digest = input_digests[item.sample_id]
                    jobs = (
                        (record, item.sample_id, config, plate_panel, well_panel)
                        for record in iter_fastq(
                            item.path,
                            source_sha256=digest,
                            allow_empty_sequence=config.library.allow_empty_reads,
                        )
                    )
                    batches = _batched(jobs, config.parallel.chunk_reads)
                    for results in _ordered_map(
                        _demux_batch,
                        batches,
                        resources.jobs,
                        backend=(
                            "process"
                            if config.parallel.backend == "process"
                            else "thread"
                        ),
                    ):
                        for row, accepted in results:
                            if writer is None:
                                fields = list(row)
                                writer = csv.DictWriter(
                                    call_handle,
                                    fieldnames=fields,
                                    lineterminator="\n",
                                )
                                writer.writeheader()
                            writer.writerow(row)
                            counts[row["call_status"]] += 1
                            if accepted is not None:
                                read_handle.write(
                                    json.dumps(
                                        accepted,
                                        sort_keys=True,
                                        separators=(",", ":"),
                                    )
                                    + "\n"
                                )
            atomic_write_json(stage.output_path("summary.json"), dict(sorted(counts.items())))

    demux_dir = run_dir / "stages" / "02_demux"
    demux_reads = demux_dir / "demuxed_reads.jsonl.gz"
    reference_collection = read_reference_libraries(
        {
            library_id: settings.fasta
            for library_id, settings in config.reference_sets.items()
        },
        transforms=flank_transforms(flanks),
    )
    references_by_library = {
        library_id: {record.id: record.sequence for record in bundle.records}
        for library_id, bundle in reference_collection.libraries
    }
    indexes_by_library = _build_reference_indexes(config, references_by_library)
    insert_gates = build_insert_gates(config, references_by_library)
    if insert_gates:
        LOGGER.info(
            "insert-scoped acceptance floors active for: %s",
            ", ".join(sorted(insert_gates)),
        )
    block_map = _build_block_map(config, references_by_library)
    compressed_plan = (
        CompressedPcrPlan(
            config.compressed_pcr.pcr_plates,
            config.compressed_pcr.blocks,
            config.compressed_pcr.clonality,
        )
        if config.compressed_pcr.enabled
        else None
    )
    reference_parameters = {
        library_id: asdict(settings)
        for library_id, settings in config.reference_sets.items()
    }
    assignment_inputs = {
        "demuxed_reads": sha256_file(demux_reads),
        "references": reference_collection.digest,
    }
    assignment_parameters = {
        "insert_thresholds": sorted(insert_gates),
        "reference_libraries": reference_parameters,
        "plate_reference_map": dict(config.plate_reference_map),
        "compressed_pcr": asdict(config.compressed_pcr),
    }
    assign_fp = _stage_fingerprint(
        "03_assignment",
        assignment_parameters,
        assignment_inputs,
        backend_versions=edlib_backend_versions,
    )
    with StageDirectory(
        run_dir,
        "03_assignment",
        assign_fp,
        pipeline_version=__version__,
        parameters=assignment_parameters,
        input_digests=assignment_inputs,
        backend_versions=edlib_backend_versions,
        resume=resume,
    ) as stage:
        if not stage.reused:
            calls_path = stage.output_path("assignment_calls.csv.gz")
            eligible_path = stage.output_path("consensus_eligible.jsonl.gz")
            counts: Counter[str] = Counter()
            with _open_gzip_text(calls_path) as call_handle, _open_gzip_text(eligible_path) as eligible_handle:
                writer: csv.DictWriter[str] | None = None
                batches = _batched(_iter_gzip_json(demux_reads), _ASSIGNMENT_BATCH_READS)

                def assign_batch(
                    batch: Sequence[Mapping[str, Any]],
                ) -> tuple[tuple[dict[str, Any], dict[str, Any] | None], ...]:
                    return tuple(
                        _assign_one(
                            read, indexes_by_library, config, compressed_plan,
                            block_map, insert_gates,
                        )
                        for read in batch
                    )

                # Threads cannot speed this stage up: its cost is Python-level
                # k-mer work, not the GIL-releasing alignment calls. Only the
                # explicitly requested process backend runs workers in parallel.
                use_processes = config.parallel.backend == "process"
                for results in _ordered_map(
                    _assign_batch_worker if use_processes else assign_batch,
                    batches,
                    resources.jobs if use_processes else 1,
                    backend="process" if use_processes else "thread",
                    initializer=_init_assignment_worker if use_processes else None,
                    initargs=(config, references_by_library) if use_processes else (),
                ):
                    for row, eligible in results:
                        if writer is None:
                            writer = csv.DictWriter(call_handle, fieldnames=list(row), lineterminator="\n")
                            writer.writeheader()
                        writer.writerow(row)
                        counts[row["assignment_status"]] += 1
                        if eligible is not None:
                            eligible_handle.write(json.dumps(eligible, sort_keys=True, separators=(",", ":")) + "\n")
            atomic_write_json(stage.output_path("summary.json"), dict(sorted(counts.items())))

    assignment_dir = run_dir / "stages" / "03_assignment"
    # A chimeric molecule is a real clone in a real well, and assignment cannot
    # describe it: an even split falls below the identity floor, a lopsided one
    # is attributed to whichever parent dominates. This stage asks a positional
    # question instead and gives those reads a consensus of their own. It sits
    # beside consensus rather than inside it so that changing a chimera setting
    # does not invalidate the reference-guided consensuses.
    if config.chimera.enabled:
        chimera_parameters = {**asdict(config.chimera), "seed": config.random_seed}
        chimera_inputs = {
            "demuxed_reads": sha256_file(demux_reads),
            "references": reference_collection.digest,
        }
        chimera_fp = _stage_fingerprint(
            "03b_chimera", chimera_parameters, chimera_inputs,
            backend_versions=edlib_backend_versions,
        )
        with StageDirectory(
            run_dir, "03b_chimera", chimera_fp,
            pipeline_version=__version__, parameters=chimera_parameters,
            input_digests=chimera_inputs, backend_versions=edlib_backend_versions,
            resume=resume,
        ) as stage:
            if not stage.reused:
                chimera_references, chimera_motifs = insert_view(
                    references_by_library, flanks
                )
                summary = _detect_chimeras(
                    config, demux_reads, chimera_references,
                    stage.output_path("clones"), compressed_plan,
                    region_motifs=chimera_motifs,
                    extractors=insert_extractors(config, flanks),
                )
                atomic_write_json(stage.output_path("summary.json"), summary)

    eligible_path = assignment_dir / "consensus_eligible.jsonl.gz"
    # A read that belongs to a chimeric clone is already building that clone's
    # consensus. Letting it also vote for one of the parent designs makes one
    # molecule appear twice, the second time under a name it does not have.
    claimed_path = run_dir / "stages" / "03b_chimera" / "clones" / "claimed_reads.csv"
    claimed_reads: set[str] = set()
    if config.chimera.enabled and config.chimera.exclude_reads != "none":
        if claimed_path.is_file():
            with claimed_path.open("r", encoding="utf-8", newline="") as handle:
                claimed_reads = {row["read_uid"] for row in csv.DictReader(handle)}
        LOGGER.info(
            "excluding %d read(s) already accounted for by a chimeric clone",
            len(claimed_reads),
        )
    consensus_inputs = {
        "eligible": sha256_file(eligible_path),
        "references": reference_collection.digest,
    }
    if claimed_path.is_file():
        consensus_inputs["chimera_claimed_reads"] = sha256_file(claimed_path)
    consensus_parameters = {
        **asdict(config.consensus),
        "seed": config.random_seed,
        "exclude_chimeric_reads": config.chimera.exclude_reads
        if config.chimera.enabled
        else "none",
    }
    consensus_fp = _stage_fingerprint(
        "04_consensus",
        consensus_parameters,
        consensus_inputs,
        backend_versions=consensus_backend_versions,
    )
    with StageDirectory(
        run_dir,
        "04_consensus",
        consensus_fp,
        pipeline_version=__version__,
        parameters=consensus_parameters,
        input_digests=consensus_inputs,
        backend_versions=consensus_backend_versions,
        resume=resume,
    ) as stage:
        if not stage.reused:
            groups: dict[
                tuple[str, str, str, str, tuple[str, ...]],
                list[tuple[int, str, ConsensusRead]],
            ] = defaultdict(list)
            available: Counter[
                tuple[str, str, str, str, tuple[str, ...]]
            ] = Counter()
            # Provenance only: the grouping key is deliberately unchanged so
            # consensus identities stay stable whether or not deconvolution runs.
            culture_plates: dict[
                tuple[str, str, str, str, tuple[str, ...]], set[str]
            ] = defaultdict(set)
            excluded = 0
            for row in _iter_gzip_json(eligible_path):
                if row["read_uid"] in claimed_reads:
                    # Counted here rather than silently skipped: a stage that
                    # drops input must say how much.
                    excluded += 1
                    continue
                aliases = tuple(row["reference_ids"])
                key = (
                    row["sample_id"],
                    row["plate_id"],
                    row["well_id"],
                    row["reference_library_id"],
                    aliases,
                )
                available[key] += 1
                if row.get("culture_plate"):
                    culture_plates[key].add(str(row["culture_plate"]))
                rank = int(hashlib.sha256(f"{config.random_seed}\0{key}\0{row['read_uid']}".encode()).hexdigest(), 16)
                qualities = None if row["quality"] is None else tuple(ord(c) - 33 for c in row["quality"])
                item = (-rank, row["read_uid"], ConsensusRead(row["read_uid"], row["sequence"], qualities))
                heap = groups[key]
                if len(heap) < config.consensus.maximum_reads:
                    heapq.heappush(heap, item)
                elif item > heap[0]:
                    heapq.heapreplace(heap, item)
            consensus_rows: list[dict[str, Any]] = []
            contributor_rows: list[dict[str, Any]] = []
            fasta_parts: list[str] = []
            for key in sorted(groups):
                sample, plate, well, library_id, aliases = key
                reads = tuple(item[2] for item in groups[key])
                ref = references_by_library[library_id][aliases[0]]
                group_id = "|".join((sample, plate, well, library_id, *aliases))
                result = build_reference_consensus(
                        ref,
                        reads,
                        group_id=group_id,
                        min_depth=config.consensus.minimum_depth,
                        # 0 means every eligible read: the cap used to discard 84%
                        # of them, and the discarded depth is what a second allele
                        # shows up in.
                        max_reads=config.consensus.maximum_reads or len(reads),
                        min_support=config.consensus.minimum_support,
                        min_base_quality=config.consensus.minimum_base_quality,
                        significance=config.consensus.significance,
                        minor_fraction=config.consensus.minimum_minor_fraction,
                        seed=config.random_seed,
                        backend=config.consensus.backend,
                    )
                consensus_id = "cons-" + canonical_digest({"group": group_id, "sequence": result.sequence})[:16]
                row = {
                    "consensus_id": consensus_id, "sample_id": sample, "plate_id": plate, "well_id": well,
                    "reference_library_id": library_id,
                    "reference_ids": "|".join(aliases), "status": result.status,
                    "culture_plate": "|".join(sorted(culture_plates.get(key, ()))),
                    "n_reads_available": available[key], "n_reads_used": result.n_reads_used,
                    "mean_depth": f"{result.mean_depth:.4f}", "min_depth": result.min_depth,
                    "ambiguous_bases": result.ambiguous_bases,
                    # Marginal calls used to leave no trace: a base decided on 61%
                    # support looked identical to one decided on 100%. These make
                    # the evidence behind a consensus auditable.
                    "mixed_positions": result.mixed_positions,
                    # Which alleles disagreed, not just how many places. An "N"
                    # cannot say whether the position carries a damage signature.
                    "mixed_alleles": "|".join(
                        f"{index}:{base}>{alleles}:"
                        + ",".join(f"{f:.4f}" for f in fractions)
                        for index, base, alleles, fractions in result.mixed_alleles
                    ),
                    "background_error_rate": f"{result.background_error_rate:.5f}",
                    "weakest_support": f"{result.weakest_support:.4f}",
                    "low_quality_bases": result.low_quality_bases,
                    "backend": result.backend,
                    "sequence_sha256": hashlib.sha256(result.sequence.encode()).hexdigest() if result.sequence else "",
                    "failure_reason": result.failure_reason or "",
                }
                consensus_rows.append(row)
                if result.sequence:
                    fasta_parts.append(
                        f">{consensus_id} reference_library={library_id} "
                        f"reference_ids={'|'.join(aliases)} sample={sample} "
                        f"plate={plate} well={well}"
                        + (
                            f" culture_plate={'|'.join(sorted(culture_plates.get(key, ())))}"
                            if culture_plates.get(key)
                            else ""
                        )
                        + f"\n{result.sequence}\n"
                    )
                for rank, read_uid in enumerate(result.contributor_ids, start=1):
                    contributor_rows.append({"consensus_id": consensus_id, "read_uid": read_uid, "selection_rank": rank})
            fields = list(consensus_rows[0]) if consensus_rows else ["consensus_id", "status"]
            _write_csv(stage.output_path("consensus.csv.gz"), consensus_rows, fields)
            _write_csv(stage.output_path("contributors.csv.gz"), contributor_rows, ["consensus_id", "read_uid", "selection_rank"])
            stage.output_path("consensus.fasta").write_text("".join(fasta_parts), encoding="ascii")
            consensus_summary = dict(sorted(Counter(row["status"] for row in consensus_rows).items()))
            if claimed_reads:
                consensus_summary["reads_excluded_as_chimeric"] = excluded
            atomic_write_json(stage.output_path("summary.json"), consensus_summary)


    consensus_dir = run_dir / "stages" / "04_consensus"
    qc_inputs = {
        "consensus": sha256_file(consensus_dir / "consensus.fasta"),
        "references": reference_collection.digest,
    }
    qc_parameters = asdict(config.qc)
    qc_fp = _stage_fingerprint(
        "05_qc",
        qc_parameters,
        qc_inputs,
        backend_versions=edlib_backend_versions,
    )
    with StageDirectory(
        run_dir,
        "05_qc",
        qc_fp,
        pipeline_version=__version__,
        parameters=qc_parameters,
        input_digests=qc_inputs,
        backend_versions=edlib_backend_versions,
        resume=resume,
    ) as stage:
        if not stage.reused:
            sequences: dict[str, str] = {}
            current = None
            for line in (consensus_dir / "consensus.fasta").read_text(encoding="ascii").splitlines():
                if line.startswith(">"):
                    current = line[1:].split()[0]; sequences[current] = ""
                elif current is not None:
                    sequences[current] += line.strip()
            qc_rows: list[dict[str, Any]] = []
            with gzip.open(consensus_dir / "consensus.csv.gz", "rt", encoding="utf-8", newline="") as handle:
                for row in csv.DictReader(handle):
                    sequence = sequences.get(row["consensus_id"])
                    if not sequence:
                        qc_rows.append(
                            {
                                "consensus_id": row["consensus_id"],
                                "reference_library_id": row["reference_library_id"],
                                "reference_ids": row["reference_ids"],
                                "full_amplicon": "not_evaluable",
                                "expected_length": "not_evaluable",
                                "reading_frame": "not_evaluable",
                                "internal_stops": "not_evaluable",
                                "overall": "not_evaluable",
                                "alignment_edit_distance": "",
                                "alignment_identity": "",
                                "alignment_query_coverage": "",
                                "alignment_reference_coverage": "",
                                "reason": "consensus sequence unavailable",
                                "protein_length": "",
                                "internal_stop_codon": "",
                            }
                        )
                        continue
                    aliases = tuple(row["reference_ids"].split("|"))
                    library_id = row["reference_library_id"]
                    reference_sequence = references_by_library[library_id][aliases[0]]
                    upstream = (
                        config.reference_sets[library_id].qc_upstream_constant
                        or config.qc.upstream_constant
                    )
                    downstream = (
                        config.reference_sets[library_id].qc_downstream_constant
                        or config.qc.downstream_constant
                    )
                    flank = flanks.get(library_id)
                    result = evaluate_consensus(
                        sequence,
                        reference_sequence,
                        min_identity=config.qc.minimum_identity,
                        min_query_coverage=config.qc.minimum_query_coverage,
                        min_reference_coverage=config.qc.minimum_reference_coverage,
                        length_tolerance=config.qc.length_tolerance,
                        # Full-length consensus: the ORF is inside it, so locate
                        # it rather than assembling constants around it.
                        upstream_constant=None if flank else upstream,
                        downstream_constant=None if flank else downstream,
                        coding_anchors=(
                            (upstream, downstream[-CODING_STOP_ANCHOR:])
                            if flank and upstream and downstream
                            else None
                        ),
                        spans=(
                            qc_regions(flank, len(reference_sequence)) if flank else None
                        ),
                    )
                    # Where each contested position sits, and whether it carries
                    # the G:C -> T:A signature of a pre-transformation lesion. Two
                    # alleles in the vector backbone and two in the reading frame
                    # mean very different things for whether the clone is usable.
                    mixed = classify_mixed_positions(
                        parse_mixed_alleles(row.get("mixed_alleles", "")),
                        reference_sequence,
                        spans=qc_regions(flank, len(reference_sequence)) if flank else None,
                        reading_frame=(
                            reading_frame_span(reference_sequence, upstream, downstream)
                            if flank and upstream and downstream
                            else None
                        ),
                    )
                    qc_rows.append(
                        {
                            "consensus_id": row["consensus_id"],
                            "reference_library_id": library_id,
                            "reference_ids": row["reference_ids"],
                            **result.to_dict(),
                            "mixed_signature": mixed_signature(mixed),
                            "mixed_detail": " ".join(p.describe() for p in mixed),
                            "mixed_in_reading_frame": sum(1 for p in mixed if p.protein_effect),
                            # A mixed clone still containing the designed base is
                            # a different screening decision from one that does
                            # not: the intended sequence is in the well.
                            "designed_allele_fraction": (
                                "" if designed_allele_fraction(mixed) is None
                                else f"{designed_allele_fraction(mixed):.4f}"
                            ),
                            "mixed_worst_effect": worst_protein_effect(mixed),
                        }
                    )
            fields = [
                "consensus_id",
                "reference_library_id",
                "reference_ids",
                "full_amplicon",
                "expected_length",
                "reading_frame",
                "internal_stops",
                "overall",
                "alignment_edit_distance",
                "alignment_identity",
                "alignment_query_coverage",
                "alignment_reference_coverage",
                "reason",
                "protein_length",
                "internal_stop_codon",
                # Empty for every clone with one allele everywhere, which is most
                # of them; present so a mixed clone can be judged without going
                # back to the reads.
                "mixed_signature",
                "mixed_detail",
                "mixed_in_reading_frame",
                "mixed_worst_effect",
                "designed_allele_fraction",
            ]
            # Full-length references report accuracy per region as well as over
            # the whole amplicon; the columns exist only when a library declares
            # flanks, so an insert-mode run keeps exactly its old schema.
            for region in ("flank_5p", "insert", "flank_3p"):
                if any(f"{region}_identity" in row for row in qc_rows):
                    fields.extend(
                        [
                            f"{region}_length",
                            f"{region}_edit_distance",
                            f"{region}_identity",
                        ]
                    )
            # Chimeric clones are sequences that are present in the well, so
            # they are graded and exported like any other consensus. They are
            # scored against the spliced parent scaffold, which is the thing a
            # correct chimeric clone should equal.
            chimera_dir = run_dir / "stages" / "03b_chimera" / "clones"
            chimera_rows: list[dict[str, Any]] = []
            if (chimera_dir / "clones.csv").exists():
                scaffolds: dict[str, str] = {}
                current = None
                scaffold_file = chimera_dir / "scaffolds.fasta"
                if scaffold_file.exists():
                    for line in scaffold_file.read_text(encoding="ascii").splitlines():
                        if line.startswith(">"):
                            current = line[1:].split()[0]; scaffolds[current] = ""
                        elif current is not None:
                            scaffolds[current] += line.strip()
                with (chimera_dir / "clones.csv").open(encoding="utf-8", newline="") as handle:
                    for row in csv.DictReader(handle):
                        path = chimera_dir / row["file"]
                        scaffold = scaffolds.get(row["chimera_id"], "")
                        if not path.is_file() or not scaffold:
                            continue
                        sequence = "".join(
                            line.strip()
                            for line in path.read_text(encoding="ascii").splitlines()
                            if not line.startswith(">")
                        )
                        library_id = row["reference_library_id"]
                        result = evaluate_consensus(
                            sequence, scaffold,
                            min_identity=config.qc.minimum_identity,
                            min_query_coverage=config.qc.minimum_query_coverage,
                            min_reference_coverage=config.qc.minimum_reference_coverage,
                            length_tolerance=config.qc.length_tolerance,
                            upstream_constant=(
                                config.reference_sets[library_id].qc_upstream_constant
                                or config.qc.upstream_constant
                            ),
                            downstream_constant=(
                                config.reference_sets[library_id].qc_downstream_constant
                                or config.qc.downstream_constant
                            ),
                        )
                        chimera_rows.append({
                            "consensus_id": row["chimera_id"],
                            "reference_library_id": library_id,
                            "reference_ids": row["parents"].replace(" >> ", "+"),
                            **result.to_dict(),
                        })
                qc_rows.extend(chimera_rows)
            _write_csv(stage.output_path("qc.csv.gz"), qc_rows, fields)
            atomic_write_json(stage.output_path("summary.json"), dict(sorted(Counter(row["overall"] for row in qc_rows).items())))
            # A browsable, graded copy of the consensuses. It lives here rather
            # than in 04_consensus because the grade needs QC, and a promoted
            # stage directory is immutable.
            with gzip.open(consensus_dir / "consensus.csv.gz", "rt", encoding="utf-8", newline="") as handle:
                consensus_rows = list(csv.DictReader(handle))
            if (chimera_dir / "clones.csv").exists():
                with (chimera_dir / "clones.csv").open(encoding="utf-8", newline="") as handle:
                    for row in csv.DictReader(handle):
                        path = chimera_dir / row["file"]
                        if not path.is_file():
                            continue
                        sequences[row["chimera_id"]] = "".join(
                            line.strip()
                            for line in path.read_text(encoding="ascii").splitlines()
                            if not line.startswith(">")
                        )
                        consensus_rows.append({
                            "consensus_id": row["chimera_id"],
                            "sample_id": "", "plate_id": row["plate_id"],
                            "well_id": row["well_id"],
                            "reference_library_id": row["reference_library_id"],
                            "reference_ids": row["parents"].replace(" >> ", "+"),
                            "status": "chimera",
                            "culture_plate": row["culture_plate"],
                            "n_reads_available": row["reads"], "n_reads_used": row["reads"],
                            "mean_depth": row["reads"], "min_depth": row["reads"],
                            "ambiguous_bases": row["ambiguous_bases"],
                            "mixed_positions": 0,
                            "background_error_rate": "",
                            "weakest_support": "",
                            "low_quality_bases": 0,
                            "backend": "portable", "sequence_sha256": "",
                            "failure_reason": "",
                        })
            tree_summary = write_consensus_tree(
                consensus_rows,
                {row["consensus_id"]: row for row in qc_rows},
                sequences,
                stage.output_path("consensus_by_plate"),
            )
            atomic_write_json(
                stage.output_path("consensus_by_plate_summary.json"), tree_summary
            )

    report_inputs = {
        "demux": sha256_file(demux_dir / "summary.json"),
        "assignment": sha256_file(assignment_dir / "summary.json"),
        "consensus": sha256_file(consensus_dir / "summary.json"),
        "qc": sha256_file(run_dir / "stages" / "05_qc" / "summary.json"),
    }
    report_parameters = {"format": "html-v1"}
    report_fp = _stage_fingerprint("06_report", report_parameters, report_inputs)
    with StageDirectory(
        run_dir,
        "06_report",
        report_fp,
        pipeline_version=__version__,
        parameters=report_parameters,
        input_digests=report_inputs,
        resume=resume,
    ) as stage:
        if not stage.reused:
            sections = {
                name: json.loads((path).read_text(encoding="utf-8"))
                for name, path in {
                    "Demultiplexing": demux_dir / "summary.json",
                    "Assignment": assignment_dir / "summary.json",
                    "Consensus": consensus_dir / "summary.json",
                    "QC": run_dir / "stages" / "05_qc" / "summary.json",
                }.items()
            }
            chimera_summary = run_dir / "stages" / "03b_chimera" / "summary.json"
            if chimera_summary.exists():
                sections["Chimeric clones"] = json.loads(
                    chimera_summary.read_text(encoding="utf-8")
                )
            if compressed_plan is not None:
                # Deconvolution is the point of compressing plates, so the
                # report says how much of each pool was actually recovered
                # rather than leaving it to the figures alone.
                from .figures import summarize_culture_plates

                recovery: dict[str, Any] = {}
                for item in summarize_culture_plates(
                    run_dir,
                    {plate: len(sources)
                     for plate, sources in compressed_plan.pcr_plates.items()},
                ):
                    per_well = sorted(item.sources_per_well.values())
                    recovery[item.plate_id] = {
                        "pooled_culture_plates": item.pooled,
                        "culture_plates_recovered": len(item.consensus_per_source),
                        "consensus_sequences": sum(item.consensus_per_source.values()),
                        "groups_below_depth": item.below_depth,
                        "median_source_plates_per_well": (
                            per_well[len(per_well) // 2] if per_well else 0
                        ),
                        "weakest_culture_plate": (
                            min(item.consensus_per_source,
                                key=item.consensus_per_source.get)
                            if item.consensus_per_source else ""
                        ),
                    }
                if recovery:
                    sections["Culture plate recovery"] = recovery
            write_html_report(stage.output_path("report.html"), title=f"Nanopore3 — {config.run_name}", sections=sections, provenance={"run_id": run_id, "pipeline_version": __version__, "config_digest": config_digest})
            # Figures need matplotlib, which is an optional extra. A run must
            # not fail because a plotting library is absent, so this is
            # best-effort and records why it was skipped.
            try:
                from .figures import write_all as _write_figures
            except ImportError as exc:
                atomic_write_json(
                    stage.output_path("figures.json"),
                    {"written": [], "skipped": f"matplotlib unavailable: {exc}"},
                )
            else:
                expected = {}
                if compressed_plan is not None:
                    for plate in compressed_plan.pcr_plates:
                        value = compressed_plan.expected_clones_per_well(plate)
                        if value is not None:
                            expected[plate] = value
                pooled = (
                    {plate: len(sources)
                     for plate, sources in compressed_plan.pcr_plates.items()}
                    if compressed_plan is not None
                    else {}
                )
                written = _write_figures(
                    run_dir, stage.output_path("figures"), expected, pooled
                )
                atomic_write_json(
                    stage.output_path("figures.json"),
                    {"written": sorted(p.name for p in written), "skipped": None},
                )

    if config.consensus_tree:
        # Built after the stages rather than inside one: it is derived entirely
        # from stage outputs, and a screening decision is made from it, so it
        # belongs where a person will find it rather than nested in 05_qc.
        from .export import write_tree_from_stages

        summary = write_tree_from_stages(
            run_dir / "stages" / "04_consensus",
            run_dir / "stages" / "05_qc" / "qc.csv.gz",
            run_dir / "consensus_by_plate",
        )
        LOGGER.info(
            "graded consensus tree: %d clone(s) in %s",
            summary["files_written"],
            run_dir / "consensus_by_plate",
        )
    return run_dir


__all__ = ["PipelineError", "run_pipeline", "validate_inputs"]
