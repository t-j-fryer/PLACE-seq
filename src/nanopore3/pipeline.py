"""Safe, portable orchestration for the conservative Nanopore3 v0.1 workflow."""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
import csv
from dataclasses import asdict
from datetime import datetime, timezone
import gzip
import hashlib
import heapq
import io
import json
from pathlib import Path
from typing import Any, Callable, Iterable, Iterator, Mapping, Sequence, TextIO, TypeVar

from . import __version__
from .assignment import AssignmentCall, ReferenceIndex, assign_sequence
from .config import BarcodeSettings, PipelineConfig
from .consensus import (
    ConsensusRead,
    build_reference_consensus,
    require_mafft_spoa,
)
from .demux import BarcodeCall, call_barcode, validate_barcodes
from .io import FastqRecord, iter_fastq
from .provenance import (
    StageDirectory,
    atomic_write_json,
    canonical_digest,
    compute_stage_fingerprint,
    git_provenance,
    sha256_file,
)
from .qc import evaluate_consensus
from .references import read_reference_libraries
from .report import write_html_report
from .runtime import doctor_report, plan_resources
from .sequence import reverse_complement


T = TypeVar("T")
U = TypeVar("U")


class PipelineError(RuntimeError):
    """A user-actionable workflow failure."""


def _trimmed_barcodes(settings: BarcodeSettings) -> dict[str, str]:
    if settings.trim_bases == 0:
        return dict(settings.sequences)
    return {
        name: sequence[settings.trim_bases : -settings.trim_bases]
        for name, sequence in settings.sequences.items()
    }


def validate_inputs(config: PipelineConfig, *, scan_fastq: bool = True) -> dict[str, Any]:
    """Perform complete preflight without creating a run directory."""

    missing = [str(item.path) for item in config.inputs if not item.path.is_file()]
    missing.extend(
        str(path)
        for settings in config.reference_sets.values()
        for path in settings.fasta
        if not path.is_file()
    )
    if missing:
        raise PipelineError("missing input file(s): " + ", ".join(missing))
    reference_collection = read_reference_libraries(
        {
            library_id: settings.fasta
            for library_id, settings in config.reference_sets.items()
        }
    )
    if config.consensus.backend == "mafft_spoa":
        require_mafft_spoa()
    for label, settings in (
        ("plate", config.plate_barcodes),
        ("well", config.well_barcodes),
    ):
        if settings.sequences:
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
        digest = sha256_file(item.path)
        count = None
        if scan_fastq:
            count = sum(1 for _ in iter_fastq(item.path, source_sha256=digest))
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


def _ordered_map(function: Callable[[T], U], values: Iterable[T], jobs: int) -> Iterator[U]:
    if jobs == 1:
        yield from map(function, values)
        return
    iterator = iter(values)
    with ThreadPoolExecutor(max_workers=jobs, thread_name_prefix="nanopore3") as pool:
        pending = deque()
        for _ in range(jobs * 4):
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


def _barcode_call(sequence: str, settings: BarcodeSettings, disabled_id: str) -> BarcodeCall:
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
    )


def _demux_one(job: tuple[FastqRecord, str, PipelineConfig]) -> tuple[dict[str, Any], dict[str, Any] | None]:
    record, sample_id, config = job
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

    plate = _barcode_call(record.sequence, config.plate_barcodes, sample_id)
    sequence, quality = record.sequence, record.quality
    if plate.status == "assigned" and plate.orientation == "reverse":
        sequence, quality = reverse_complement(sequence), quality[::-1]
    well = (
        _barcode_call(sequence, config.well_barcodes, "well")
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


def _assignment_row(
    read: Mapping[str, Any],
    call: AssignmentCall,
    k: int,
    reference_library_id: str,
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
    }


def _assign_one(
    read: Mapping[str, Any],
    indexes_by_library: Mapping[str, Sequence[ReferenceIndex]],
    config: PipelineConfig,
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
            },
            None,
        )
    indexes = indexes_by_library[library_id]
    settings = config.reference_sets[library_id]
    final: AssignmentCall | None = None
    used_k = indexes[-1].k
    for index in indexes:
        call = assign_sequence(
            str(read["sequence"]),
            index,
            left_motif=config.library.forward_motif,
            right_motif=config.library.reverse_motif,
            motif_max_edits=config.library.motif_max_edits,
            top_n=settings.candidate_count,
            min_kmer_score=settings.minimum_kmer_score,
            min_identity=settings.minimum_identity,
            min_query_coverage=settings.minimum_query_coverage,
            min_reference_coverage=settings.minimum_reference_coverage,
            min_identity_margin=settings.minimum_identity_margin,
            fallback_align_all=True,
        )
        final, used_k = call, index.k
        if call.status in {"assigned_unique", "assigned_alias_set", "ambiguous", "motif_missing"}:
            break
    assert final is not None
    row = _assignment_row(read, final, used_k, library_id)
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
        }
        if len(eligible["quality"]) != len(eligible["sequence"]):
            eligible["quality"] = None
    return row, eligible


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
) -> Path:
    """Execute the v0.1 portable workflow into a new immutable run directory."""

    preflight = validate_inputs(config, scan_fastq=True)
    config_digest = canonical_digest(config.as_dict())
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
    metadata_path = run_dir / "run.json"
    if metadata_path.exists():
        existing = json.loads(metadata_path.read_text(encoding="utf-8"))
        if existing.get("config_digest") != config_digest:
            raise PipelineError("resume refused: run configuration differs")
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
            with _open_gzip_text(calls_path) as call_handle, _open_gzip_text(reads_path) as read_handle:
                writer: csv.DictWriter[str] | None = None
                for item in config.inputs:
                    digest = input_digests[item.sample_id]
                    jobs = ((record, item.sample_id, config) for record in iter_fastq(item.path, source_sha256=digest))
                    for row, accepted in _ordered_map(_demux_one, jobs, resources.jobs):
                        if writer is None:
                            fields = list(row)
                            writer = csv.DictWriter(call_handle, fieldnames=fields, lineterminator="\n")
                            writer.writeheader()
                        writer.writerow(row)
                        counts[row["call_status"]] += 1
                        if accepted is not None:
                            read_handle.write(json.dumps(accepted, sort_keys=True, separators=(",", ":")) + "\n")
            atomic_write_json(stage.output_path("summary.json"), dict(sorted(counts.items())))

    demux_dir = run_dir / "stages" / "02_demux"
    demux_reads = demux_dir / "demuxed_reads.jsonl.gz"
    reference_collection = read_reference_libraries(
        {
            library_id: settings.fasta
            for library_id, settings in config.reference_sets.items()
        }
    )
    references_by_library = {
        library_id: {record.id: record.sequence for record in bundle.records}
        for library_id, bundle in reference_collection.libraries
    }
    indexes_by_library = {
        library_id: tuple(
            ReferenceIndex(
                references_by_library[library_id],
                k=k,
                max_kmer_owners=config.reference_sets[library_id].max_kmer_owners,
            )
            for k in config.reference_sets[library_id].kmer_sizes
        )
        for library_id in reference_collection.ids
    }
    reference_parameters = {
        library_id: asdict(settings)
        for library_id, settings in config.reference_sets.items()
    }
    assignment_inputs = {
        "demuxed_reads": sha256_file(demux_reads),
        "references": reference_collection.digest,
    }
    assignment_parameters = {
        "reference_libraries": reference_parameters,
        "plate_reference_map": dict(config.plate_reference_map),
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
                jobs = _iter_gzip_json(demux_reads)
                fn = lambda read: _assign_one(read, indexes_by_library, config)
                for row, eligible in _ordered_map(fn, jobs, resources.jobs):
                    if writer is None:
                        writer = csv.DictWriter(call_handle, fieldnames=list(row), lineterminator="\n")
                        writer.writeheader()
                    writer.writerow(row)
                    counts[row["assignment_status"]] += 1
                    if eligible is not None:
                        eligible_handle.write(json.dumps(eligible, sort_keys=True, separators=(",", ":")) + "\n")
            atomic_write_json(stage.output_path("summary.json"), dict(sorted(counts.items())))

    assignment_dir = run_dir / "stages" / "03_assignment"
    eligible_path = assignment_dir / "consensus_eligible.jsonl.gz"
    consensus_inputs = {
        "eligible": sha256_file(eligible_path),
        "references": reference_collection.digest,
    }
    consensus_parameters = {**asdict(config.consensus), "seed": config.random_seed}
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
            for row in _iter_gzip_json(eligible_path):
                aliases = tuple(row["reference_ids"])
                key = (
                    row["sample_id"],
                    row["plate_id"],
                    row["well_id"],
                    row["reference_library_id"],
                    aliases,
                )
                available[key] += 1
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
                    max_reads=config.consensus.maximum_reads,
                    min_support=config.consensus.minimum_support,
                    seed=config.random_seed,
                    backend=config.consensus.backend,
                    threads=config.parallel.threads_per_job,
                )
                consensus_id = "cons-" + canonical_digest({"group": group_id, "sequence": result.sequence})[:16]
                row = {
                    "consensus_id": consensus_id, "sample_id": sample, "plate_id": plate, "well_id": well,
                    "reference_library_id": library_id,
                    "reference_ids": "|".join(aliases), "status": result.status,
                    "n_reads_available": available[key], "n_reads_used": result.n_reads_used,
                    "mean_depth": f"{result.mean_depth:.4f}", "min_depth": result.min_depth,
                    "ambiguous_bases": result.ambiguous_bases, "backend": result.backend,
                    "sequence_sha256": hashlib.sha256(result.sequence.encode()).hexdigest() if result.sequence else "",
                    "failure_reason": result.failure_reason or "",
                }
                consensus_rows.append(row)
                if result.sequence:
                    fasta_parts.append(
                        f">{consensus_id} reference_library={library_id} "
                        f"reference_ids={'|'.join(aliases)} sample={sample} "
                        f"plate={plate} well={well}\n{result.sequence}\n"
                    )
                for rank, read_uid in enumerate(result.contributor_ids, start=1):
                    contributor_rows.append({"consensus_id": consensus_id, "read_uid": read_uid, "selection_rank": rank})
            fields = list(consensus_rows[0]) if consensus_rows else ["consensus_id", "status"]
            _write_csv(stage.output_path("consensus.csv.gz"), consensus_rows, fields)
            _write_csv(stage.output_path("contributors.csv.gz"), contributor_rows, ["consensus_id", "read_uid", "selection_rank"])
            stage.output_path("consensus.fasta").write_text("".join(fasta_parts), encoding="ascii")
            atomic_write_json(stage.output_path("summary.json"), dict(sorted(Counter(row["status"] for row in consensus_rows).items())))

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
                            }
                        )
                        continue
                    aliases = tuple(row["reference_ids"].split("|"))
                    library_id = row["reference_library_id"]
                    result = evaluate_consensus(
                        sequence,
                        references_by_library[library_id][aliases[0]],
                        min_identity=config.qc.minimum_identity,
                        min_query_coverage=config.qc.minimum_query_coverage,
                        min_reference_coverage=config.qc.minimum_reference_coverage,
                        length_tolerance=config.qc.length_tolerance,
                    )
                    qc_rows.append(
                        {
                            "consensus_id": row["consensus_id"],
                            "reference_library_id": library_id,
                            "reference_ids": row["reference_ids"],
                            **result.to_dict(),
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
            ]
            _write_csv(stage.output_path("qc.csv.gz"), qc_rows, fields)
            atomic_write_json(stage.output_path("summary.json"), dict(sorted(Counter(row["overall"] for row in qc_rows).items())))

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
            write_html_report(stage.output_path("report.html"), title=f"Nanopore3 — {config.run_name}", sections=sections, provenance={"run_id": run_id, "pipeline_version": __version__, "config_digest": config_digest})
    return run_dir


__all__ = ["PipelineError", "run_pipeline", "validate_inputs"]
