#!/usr/bin/env python3
"""Benchmark reference assignment without writing scientific outputs.

The stage reads demultiplexed reads produced by ``02_demux`` so that the
measurement covers exactly the work the pipeline performs, including the
k-mer cascade and the rescue policy declared by each reference library.
"""

from __future__ import annotations

import argparse
import gzip
import json
from collections import Counter
from multiprocessing import freeze_support
from pathlib import Path
from time import perf_counter

from nanopore3.config import load_config
from nanopore3.pipeline import (
    _assign_batch_worker,
    _assign_one,
    _batched,
    _build_reference_indexes,
    _init_assignment_worker,
    _ordered_map,
)
from nanopore3.references import read_reference_libraries


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument(
        "--demuxed-reads",
        type=Path,
        required=True,
        help="demuxed_reads.jsonl.gz from a completed 02_demux stage",
    )
    parser.add_argument("--backend", choices=("serial", "thread", "process"), required=True)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--batch-reads", type=int, default=32)
    parser.add_argument(
        "--rescue",
        choices=("none", "kmer", "all"),
        help="override every library's rescue policy for this measurement",
    )
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()
    if args.jobs < 1 or args.batch_reads < 1:
        parser.error("--jobs and --batch-reads must be positive")

    config = load_config(args.config)
    if args.rescue is not None:
        config = config.with_rescue_policy(args.rescue)

    collection = read_reference_libraries(
        {library_id: settings.fasta for library_id, settings in config.reference_sets.items()}
    )
    references = {
        library_id: {record.id: record.sequence for record in bundle.records}
        for library_id, bundle in collection.libraries
    }
    index_started = perf_counter()
    indexes = _build_reference_indexes(config, references)
    index_seconds = perf_counter() - index_started

    def reads():
        with gzip.open(args.demuxed_reads, "rt", encoding="utf-8", newline="") as handle:
            for position, line in enumerate(handle):
                if args.limit is not None and position >= args.limit:
                    break
                yield json.loads(line)

    def assign_batch(batch):
        return tuple(_assign_one(read, indexes, config) for read in batch)

    worker_count = 1 if args.backend == "serial" else args.jobs
    use_processes = args.backend == "process"
    started = perf_counter()
    counts: Counter[str] = Counter()
    total = 0
    for results in _ordered_map(
        _assign_batch_worker if use_processes else assign_batch,
        _batched(reads(), args.batch_reads),
        worker_count,
        backend="process" if use_processes else "thread",
        initializer=_init_assignment_worker if use_processes else None,
        initargs=(config, references) if use_processes else (),
    ):
        total += len(results)
        counts.update(row["assignment_status"] for row, _ in results)
    elapsed = perf_counter() - started

    print(
        json.dumps(
            {
                "backend": args.backend,
                "jobs": worker_count,
                "batch_reads": args.batch_reads,
                "rescue": args.rescue or "per-library",
                "records": total,
                "index_seconds": round(index_seconds, 6),
                "seconds": round(elapsed, 6),
                "reads_per_second": round(total / elapsed, 3) if elapsed else None,
                "call_counts": dict(sorted(counts.items())),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    freeze_support()
    raise SystemExit(main())
