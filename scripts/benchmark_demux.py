#!/usr/bin/env python3
"""Benchmark deterministic demultiplexing without writing scientific outputs."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from multiprocessing import freeze_support
from pathlib import Path
from time import perf_counter

from nanopore3.config import load_config
from nanopore3.demux import prepare_barcode_panel
from nanopore3.io import iter_fastq
from nanopore3.pipeline import (
    _batched,
    _demux_batch,
    _ordered_map,
    _trimmed_barcodes,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--fastq", type=Path, required=True)
    parser.add_argument("--backend", choices=("serial", "thread", "process"), required=True)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--chunk-reads", type=int, default=1000)
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()
    if args.jobs < 1 or args.chunk_reads < 1:
        parser.error("--jobs and --chunk-reads must be positive")

    config = load_config(args.config)
    plate = (
        prepare_barcode_panel(_trimmed_barcodes(config.plate_barcodes))
        if config.plate_barcodes.sequences
        else None
    )
    well = (
        prepare_barcode_panel(_trimmed_barcodes(config.well_barcodes))
        if config.well_barcodes.sequences
        else None
    )

    def jobs():
        for index, record in enumerate(iter_fastq(args.fastq)):
            if args.limit is not None and index >= args.limit:
                break
            yield record, "BENCHMARK", config, plate, well

    started = perf_counter()
    counts: Counter[str] = Counter()
    total = 0
    worker_count = 1 if args.backend == "serial" else args.jobs
    executor = "thread" if args.backend == "serial" else args.backend
    for batch in _ordered_map(
        _demux_batch,
        _batched(jobs(), args.chunk_reads),
        worker_count,
        backend=executor,
    ):
        total += len(batch)
        counts.update(row["call_status"] for row, _ in batch)
    elapsed = perf_counter() - started
    print(
        json.dumps(
            {
                "backend": args.backend,
                "jobs": worker_count,
                "chunk_reads": args.chunk_reads,
                "records": total,
                "seconds": round(elapsed, 6),
                "reads_per_second": round(total / elapsed, 3),
                "call_counts": dict(sorted(counts.items())),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    freeze_support()
    raise SystemExit(main())
