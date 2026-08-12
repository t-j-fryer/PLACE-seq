#!/usr/bin/env python3
"""Create a fast, deterministic diagnostic sample from a large plain FASTQ.

Offsets are uniform over bytes, so long records are overrepresented.  This is
appropriate for a bounded workflow smoke test, not for estimating biological
read frequencies.  The source file is opened read-only and never modified.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import random


DNA = frozenset(b"ACGTRYSWKMBDHVNacgtryswkmbdhvn")


def _candidate(handle, offset: int) -> tuple[int, bytes] | None:
    handle.seek(offset)
    if offset:
        handle.readline()
    for _ in range(32):
        position = handle.tell()
        header = handle.readline()
        if not header:
            return None
        if not header.startswith(b"@"):
            continue
        sequence = handle.readline()
        plus = handle.readline()
        quality = handle.readline()
        seq = sequence.rstrip(b"\r\n")
        qual = quality.rstrip(b"\r\n")
        if (
            seq
            and plus.startswith(b"+")
            and len(seq) == len(qual)
            and not (set(seq) - DNA)
        ):
            return position, header + sequence + plus + quality
    return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--records", type=int, default=20_000)
    parser.add_argument("--seed", type=int, default=141142)
    parser.add_argument("--source-sha256")
    args = parser.parse_args()
    source = args.source.expanduser().resolve(strict=True)
    output = args.output.expanduser().resolve(strict=False)
    if output.exists() or output.with_suffix(output.suffix + ".sample.json").exists():
        parser.error("output or sample manifest already exists")
    if args.records < 1:
        parser.error("--records must be positive")
    size = source.stat().st_size
    rng = random.Random(args.seed)
    selected: dict[int, bytes] = {}
    attempts = 0
    limit = max(args.records * 20, 1_000)
    with source.open("rb") as handle:
        while len(selected) < args.records and attempts < limit:
            attempts += 1
            record = _candidate(handle, rng.randrange(size))
            if record is not None:
                selected.setdefault(*record)
    if len(selected) != args.records:
        raise RuntimeError(
            f"found only {len(selected):,} unique records after {attempts:,} attempts"
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256()
    with output.open("xb") as destination:
        for position in sorted(selected):
            record = selected[position]
            destination.write(record)
            digest.update(record)
    manifest = {
        "schema_version": 1,
        "method": "deterministic_uniform_byte_offset_resynchronization",
        "statistical_caveat": "length-biased; diagnostic use only",
        "source": str(source),
        "source_size_bytes": size,
        "source_sha256": args.source_sha256,
        "output": str(output),
        "output_sha256": digest.hexdigest(),
        "records": len(selected),
        "seed": args.seed,
        "attempts": attempts,
        "source_byte_offsets_sha256": hashlib.sha256(
            "\n".join(str(value) for value in sorted(selected)).encode("ascii")
        ).hexdigest(),
    }
    manifest_path = output.with_suffix(output.suffix + ".sample.json")
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
