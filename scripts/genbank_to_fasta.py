#!/usr/bin/env python3
"""Convert GenBank/ApE records to FASTA, optionally extracting a motif region.

Plasmid maps often arrive as ``.gb`` files, which the pipeline cannot read.  With
``--forward-motif``/``--reverse-motif`` this writes the region between the two
motifs instead of the whole record, which is what the reference-guided path
should be compared against: scoring identity across a whole 6 kb vector measures
the invariant backbone rather than the variable insert.

Records are treated as circular, so a region spanning the origin is still found.
"""

from __future__ import annotations

import argparse
import re
import textwrap
from pathlib import Path

_ORIGIN = re.compile(r"^ORIGIN.*?$(.*?)^//", re.M | re.S)
_NOT_BASE = re.compile(r"[^acgtACGTnN]")
_UNSAFE = re.compile(r"[^A-Za-z0-9._-]+")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def read_genbank(path: Path) -> str | None:
    """Return the ORIGIN sequence of a GenBank record, or None if absent."""

    match = _ORIGIN.search(path.read_text(errors="replace"))
    if match is None:
        return None
    return "".join(_NOT_BASE.sub("", line) for line in match.group(1).splitlines()).upper()


def extract_region(sequence: str, forward: str, reverse: str, limit: int) -> str | None:
    """Sequence strictly between two motifs, searching both strands circularly."""

    for strand in (sequence, reverse_complement(sequence)):
        doubled = strand + strand
        start = doubled.find(forward)
        if start < 0:
            continue
        start += len(forward)
        end = doubled.find(reverse, start)
        if 0 <= end and end - start <= limit:
            return doubled[start:end]
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path, required=True,
        help="a .gb file or a directory of them",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--forward-motif")
    parser.add_argument("--reverse-motif")
    parser.add_argument("--max-region", type=int, default=20000)
    args = parser.parse_args()
    if (args.forward_motif is None) != (args.reverse_motif is None):
        parser.error("--forward-motif and --reverse-motif must be given together")

    sources = (
        sorted(args.input.glob("*.gb")) if args.input.is_dir() else [args.input]
    )
    if not sources:
        parser.error(f"no .gb files found in {args.input}")

    records: dict[str, str] = {}
    skipped: list[str] = []
    for source in sources:
        sequence = read_genbank(source)
        name = _UNSAFE.sub("_", source.stem).strip("_")
        if not sequence:
            skipped.append(f"{name} (no ORIGIN block)")
            continue
        if args.forward_motif:
            region = extract_region(
                sequence, args.forward_motif.upper(), args.reverse_motif.upper(),
                args.max_region,
            )
            if region is None:
                skipped.append(f"{name} (motif pair not found)")
                continue
            sequence = region
        records[name] = sequence

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="ascii") as handle:
        for name, sequence in sorted(records.items()):
            handle.write(f">{name}\n" + "\n".join(textwrap.wrap(sequence, 60)) + "\n")
    lengths = sorted(len(s) for s in records.values())
    print(f"wrote {len(records)} records to {args.output}")
    if lengths:
        print(f"  length min={lengths[0]} median={lengths[len(lengths)//2]} max={lengths[-1]}")
        print(f"  distinct sequences: {len(set(records.values()))}/{len(records)}")
    if skipped:
        print(f"  skipped {len(skipped)}:")
        for item in skipped:
            print(f"    {item}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
