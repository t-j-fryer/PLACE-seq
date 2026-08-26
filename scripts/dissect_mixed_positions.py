#!/usr/bin/env python3
"""Characterise every contested position in the clones flagged as two-allele.

A clone graded `mixed_variants` says the reads for one design disagree at a
position beyond the minor-allele floor.  The interesting question is *why*, and
the candidate explanations make different predictions about where those positions
sit and what the minor fraction is:

* an early-cycle colony-PCR error  -> minor fractions clustered near 1/2, 1/4, 1/8,
  positions scattered, no reuse of the same coordinate between wells
* a systematic basecalling error   -> the same coordinate contested in many wells,
  and a distinctive sequence context (homopolymer, Dam `GATC`, Dcm `CCWGG`)
* a real mixed plasmid population  -> arbitrary fractions, scattered positions, and
  the same allele seen again when the same well is sequenced twice

So this reports, per contested position: where it is (insert or which flank, and
how far into it), the two alleles and their read counts, and the local sequence
context.  Aggregation is left to the caller - the JSON is the point.

    python scripts/dissect_mixed_positions.py --run-dir runs/260608-full-length-v8c
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import edlib  # noqa: E402

from nanopore3 import pipeline  # noqa: E402
from nanopore3.config import load_config  # noqa: E402
from nanopore3.references import read_reference_libraries  # noqa: E402

_CIGAR_TOKEN = re.compile(r"(\d+)([=XID])")
_DCM = re.compile(r"CC[AT]GG")


def load_run_config(run_dir: Path, override: Path | None):
    if override:
        return load_config(override)
    name = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))["config"]["run_name"]
    for candidate in sorted(Path("configs/runs").glob("*.yaml")):
        if load_config(candidate).run_name == name:
            return load_config(candidate)
    raise SystemExit("could not find the run's config; pass --config")


def homopolymer_run(sequence: str, index: int) -> int:
    """Length of the homopolymer run the position sits in."""

    base = sequence[index]
    left = index
    while left > 0 and sequence[left - 1] == base:
        left -= 1
    right = index
    while right + 1 < len(sequence) and sequence[right + 1] == base:
        right += 1
    return right - left + 1


def pileup(reference: str, reads: list[tuple[str, str, list[int] | None]], min_quality: int):
    """Base counts, spanning depth and deletion counts per reference position.

    The same CIGAR walk the portable consensus backend uses: quality gates
    whether a base is *counted*, never how many reads were *there*.
    """

    votes: list[Counter[str]] = [Counter() for _ in reference]
    spanning = [0] * len(reference)
    deletions = [0] * len(reference)
    for _uid, sequence, qualities in reads:
        query = sequence.upper()
        result = edlib.align(query, reference, mode="NW", task="path")
        cigar = result.get("cigar")
        if not cigar:
            continue
        qi = ri = 0
        for count_text, operation in _CIGAR_TOKEN.findall(cigar):
            count = int(count_text)
            if operation in "=X":
                for _ in range(count):
                    spanning[ri] += 1
                    if qualities is None or qualities[qi] >= min_quality:
                        votes[ri][query[qi]] += 1
                    qi += 1
                    ri += 1
            elif operation == "I":
                qi += count
            else:
                for _ in range(count):
                    spanning[ri] += 1
                    deletions[ri] += 1
                    ri += 1
    return votes, spanning, deletions


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--minor-fraction", type=float, default=0.15)
    parser.add_argument("--min-quality", type=int, default=10)
    parser.add_argument("--context", type=int, default=10)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    config = load_run_config(args.run_dir, args.config)
    flanks = pipeline.resolve_flanks(config)
    config = pipeline.apply_flanks(config, flanks)
    collection = read_reference_libraries(
        {k: s.fasta for k, s in config.reference_sets.items()},
        transforms=pipeline.flank_transforms(flanks),
    )
    references: dict[str, dict[str, str]] = {}
    spans: dict[str, dict[str, tuple[int, int]]] = {}
    for library_id in config.reference_sets:
        records = collection.get(library_id).records
        references[library_id] = {r.id: r.sequence.upper() for r in records}
        # The insert span is a property of the flank model and the reference
        # length, not of the record: ask the pipeline for it, so "insert" here
        # means exactly what it means in the QC table.
        flank = flanks.get(library_id) if isinstance(flanks, dict) else flanks
        spans[library_id] = {
            r.id: pipeline.qc_regions(flank, len(r.sequence))["insert"] for r in records
        }

    stages = args.run_dir / "stages"
    flagged: dict[str, dict[str, str]] = {}
    with gzip.open(stages / "04_consensus" / "consensus.csv.gz", "rt") as handle:
        for row in csv.DictReader(handle):
            if (row.get("mixed_positions") or "0") not in ("0", ""):
                flagged[row["consensus_id"]] = row
    print(f"{len(flagged)} clones flagged with a contested position")

    wanted: dict[str, str] = {}
    with gzip.open(stages / "04_consensus" / "contributors.csv.gz", "rt") as handle:
        for row in csv.DictReader(handle):
            if row["consensus_id"] in flagged:
                wanted[row["read_uid"]] = row["consensus_id"]
    print(f"{len(wanted)} contributing reads to recover")

    by_clone: dict[str, list[tuple[str, str, list[int] | None]]] = defaultdict(list)
    with gzip.open(stages / "03_assignment" / "consensus_eligible.jsonl.gz", "rt") as handle:
        for line in handle:
            uid_at = line.find('"read_uid"')
            if uid_at < 0:
                continue
            record = json.loads(line)
            clone = wanted.get(record["read_uid"])
            if clone is None:
                continue
            quality = record.get("quality")
            by_clone[clone].append(
                (
                    record["read_uid"],
                    record["sequence"],
                    # Phred+33, the same decoding the consensus stage uses.
                    [ord(symbol) - 33 for symbol in quality] if quality else None,
                )
            )
    print(f"recovered reads for {len(by_clone)} clones")

    findings = []
    for consensus_id, row in flagged.items():
        reads = by_clone.get(consensus_id)
        if not reads:
            continue
        library = row["reference_library_id"]
        design = row["reference_ids"].split("|")[0]
        reference = references.get(library, {}).get(design)
        if not reference:
            continue
        start, end = spans[library][design]
        votes, spanning, deletions = pileup(reference, reads, args.min_quality)
        for index, counter in enumerate(votes):
            depth = sum(counter.values())
            if depth < 10:
                continue
            ranked = counter.most_common()
            if len(ranked) < 2:
                continue
            (major, major_n), (minor, minor_n) = ranked[0], ranked[1]
            fraction = minor_n / depth
            if fraction < args.minor_fraction:
                continue
            lo = max(0, index - args.context)
            hi = min(len(reference), index + args.context + 1)
            if index < start:
                region, offset = "flank_5", index
            elif index >= end:
                region, offset = "flank_3", index - end
            else:
                region, offset = "insert", index - start
            findings.append(
                {
                    "consensus_id": consensus_id,
                    "plate_id": row["plate_id"],
                    "well_id": row["well_id"],
                    "culture_plate": row.get("culture_plate", ""),
                    "design": design,
                    "library": library,
                    "position": index,
                    "region": region,
                    "region_offset": offset,
                    "reference_base": reference[index],
                    "major": major,
                    "major_reads": major_n,
                    "minor": minor,
                    "minor_reads": minor_n,
                    "depth": depth,
                    "minor_fraction": round(fraction, 4),
                    "deletion_reads": deletions[index],
                    "spanning": spanning[index],
                    "matches_reference": major == reference[index],
                    "context": reference[lo:hi],
                    "homopolymer": homopolymer_run(reference, index),
                    "dam": "GATC" in reference[max(0, index - 3) : index + 4],
                    "dcm": bool(_DCM.search(reference[max(0, index - 4) : index + 5])),
                }
            )

    print(f"{len(findings)} contested positions at minor fraction >= {args.minor_fraction}")
    out = args.out or args.run_dir / "figures" / "mixed_positions.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(findings, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
