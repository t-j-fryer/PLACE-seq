#!/usr/bin/env python3
"""Ask whether a contested position reappears in an independent PCR of the same well.

Two culture plates were sequenced twice: once on a dedicated barcode and once
inside the pooled barcode.  Those are separate colony PCRs and separate library
preparations from the *same* well, so a variant seen in both cannot have been made
by either PCR or by the basecaller - it was in the well.  A variant seen in only
one was made after the well.

For every clone flagged as two-allele in one member of a replicate pair, this
reports the allele counts at the same reference position in the other member -
including when the other member is not flagged there, which is the case that
matters (a variant at 96% is not "contested", it is called).

    python scripts/replicate_alleles.py --run-dir runs/260608-full-length-v8c
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dissect_mixed_positions import load_run_config, pileup  # noqa: E402

from nanopore3 import pipeline  # noqa: E402
from nanopore3.references import read_reference_libraries  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--min-quality", type=int, default=10)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    config = load_run_config(args.run_dir, args.config)
    flanks = pipeline.resolve_flanks(config)
    config = pipeline.apply_flanks(config, flanks)
    collection = read_reference_libraries(
        {k: s.fasta for k, s in config.reference_sets.items()},
        transforms=pipeline.flank_transforms(flanks),
    )
    references = {
        library: {r.id: r.sequence.upper() for r in collection.get(library).records}
        for library in config.reference_sets
    }

    stages = args.run_dir / "stages"
    pairs = list(csv.DictReader(
        (args.run_dir / "figures" / "replicate_concordance_wells.csv").open(encoding="utf-8")
    ))
    # (culture plate, well) -> the two barcodes that carry it
    partner: dict[tuple[str, str], tuple[str, str]] = {
        (row["culture_plate"], row["well_id"]): (row["dedicated_plate"], row["pooled_plate"])
        for row in pairs
    }

    contested = json.load((args.run_dir / "figures" / "mixed_positions.json").open())

    # Which (plate, well, design) groups do we need reads for? The flagged clone
    # and its partner in the other barcode.
    consensus_rows: dict[str, dict[str, str]] = {}
    with gzip.open(stages / "04_consensus" / "consensus.csv.gz", "rt") as handle:
        for row in csv.DictReader(handle):
            consensus_rows[row["consensus_id"]] = row

    def normalise(well: str) -> str:
        letter, digits = well[0], well[1:]
        return f"{letter}{int(digits):02d}" if digits.isdigit() else well

    targets = []
    for item in contested:
        row = consensus_rows[item["consensus_id"]]
        culture, well = row.get("culture_plate", ""), normalise(row["well_id"])
        if "|" in culture or (culture, well) not in partner:
            continue
        dedicated, pooled = partner[(culture, well)]
        other = dedicated if row["plate_id"] == pooled else pooled
        targets.append((item, other, culture, well))
    print(f"{len(targets)} contested positions sit in a well that was sequenced twice")
    if not targets:
        return 0

    # Reads for the partner group, straight from the assignment stage: the partner
    # clone may not exist as a consensus at all, so contributors are no use here.
    need = {
        (other, well, item["design"]) for item, other, _culture, well in targets
    }
    reads: dict[tuple[str, str, str], list[tuple[str, str, list[int] | None]]] = defaultdict(list)
    with gzip.open(stages / "03_assignment" / "consensus_eligible.jsonl.gz", "rt") as handle:
        for line in handle:
            record = json.loads(line)
            key = (
                record["plate_id"],
                normalise(record["well_id"] or ""),
                # reference_ids is a list here and a "|"-joined string in the
                # consensus table; normalise to the first alias either way.
                (record["reference_ids"] or [""])[0]
                if isinstance(record["reference_ids"], list)
                else (record["reference_ids"] or "").split("|")[0],
            )
            if key not in need:
                continue
            quality = record.get("quality")
            reads[key].append(
                (
                    record["read_uid"],
                    record["sequence"],
                    [ord(s) - 33 for s in quality] if quality else None,
                )
            )
    print(f"recovered reads for {len(reads)} of {len(need)} partner groups")

    findings = []
    cache: dict[tuple[str, str, str], tuple] = {}
    for item, other, _culture, well in targets:
        key = (other, well, item["design"])
        group = reads.get(key)
        if not group:
            findings.append({**item, "partner_plate": other, "partner_reads": 0,
                             "partner_status": "no reads"})
            continue
        if key not in cache:
            reference = references[item["library"]][item["design"]]
            cache[key] = pileup(reference, group, args.min_quality)
        votes, _spanning, _deletions = cache[key]
        counter = votes[item["position"]]
        depth = sum(counter.values())
        same = counter.get(item["minor"], 0)
        findings.append(
            {
                **item,
                "partner_plate": other,
                "partner_reads": len(group),
                "partner_depth": depth,
                "partner_minor_allele_reads": same,
                "partner_minor_fraction": round(same / depth, 4) if depth else None,
                "partner_counts": dict(counter.most_common()),
                "partner_status": "ok",
            }
        )

    out = args.out or args.run_dir / "figures" / "replicate_alleles.json"
    out.write_text(json.dumps(findings, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
