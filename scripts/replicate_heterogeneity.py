#!/usr/bin/env python3
"""Separate variation that was in the source well from variation made after it.

Two culture plates were sequenced twice, and the two routes were grown as separate
cultures.  That makes them an experiment on where a minority allele comes from:

* a minority allele **present in both** outgrowths was in whatever was picked -
  a heterogeneous colony, or a mutation that arose before the two were split.
  Independent cultures cannot invent the same base at the same position.
* a minority allele **present in one** was made afterwards - during that culture,
  that colony PCR, that library prep, or that basecall.

Run over every position of every clone the two routes share, this measures how
much of the observed heterogeneity is biological, and - by splitting on the
substitution - tests the prediction that the G:C -> T:A class behaves differently
from the rest.  8-oxoguanine damage happens to a DNA sample, so it should be
unshared; a real mutation in the picked cell should be shared.

    python scripts/replicate_heterogeneity.py --run-dir runs/260608-full-length-v8c
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
sys.path.insert(0, str(Path(__file__).resolve().parent))

from dissect_mixed_positions import load_run_config, pileup  # noqa: E402

from nanopore3 import pipeline  # noqa: E402
from nanopore3.references import read_reference_libraries  # noqa: E402

TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
OXO = {("G", "T"), ("C", "A")}


def normalise(well: str) -> str:
    letter, digits = well[:1], well[1:]
    return f"{letter}{int(digits):02d}" if digits.isdigit() else well


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--min-depth", type=int, default=20)
    parser.add_argument("--present", type=float, default=0.05,
                        help="minority fraction counted as present")
    parser.add_argument("--absent", type=float, default=0.02,
                        help="minority fraction counted as absent; between the two is unscored")
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
        lib: {r.id: r.sequence.upper() for r in collection.get(lib).records}
        for lib in config.reference_sets
    }

    pairs = list(csv.DictReader(
        (args.run_dir / "figures" / "replicate_concordance_wells.csv").open(encoding="utf-8")
    ))
    # Only the clones both routes actually recovered can be compared.
    groups: list[tuple[str, str, str, str, str]] = []
    for row in pairs:
        for design in filter(None, row["shared_clones"].split("|")):
            groups.append(
                (row["culture_plate"], normalise(row["well_id"]),
                 design, row["dedicated_plate"], row["pooled_plate"])
            )
    print(f"{len(groups)} clones recovered by both routes, across {len(pairs)} wells")

    need = set()
    for _culture, well, design, dedicated, pooled in groups:
        need.add((dedicated, well, design))
        need.add((pooled, well, design))

    reads: dict[tuple[str, str, str], list] = defaultdict(list)
    with gzip.open(
        args.run_dir / "stages" / "03_assignment" / "consensus_eligible.jsonl.gz", "rt"
    ) as handle:
        for line in handle:
            record = json.loads(line)
            ids = record["reference_ids"]
            design = (ids or [""])[0] if isinstance(ids, list) else (ids or "").split("|")[0]
            key = (record["plate_id"], normalise(record["well_id"] or ""), design)
            if key not in need:
                continue
            quality = record.get("quality")
            reads[key].append(
                (record["read_uid"], record["sequence"],
                 [ord(s) - 33 for s in quality] if quality else None)
            )
    print(f"recovered reads for {len(reads)} of {len(need)} groups")

    findings = []
    for culture, well, design, dedicated, pooled in groups:
        library = next(
            (lib for lib in references if design in references[lib]), None
        )
        if library is None:
            continue
        reference = references[library][design]
        left, right = reads.get((dedicated, well, design)), reads.get((pooled, well, design))
        if not left or not right:
            continue
        votes_l, _s, _d = pileup(reference, left, args.min_quality)
        votes_r, _s2, _d2 = pileup(reference, right, args.min_quality)
        for index, base in enumerate(reference):
            counter_l, counter_r = votes_l[index], votes_r[index]
            depth_l, depth_r = sum(counter_l.values()), sum(counter_r.values())
            if depth_l < args.min_depth or depth_r < args.min_depth:
                continue
            for allele in "ACGT":
                if allele == base:
                    continue
                fl_, fr_ = counter_l.get(allele, 0) / depth_l, counter_r.get(allele, 0) / depth_r
                if max(fl_, fr_) < args.present:
                    continue
                if min(fl_, fr_) >= args.present:
                    verdict = "shared"
                elif min(fl_, fr_) < args.absent:
                    verdict = "one_route_only"
                else:
                    verdict = "ambiguous"
                findings.append({
                    "culture_plate": culture, "well_id": well, "design": design,
                    "position": index, "reference_base": base, "allele": allele,
                    "dedicated_fraction": round(fl_, 4), "pooled_fraction": round(fr_, 4),
                    "dedicated_depth": depth_l, "pooled_depth": depth_r,
                    "verdict": verdict,
                    "class": ("oxo_GC_to_TA" if (base, allele) in OXO
                              else "transition" if (base, allele) in TRANSITIONS
                              else "other_transversion"),
                })

    out = args.out or args.run_dir / "figures" / "replicate_heterogeneity.json"
    out.write_text(json.dumps(findings, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"{len(findings)} minority alleles at >= {args.present:.0%} in at least one route")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
