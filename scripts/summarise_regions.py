#!/usr/bin/env python3
"""Per-region accuracy for a full-length run, and what an insert-only view hides.

A full-length consensus is scored over the whole amplicon, but the amplicon is
not one thing: the insert is the designed, synthesised part and the flanks are
the vector every clone shares.  This reports them separately, and counts the
clones that look perfect at the insert alone yet carry differences elsewhere.

    python scripts/summarise_regions.py --run-dir runs/260608-RP05-RP08-v5
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path

REGIONS = ("flank_5p", "insert", "flank_3p")


def load_rows(run_dir: Path) -> list[dict[str, str]]:
    with gzip.open(run_dir / "stages" / "05_qc" / "qc.csv.gz", "rt", encoding="utf-8") as handle:
        return [row for row in csv.DictReader(handle) if row.get("insert_identity")]


def library_of(run_dir: Path) -> dict[str, str]:
    """Consensus id to plate, so the two libraries can be reported apart."""

    path = run_dir / "stages" / "04_consensus" / "consensus.csv.gz"
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return {row["consensus_id"]: row["plate_id"] for row in csv.DictReader(handle)}


def summarise(rows: list[dict[str, str]]) -> dict[str, object]:
    out: dict[str, object] = {"consensuses": len(rows)}
    for region in REGIONS:
        identities = [float(row[f"{region}_identity"]) for row in rows]
        edits = [int(row[f"{region}_edit_distance"]) for row in rows]
        lengths = {int(row[f"{region}_length"]) for row in rows}
        out[region] = {
            "reference_length_median": int(statistics.median(lengths)),
            "mean_identity": round(100 * statistics.fmean(identities), 4),
            "median_identity": round(100 * statistics.median(identities), 4),
            "exact_percent": round(100 * sum(1 for e in edits if e == 0) / len(edits), 2),
            "errors_per_kb": round(
                1000 * sum(edits) / sum(int(r[f"{region}_length"]) for r in rows), 3
            ),
        }
    whole = [int(row["alignment_edit_distance"]) for row in rows]
    out["whole_amplicon"] = {
        "mean_identity": round(
            100 * statistics.fmean(float(row["alignment_identity"]) for row in rows), 4
        ),
        "exact_percent": round(100 * sum(1 for e in whole if e == 0) / len(whole), 2),
    }
    insert_perfect = [row for row in rows if int(row["insert_edit_distance"]) == 0]
    hidden = [
        row
        for row in insert_perfect
        if int(row["flank_5p_edit_distance"]) + int(row["flank_3p_edit_distance"]) > 0
    ]
    out["insert_perfect"] = len(insert_perfect)
    out["insert_perfect_with_flank_differences"] = len(hidden)
    out["insert_perfect_with_flank_differences_percent"] = (
        round(100 * len(hidden) / len(insert_perfect), 2) if insert_perfect else 0.0
    )
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, help="default <run-dir>/figures/region_summary.json")
    args = parser.parse_args()

    rows = load_rows(args.run_dir)
    if not rows:
        print("no per-region columns in this run; it was not a full-length run")
        return 1
    plates = library_of(args.run_dir)
    by_plate: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_plate[plates.get(row["consensus_id"], "?")].append(row)

    report = {"all": summarise(rows)}
    for plate, subset in sorted(by_plate.items()):
        report[plate] = summarise(subset)

    for scope, values in report.items():
        print(f"\n=== {scope} ({values['consensuses']:,} consensuses) ===")
        print(f"  {'region':12s} {'len':>6s} {'mean id':>9s} {'exact':>8s} {'err/kb':>8s}")
        for region in (*REGIONS, "whole_amplicon"):
            data = values[region]
            length = data.get("reference_length_median", "")
            print(
                f"  {region:12s} {str(length):>6s} {data['mean_identity']:>8.3f}% "
                f"{data['exact_percent']:>7.1f}% {str(data.get('errors_per_kb','')):>8s}"
            )
        print(
            f"  insert-perfect clones: {values['insert_perfect']:,}; of those, "
            f"{values['insert_perfect_with_flank_differences']:,} "
            f"({values['insert_perfect_with_flank_differences_percent']}%) differ from the "
            "reference somewhere in the constant region"
        )

    out = args.out or args.run_dir / "figures" / "region_summary.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
