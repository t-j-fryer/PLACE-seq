#!/usr/bin/env python3
"""Separate the oversampling effect from the benefit of a second DNA encoding.

Two encodings of one design set recover more designs than either alone, but they
also spend twice the colony-picking effort, so the headline comparison confounds
the two. This separates them three ways:

1. **As sampled** - what each set actually recovered, at its own effort.
2. **At matched effort** - the union interpolated down to the effort one encoding
   spent. Interpolation needs no model: it re-samples the wells that exist.
3. **The ceiling** - the designs each encoding produced *at all*, at any depth.
   A design that never appeared in a single read cannot be found by picking more
   colonies, so this bounds what any amount of sampling could reach.

Extrapolating one encoding upwards is deliberately not the headline. Chao1 and
its relatives assume unseen designs are rare rather than absent, and here that is
false: 19 designs never appear in encoding A at all, so Chao1's 98.3% asymptote
for A is 4 points above what A can actually make.

    python scripts/rarefy_encodings.py --run-dir runs/260608-full-length-v7
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
from collections import defaultdict
from math import exp, lgamma
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from nanopore3.figures import MM, save_figure, use_journal_style  # noqa: E402

PREFIX = re.compile(r"^(?:([AB])_)?Block_(\d+)_")
COLOURS = {"A": "#0072B2", "B": "#D55E00", "A+B": "#009E73"}
NARROW = 89 * MM


def log_choose(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)


def rarefy(well_counts: list[int], total_wells: int, n: int) -> float:
    """Expected designs recovered when n of the wells are sampled, without replacement."""

    if n >= total_wells:
        return float(len(well_counts))
    denominator = log_choose(total_wells, n)
    expected = 0.0
    for count in well_counts:
        if total_wells - count < n:
            expected += 1.0
        else:
            expected += 1.0 - exp(log_choose(total_wells - count, n) - denominator)
    return expected


def collect(run_dir: Path, plate: str, universe: int):
    """Per encoding: wells per recovered design, and every design seen at all."""

    stages = run_dir / "stages"
    with gzip.open(stages / "05_qc" / "qc.csv.gz", "rt", encoding="utf-8") as handle:
        exact = {
            row["consensus_id"]
            for row in csv.DictReader(handle)
            if row.get("insert_edit_distance") == "0"
        }
    wells: dict[str, dict[str, set]] = defaultdict(lambda: defaultdict(set))
    with gzip.open(stages / "04_consensus" / "consensus.csv.gz", "rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["plate_id"] != plate or row["status"] not in (
                "consensus_pass",
                "mixed_variants",
            ):
                continue
            match = PREFIX.match(row["reference_ids"])
            if not match or match.group(1) not in ("A", "B"):
                continue
            encoding = match.group(1)
            if not row["culture_plate"].startswith(f"SUMO_{encoding}_"):
                continue
            if row["consensus_id"] in exact:
                key = row["reference_ids"][match.end() :]
                wells[encoding][key].add((row["culture_plate"], row["well_id"]))

    produced: dict[str, set] = defaultdict(set)
    with gzip.open(
        stages / "03_assignment" / "assignment_calls.csv.gz", "rt", encoding="utf-8"
    ) as handle:
        for row in csv.DictReader(handle):
            if row["plate_id"] != plate:
                continue
            match = PREFIX.match(row["reference_ids"] or "")
            if not match or match.group(1) not in ("A", "B"):
                continue
            encoding = match.group(1)
            if not row["culture_plate"].startswith(f"SUMO_{encoding}_"):
                continue
            produced[encoding].add(row["reference_ids"][match.end() :])
    return wells, produced


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--plate", default="RP05")
    parser.add_argument("--universe", type=int, default=342, help="designs in the library")
    parser.add_argument("--wells-per-encoding", type=int, default=11 * 95)
    parser.add_argument("--out-dir", type=Path)
    args = parser.parse_args()
    out_dir = args.out_dir or args.run_dir / "figures"

    wells, produced = collect(args.run_dir, args.plate, args.universe)
    per_encoding = args.wells_per_encoding
    sets = {
        "A": ({d: len(w) for d, w in wells["A"].items()}, per_encoding, produced["A"]),
        "B": ({d: len(w) for d, w in wells["B"].items()}, per_encoding, produced["B"]),
    }
    union: dict[str, int] = defaultdict(int)
    for encoding in ("A", "B"):
        for design, ws in wells[encoding].items():
            union[design] += len(ws)
    sets["A+B"] = (dict(union), 2 * per_encoding, produced["A"] | produced["B"])

    report: dict[str, object] = {"universe": args.universe, "sets": {}}
    print(f"library: {args.universe} designs;  wells picked per encoding: {per_encoding}\n")
    print(f"{'set':5s} {'effort':>7s} {'recovered':>10s} {'as sampled':>12s} "
          f"{'at 1x effort':>13s} {'ceiling':>9s}")
    for name, (counts, effort, ever) in sets.items():
        values = list(counts.values())
        matched = rarefy(values, effort, per_encoding)
        report["sets"][name] = {
            "wells": effort,
            "recovered": len(counts),
            "recovered_percent": round(100 * len(counts) / args.universe, 2),
            "at_matched_effort": round(matched, 1),
            "at_matched_effort_percent": round(100 * matched / args.universe, 2),
            "ceiling": len(ever),
            "ceiling_percent": round(100 * len(ever) / args.universe, 2),
        }
        print(f"{name:5s} {effort:7d} {len(counts):10d} "
              f"{100 * len(counts) / args.universe:11.1f}% "
              f"{100 * matched / args.universe:12.1f}% "
              f"{100 * len(ever) / args.universe:8.1f}%")

    a, b = sets["A"][2], sets["B"][2]
    report["structural"] = {
        "produced_only_by_A": len(a - b),
        "produced_only_by_B": len(b - a),
        "produced_by_neither": args.universe - len(a | b),
    }
    print(f"\ndesigns produced by only one encoding: A {len(a - b)}, B {len(b - a)}; "
          f"by neither {args.universe - len(a | b)}")
    headline = 100 * (sets["A+B"][0] and len(sets["A+B"][0]) - len(sets["A"][0])) / args.universe
    ceiling_gain = 100 * (len(a | b) - len(a)) / args.universe
    print(f"\nheadline gain of A+B over A as sampled : {headline:.1f} points")
    print(f"  of which the ceiling can explain      : {ceiling_gain:.1f} points (structural)")
    print(f"  the remainder is oversampling         : {headline - ceiling_gain:.1f} points")
    report["decomposition"] = {
        "headline_gain_points": round(headline, 2),
        "structural_points": round(ceiling_gain, 2),
        "oversampling_points": round(headline - ceiling_gain, 2),
    }

    use_journal_style()
    figure, axes = plt.subplots(figsize=(NARROW, 2.1))
    grid = list(range(1, 2 * per_encoding + 1, 10))
    for name, (counts, effort, ever) in sets.items():
        values = list(counts.values())
        xs = [n for n in grid if n <= effort]
        ys = [100 * rarefy(values, effort, n) / args.universe for n in xs]
        axes.plot(
            xs, ys, color=COLOURS[name], linewidth=1.1, zorder=3,
            label=f"Library 3 ({name})",
        )
        axes.axhline(100 * len(ever) / args.universe, color=COLOURS[name],
                     linewidth=0.6, linestyle=(0, (3, 2)), zorder=2)
    axes.axvline(per_encoding, color="black", linewidth=0.6, linestyle=(0, (1, 2)), zorder=1)
    axes.annotate(
        "effort of one\nencoding", xy=(per_encoding, 12),
        xytext=(per_encoding * 0.62, 12),
        fontsize=5.6, ha="right", va="center", color="black",
    )
    axes.annotate("dashed: designs each route ever produced", xy=(2 * per_encoding, 99),
                  xytext=(2 * per_encoding, 60), fontsize=5.6, ha="right", va="center")
    axes.set_xlabel("Wells picked")
    axes.set_ylabel("Designs recovered (%)")
    axes.set_ylim(0, 105)
    axes.set_yticks([0, 25, 50, 75, 100])
    axes.set_xlim(0, 2 * per_encoding)
    axes.tick_params(axis="both", direction="in", length=3, width=0.6)
    for side in ("top", "right"):
        axes.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        axes.spines[side].set_linewidth(0.6)
    axes.legend(loc="lower right", frameon=False, fontsize=6, handlelength=1.4,
                borderpad=0.2, labelspacing=0.3)
    figure.subplots_adjust(left=0.145, right=0.99, top=0.97, bottom=0.20)
    path = save_figure(figure, out_dir / "encoding_rarefaction")

    (out_dir / "encoding_rarefaction.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"\nwrote {path}, {path.with_suffix('.png')}, {path.with_suffix('.svg')},")
    print(f"      {out_dir / 'encoding_rarefaction.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
