#!/usr/bin/env python3
"""Compare two DNA encodings of the same protein designs.

Every design in this library was ordered twice, as two different DNA sequences
(A and B), assembled separately and picked separately. That makes the pair a
controlled experiment on the DNA rather than the protein: if a design fails in one
encoding and succeeds in the other, the protein cannot be the reason.

Reports:

* the outcome matrix - what happened to each design in A against what happened in
  B, over four states: recovered perfectly, recovered imperfectly, seen only in
  reads, never seen at all;
* whether failures are **independent** between encodings, which distinguishes a
  per-molecule accident from a property of the design;
* whether anything about the DNA predicts failure - GC content, longest
  homopolymer, longest G-run, and the number of synthesised fragments the design
  needed.

    python scripts/encoding_defects.py --run-dir runs/260608-full-length-v7 \
        --fragments "<...>/SUMO_AB_FULL_INFO.csv"
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import random
import re
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nanopore3 import pipeline  # noqa: E402
from nanopore3.config import load_config  # noqa: E402
from nanopore3.references import read_reference_libraries  # noqa: E402

PREFIX = re.compile(r"^([AB])_Block_(\d+)_")
STATES = ("perfect", "imperfect", "reads only", "absent")


def gc_percent(sequence: str) -> float:
    return 100 * (sequence.count("G") + sequence.count("C")) / len(sequence)


def longest_homopolymer(sequence: str) -> int:
    best = current = 1
    for previous, base in zip(sequence, sequence[1:], strict=False):
        current = current + 1 if base == previous else 1
        best = max(best, current)
    return best


def longest_run(sequence: str, base: str) -> int:
    return max((len(m) for m in re.findall(f"{base}+", sequence)), default=0)


def permutation_p(
    failures: list[float],
    recoveries: list[float],
    *,
    trials: int,
    seed: int,
) -> float:
    """One-sided p for failures scoring higher, by shuffling the labels."""

    observed = statistics.mean(failures) - statistics.mean(recoveries)
    pool = failures + recoveries
    rng = random.Random(seed)
    hits = 0
    for _ in range(trials):
        rng.shuffle(pool)
        shuffled = statistics.mean(pool[: len(failures)]) - statistics.mean(
            pool[len(failures) :]
        )
        if shuffled >= observed:
            hits += 1
    return hits / trials


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=Path("configs/runs/260608_full_length.yaml"))
    parser.add_argument("--library", default="sumo_ab")
    parser.add_argument("--plate", default="RP05")
    parser.add_argument("--fragments", type=Path, help="oPool FULL_INFO table, for fragment counts")
    parser.add_argument("--trials", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    config = load_config(args.config)
    flanks = pipeline.resolve_flanks(config)
    config = pipeline.apply_flanks(config, flanks)
    collection = read_reference_libraries(
        {k: s.fasta for k, s in config.reference_sets.items()},
        transforms=pipeline.flank_transforms(flanks),
    )
    full = {r.id: r.sequence for r in collection.get(args.library).records}
    inserts = pipeline.insert_view({args.library: full}, flanks)[0][args.library]
    by_key: dict[tuple[str, str], str] = {}
    for reference_id in inserts:
        match = PREFIX.match(reference_id)
        if match:
            by_key[(match.group(1), reference_id[match.end() :])] = reference_id
    designs = sorted({key for _, key in by_key})

    stages = args.run_dir / "stages"
    with gzip.open(stages / "05_qc" / "qc.csv.gz", "rt", encoding="utf-8") as handle:
        exact = {
            row["consensus_id"]
            for row in csv.DictReader(handle)
            if row.get("insert_edit_distance") == "0"
        }
    state: dict[str, dict[str, str]] = defaultdict(dict)
    with gzip.open(stages / "04_consensus" / "consensus.csv.gz", "rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["plate_id"] != args.plate or row["status"] not in (
                "consensus_pass",
                "mixed_variants",
            ):
                continue
            match = PREFIX.match(row["reference_ids"])
            if not match or not row["culture_plate"].startswith(f"SUMO_{match.group(1)}_"):
                continue
            encoding, key = match.group(1), row["reference_ids"][match.end() :]
            if row["consensus_id"] in exact or state[key].get(encoding) == "perfect":
                state[key][encoding] = "perfect"
            else:
                state[key].setdefault(encoding, "imperfect")
    seen: dict[str, set[str]] = defaultdict(set)
    with gzip.open(
        stages / "03_assignment" / "assignment_calls.csv.gz", "rt", encoding="utf-8"
    ) as handle:
        for row in csv.DictReader(handle):
            if row["plate_id"] != args.plate:
                continue
            match = PREFIX.match(row["reference_ids"] or "")
            if not match or not row["culture_plate"].startswith(f"SUMO_{match.group(1)}_"):
                continue
            seen[match.group(1)].add(row["reference_ids"][match.end() :])

    def outcome(key: str, encoding: str) -> str:
        if encoding in state[key]:
            return state[key][encoding]
        return "reads only" if key in seen[encoding] else "absent"

    matrix = Counter((outcome(k, "A"), outcome(k, "B")) for k in designs)
    total = len(designs)
    perfect_a = sum(1 for k in designs if outcome(k, "A") == "perfect")
    perfect_b = sum(1 for k in designs if outcome(k, "B") == "perfect")
    both = matrix[("perfect", "perfect")]
    either = sum(1 for k in designs if "perfect" in (outcome(k, "A"), outcome(k, "B")))
    expected = perfect_a * perfect_b / total

    print(f"{total} protein designs, each ordered as two independent DNA sequences\n")
    header = "A \\ B"
    print(f"{header:>12s}" + "".join(f"{s:>12s}" for s in STATES) + f"{'total':>8s}")
    for a in STATES:
        row = [matrix[(a, b)] for b in STATES]
        print(f"{a:>12s}" + "".join(f"{v:12d}" for v in row) + f"{sum(row):8d}")
    print(f"{'total':>12s}" + "".join(
        f"{sum(matrix[(a, b)] for a in STATES):12d}" for b in STATES))
    print(f"\nperfect in A {perfect_a} ({100*perfect_a/total:.1f}%), "
          f"B {perfect_b} ({100*perfect_b/total:.1f}%), "
          f"both {both} ({100*both/total:.1f}%), either {either} ({100*either/total:.1f}%)")
    print(f"expected in both if the encodings failed independently: {expected:.0f}")
    verdict = (
        "correlated: the design matters"
        if abs(both - expected) > 15
        else "independent: the molecule matters, not the design"
    )
    print(f"  observed {both}, a difference of {both - expected:+.0f} - {verdict}")

    report: dict[str, object] = {
        "designs": total,
        "matrix": {
            f"A={a}|B={b}": matrix[(a, b)]
            for a in STATES
            for b in STATES
            if matrix[(a, b)]
        },
        "perfect": {"A": perfect_a, "B": perfect_b, "both": both, "either": either,
                    "expected_both_if_independent": round(expected, 1)},
    }

    print("\ndoes anything about the DNA predict failure?")
    report["predictors"] = {}
    for encoding in ("A", "B"):
        good = [
            inserts[by_key[(encoding, k)]]
            for k in designs
            if outcome(k, encoding) == "perfect"
        ]
        bad = [inserts[by_key[(encoding, k)]] for k in designs if outcome(k, encoding) != "perfect"]
        entry = {"recovered": len(good), "failed": len(bad)}
        for label, fn in (("gc_percent", gc_percent),
                          ("longest_homopolymer", longest_homopolymer),
                          ("longest_g_run", lambda s: longest_run(s, "G"))):
            g = [fn(s) for s in good]
            b = [fn(s) for s in bad]
            p = permutation_p(b, g, trials=args.trials, seed=args.seed)
            entry[label] = {
                "recovered": round(statistics.mean(g), 2),
                "failed": round(statistics.mean(b), 2),
                "p": p,
            }
            print(
                f"   {encoding}  {label:20s} recovered {statistics.mean(g):6.2f}  "
                f"failed {statistics.mean(b):6.2f}   p = {p:.3f}"
            )
        report["predictors"][encoding] = entry

    if args.fragments and args.fragments.is_file():
        parts: dict[tuple[str, str], int] = {}
        with args.fragments.open("r", encoding="utf-8-sig") as handle:
            for row in csv.DictReader(handle):
                encoding = (row.get("Reference Source Label") or "").strip()
                count = sum(1 for i in (1, 2, 3) if (row.get(f"DNA Fragment {i}") or "").strip())
                parts[(encoding, row["Sequence Name"])] = count
        table: dict[int, list[int]] = defaultdict(lambda: [0, 0])
        for encoding in ("A", "B"):
            for key in designs:
                count = parts.get((encoding, key))
                if count is None:
                    continue
                table[count][0] += 1
                table[count][1] += outcome(key, encoding) == "perfect"
        print("\nrecovery by the number of synthesised fragments the design needed:")
        report["fragments"] = {}
        for count in sorted(table):
            picked, ok = table[count]
            print(f"   {count} fragment(s): {ok}/{picked} = {100 * ok / picked:.1f}%")
            report["fragments"][str(count)] = {
                "designs": picked,
                "recovered": ok,
                "percent": round(100 * ok / picked, 1),
            }

    out = args.out or args.run_dir / "figures" / "encoding_defects.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
