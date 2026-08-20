#!/usr/bin/env python3
"""Compare dedicated plate barcodes against the same culture plates in a pool.

Two culture plates in the 260608 run were sequenced twice: once on a barcode of
their own (RP06, RP07) and once inside the compressed RP05 pool, where 22 culture
plates share every well.  Comparing them measures demultiplexing, assignment,
consensus and compressed-PCR deconvolution against independent data.

Reproduce the published figure with:

    python scripts/compare_replicate_plates.py \
        --run-dir runs/260608-AI-DBTL-v4 \
        --pair RP06=RP05/SUMO_A_P1 \
        --pair RP07=RP05/SUMO_A_P2

Each pair is verified against the run's own block-to-culture-plate map before it
is used, so a mistaken pairing fails loudly instead of producing a plausible
figure.  Outputs land in <run-dir>/figures unless --out-dir says otherwise.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.transforms import blended_transform_factory  # noqa: E402

from nanopore3 import replicate  # noqa: E402
from nanopore3.figures import (  # noqa: E402
    INK,
    MUTED,
    SINGLE_COLUMN,
    use_print_style,
)

# Teal for agreement, warm hues for the two ways a well can be seen only once,
# brick for outright disagreement, neutral for "no data".
#
#   scripts/validate_palette.py "#2E6F6A,#6FA9A0,#A8442E,#E0A33A" --pairs all
#
# reports worst normal-vision dE 17.8 (floor 15) and worst CVD dE 9.5 (target 8),
# so every pair stays separable including under protanopia and deuteranopia.  The
# per-colour chroma and contrast policies are relaxed here on purpose: these are
# large area fills carrying a legend and white in-bar text, not thin line series,
# and the neutral is exempt by the convention in figures.py.  An earlier draft
# copied the source figure's pale teal against pale grey, which the validator
# rejected at 13.6 dE - the two thin slivers really are hard to tell apart.
CATEGORY_COLOURS = {
    replicate.MATCHED_EXACT: "#2E6F6A",
    replicate.MATCHED_DIFFERENT_SEQUENCE: "#6FA9A0",
    replicate.DIFFERENT_DESIGN: "#A8442E",
    replicate.POOL_ONLY: "#E0A33A",
    replicate.DEDICATED_ONLY: "#EEF1F3",
}


def parse_pair(text: str) -> tuple[str, str, str]:
    """`RP06=RP05/SUMO_A_P1` -> (dedicated, pooled barcode, culture plate)."""

    dedicated, _, pooled = text.partition("=")
    barcode, _, culture_plate = pooled.partition("/")
    if not (dedicated and barcode and culture_plate):
        raise argparse.ArgumentTypeError(
            f"expected DEDICATED=POOLED/CULTURE_PLATE, got {text!r}"
        )
    return dedicated, barcode, culture_plate


def write_wells_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = [
        "label",
        "dedicated_plate",
        "pooled_plate",
        "culture_plate",
        "well_id",
        "category",
        "dedicated_clones",
        "pooled_clones",
        "shared_clones",
        "exact_clones",
        "differing_clones",
        "edit_distances",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def concordance_figure(results: list[dict[str, object]], path: Path) -> Path:
    """Stacked outcome bars, one per plate pair.

    Each bar spans every well with a sequenced clone in either dataset, so a well
    only the pool recovered widens the bar rather than silently vanishing.
    """

    use_print_style()
    present = [
        category
        for category in replicate.CATEGORY_ORDER
        if any(r["summary"]["counts"][category] for r in results)
    ]
    columns = 2 if len(present) <= 4 else 3
    legend_rows = -(-len(present) // columns)
    # Height is built from the parts rather than guessed, so adding a plate pair
    # or a category does not push the legend into the axis label.
    bars = 0.40 * len(results) + 0.16
    below = 0.42 + 0.18 * legend_rows
    height = bars + below
    figure, axes = plt.subplots(figsize=(SINGLE_COLUMN, height))
    label_transform = blended_transform_factory(axes.transAxes, axes.transData)

    positions = list(range(len(results)))[::-1]
    for position, result in zip(positions, results, strict=True):
        summary = result["summary"]
        counts = summary["counts"]
        total = summary["wells_total"] or 1
        left = 0.0
        for category in replicate.CATEGORY_ORDER:
            count = counts[category]
            if not count:
                continue
            width = 100.0 * count / total
            axes.barh(
                position,
                width,
                left=left,
                height=0.60,
                color=CATEGORY_COLOURS[category],
                edgecolor="white",
                linewidth=0.4,
                zorder=2,
            )
            # Name a minority segment in place where it is wide enough to read;
            # the legend carries the rest.
            if category != replicate.MATCHED_EXACT and width >= 7.0:
                axes.text(
                    left + width / 2,
                    position,
                    f"+{count}",
                    va="center",
                    ha="center",
                    fontsize=6,
                    fontweight="bold",
                    color="white",
                    zorder=3,
                )
            left += width

        axes.text(
            1.8,
            position,
            f"{summary['wells_matched']}/{summary['wells_dedicated']} matched "
            f"({100 * summary['matched_fraction']:.1f}%)\n"
            f"{summary['clones_exact']}/{summary['clones_shared']} exact seq. "
            f"({100 * summary['exact_fraction']:.1f}%)",
            va="center",
            ha="left",
            fontsize=6.5,
            fontweight="bold",
            color="white",
            linespacing=1.35,
            zorder=3,
        )

        # Bar labels are drawn by hand so the plate name can carry weight while
        # the barcode-to-culture-plate mapping stays recessive but present.
        axes.text(
            -0.02,
            position + 0.11,
            result["label"],
            transform=label_transform,
            ha="right",
            va="center",
            fontsize=7,
            color=INK,
        )
        axes.text(
            -0.02,
            position - 0.14,
            f"{result['dedicated_plate']} \u00b7 {result['culture_plate']}",
            transform=label_transform,
            ha="right",
            va="center",
            fontsize=5.8,
            color=MUTED,
        )

    axes.set_yticks([])
    axes.set_xlim(0, 100)
    axes.set_xticks([0, 20, 40, 60, 80, 100])
    axes.set_xticklabels([f"{v}%" for v in (0, 20, 40, 60, 80, 100)])
    axes.set_xlabel("Outcome across wells with a sequenced clone")
    axes.set_ylim(-0.6, len(results) - 0.4)
    for side in ("top", "right", "left"):
        axes.spines[side].set_visible(False)
    axes.spines["bottom"].set_color(MUTED)

    handles = [
        Patch(facecolor=CATEGORY_COLOURS[c], edgecolor="white", linewidth=0.4)
        for c in present
    ]
    figure.legend(
        handles,
        [replicate.CATEGORY_LABELS[c] for c in present],
        loc="lower center",
        bbox_to_anchor=(0.5, 0.0),
        ncol=columns,
        frameon=False,
        handlelength=1.0,
        handleheight=0.85,
        columnspacing=1.1,
        fontsize=6.5,
        labelcolor=INK,
    )
    figure.subplots_adjust(left=0.30, right=0.99, top=0.98, bottom=below / height)

    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path.with_suffix(".pdf"))
    figure.savefig(path.with_suffix(".png"), dpi=600)
    plt.close(figure)
    return path.with_suffix(".pdf")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument(
        "--pair",
        action="append",
        required=True,
        metavar="DEDICATED=POOLED/CULTURE_PLATE",
        help="repeatable, e.g. RP06=RP05/SUMO_A_P1",
    )
    parser.add_argument("--label", action="append", help="bar label per pair, in order")
    parser.add_argument("--out-dir", type=Path, help="default <run-dir>/figures")
    parser.add_argument("--stem", default="replicate_concordance")
    parser.add_argument(
        "--no-chimeras",
        action="store_true",
        help="exclude chimeric clones (they are included by default: real sequence in a real well)",
    )
    parser.add_argument(
        "--allow-unconfirmed-pairing",
        action="store_true",
        help="proceed even if the dedicated barcode's blocks do not match the named culture plate",
    )
    args = parser.parse_args()

    out_dir = args.out_dir or args.run_dir / "figures"
    pairs = [parse_pair(text) for text in args.pair]
    labels = args.label or [f"Plate {i + 1}" for i in range(len(pairs))]
    if len(labels) != len(pairs):
        parser.error(f"{len(labels)} labels for {len(pairs)} pairs")

    results: list[dict[str, object]] = []
    well_rows: list[dict[str, object]] = []
    failed = False

    for (dedicated, barcode, culture_plate), label in zip(pairs, labels, strict=True):
        check = replicate.verify_pairing(args.run_dir, dedicated, culture_plate)
        state = "confirmed" if check["confirmed"] else "NOT CONFIRMED"
        print(
            f"{label}: {dedicated} vs {barcode}/{culture_plate} - pairing {state} "
            f"({check['clones_consistent']}/{check['clones_examined']} clones, "
            f"sources seen {check['sources_seen']})"
        )
        if not check["confirmed"] and not args.allow_unconfirmed_pairing:
            failed = True
            continue

        left = replicate.load_calls(
            args.run_dir, dedicated, include_chimeras=not args.no_chimeras
        )
        right = replicate.load_calls(
            args.run_dir,
            barcode,
            culture_plate=culture_plate,
            include_chimeras=not args.no_chimeras,
        )
        outcomes = replicate.compare_wells(left, right)
        summary = replicate.summarise(outcomes)
        results.append(
            {
                "label": label,
                "dedicated_plate": dedicated,
                "pooled_plate": barcode,
                "culture_plate": culture_plate,
                "pairing_check": check,
                "summary": summary,
            }
        )
        print(
            f"  {summary['wells_matched']}/{summary['wells_dedicated']} dedicated wells matched; "
            f"{summary['clones_exact']}/{summary['clones_shared']} shared clones byte-identical; "
            f"{summary['counts'][replicate.POOL_ONLY]} pool-only, "
            f"{summary['counts'][replicate.DEDICATED_ONLY]} dedicated-only, "
            f"{summary['chimeras_matched']} chimeric clones matched"
        )
        for outcome in outcomes:
            well_rows.append(
                {
                    "label": label,
                    "dedicated_plate": dedicated,
                    "pooled_plate": barcode,
                    "culture_plate": culture_plate,
                    "well_id": outcome.well_id,
                    "category": outcome.category,
                    "dedicated_clones": ";".join(outcome.dedicated),
                    "pooled_clones": ";".join(outcome.pooled),
                    "shared_clones": ";".join(outcome.shared),
                    "exact_clones": ";".join(outcome.exact),
                    "differing_clones": ";".join(outcome.differing),
                    "edit_distances": ";".join(f"{k}={v}" for k, v in outcome.edit_distances),
                }
            )

    if failed:
        print(
            "\nRefusing to plot an unverified pairing; "
            "pass --allow-unconfirmed-pairing to override."
        )
        return 1
    if not results:
        print("No pairs to plot.")
        return 1

    write_wells_csv(out_dir / f"{args.stem}_wells.csv", well_rows)
    (out_dir / f"{args.stem}.json").write_text(
        json.dumps(
            {
                "run_dir": str(args.run_dir),
                "include_chimeras": not args.no_chimeras,
                "pairs": results,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    figure_path = concordance_figure(results, out_dir / args.stem)
    print(f"\nwrote {figure_path}, {figure_path.with_suffix('.png')},")
    print(f"      {out_dir / f'{args.stem}_wells.csv'}, {out_dir / f'{args.stem}.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
