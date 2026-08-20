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

from nanopore3 import replicate  # noqa: E402
from nanopore3.figures import SINGLE_COLUMN, use_print_style  # noqa: E402

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
    # Arial throughout, black text and furniture, white only where text sits on
    # a filled segment.
    plt.rcParams.update(
        {
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "text.color": "black",
            "axes.labelcolor": "black",
            "axes.edgecolor": "black",
            "xtick.color": "black",
            "ytick.color": "black",
            "xtick.direction": "in",
            "ytick.direction": "in",
            # Keep SVG text as text so the vector stays editable in Illustrator
            # and the Arial the figure asks for is the Arial that renders.
            "svg.fonttype": "none",
            # figures.py sets savefig.bbox="tight" globally, which trims the
            # canvas: this figure was coming out 87.0 mm instead of the 89.0 mm
            # single column, and a journal rescaling it to fit would change the
            # effective type size - the one thing these point sizes exist to
            # control.  Margins are set by subplots_adjust below instead.
            "savefig.bbox": "standard",
        }
    )

    present = [
        category
        for category in replicate.CATEGORY_ORDER
        if any(r["summary"]["counts"][category] for r in results)
    ]
    # One row, as long as one row is legible.  Three short labels fit 89 mm at
    # 6 pt almost exactly, so the fit is measured rather than assumed: below 5 pt
    # the legend wraps to two rows instead of running off the page.
    labels = [replicate.CATEGORY_LABELS[c] for c in present]
    width_pt = SINGLE_COLUMN * 72
    per_entry = 2.4  # handle + padding + column gap, in font-size units
    required = sum(len(label) for label in labels) * 0.52 + len(present) * per_entry
    legend_size = min(6.0, width_pt / required)
    if legend_size < 5.0:
        columns = -(-len(present) // 2)
        required = (
            max(sum(len(label) for label in labels[:columns]), 1) * 0.52
            + columns * per_entry
        )
        legend_size = min(6.0, width_pt / required)
        legend_rows = 2
    else:
        columns = len(present)
        legend_rows = 1

    bars = 0.40 * len(results) + 0.16
    below = 0.42 + 0.18 * legend_rows
    height = bars + below
    figure, axes = plt.subplots(figsize=(SINGLE_COLUMN, height))

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
                edgecolor="black",
                linewidth=0.5,
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

        # Two ratios over two different units, each naming its own, centred in
        # the bar.  White because it sits on the filled segment.
        axes.text(
            50.0,
            position,
            f"{summary['wells_matched']}/{summary['wells_dedicated']} wells with "
            f"same content ({100 * summary['matched_fraction']:.1f}%)\n"
            f"{summary['clones_exact']}/{summary['clones_shared']} clones with "
            f"identical sequence ({100 * summary['exact_fraction']:.1f}%)",
            va="center",
            ha="center",
            fontsize=6.5,
            fontweight="bold",
            color="white",
            linespacing=1.35,
            zorder=3,
        )

    axes.set_yticks(positions)
    axes.set_yticklabels([r["label"] for r in results])
    axes.tick_params(axis="y", length=0)
    axes.tick_params(axis="x", direction="in", length=3, width=0.6, color="black")
    axes.set_xlim(0, 100)
    axes.set_xticks([0, 20, 40, 60, 80, 100])
    axes.set_xticklabels([f"{v}%" for v in (0, 20, 40, 60, 80, 100)])
    axes.set_xlabel("Outcome across wells with a sequenced clone")
    axes.set_ylim(-0.6, len(results) - 0.4)
    for side in ("top", "right"):
        axes.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        axes.spines[side].set_color("black")
        axes.spines[side].set_linewidth(0.6)

    handles = [
        Patch(facecolor=CATEGORY_COLOURS[c], edgecolor="black", linewidth=0.5)
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
        columnspacing=1.0,
        handletextpad=0.4,
        fontsize=legend_size,
        labelcolor="black",
    )
    # right < 1 leaves room for the "100%" tick label, which the tight bbox used
    # to absorb silently.
    figure.subplots_adjust(left=0.16, right=0.955, top=0.97, bottom=below / height)

    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path.with_suffix(".pdf"))
    figure.savefig(path.with_suffix(".png"), dpi=600)
    # Transparent so the SVG drops onto any figure background; the in-bar text is
    # white and only ever sits on a filled segment, so nothing vanishes.
    figure.savefig(path.with_suffix(".svg"), transparent=True)
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
    print(
        f"\nwrote {figure_path}, {figure_path.with_suffix('.png')}, "
        f"{figure_path.with_suffix('.svg')} (transparent),"
    )
    print(f"      {out_dir / f'{args.stem}_wells.csv'}, {out_dir / f'{args.stem}.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
