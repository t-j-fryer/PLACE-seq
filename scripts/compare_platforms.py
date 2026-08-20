#!/usr/bin/env python3
"""Nanopore clone picking against pooled Illumina, at matched sampling effort.

Three figures over the same library sets:

1. `sequences_per_well` - distinct sequences recovered per allocated well, which
   says how often a picked well held exactly one clone.
2. `reference_recovery` - percentage of each designed library recovered, split
   by the best outcome seen, nanopore against Illumina.
3. `read_accuracy` - identity to the design for nanopore reads, nanopore
   consensuses and Illumina reads.

Illumina is subsampled per block to the number of wells the nanopore run
sequenced of that block, so both platforms spend the same sampling effort: N
colonies picked from a block against N reads taken of it.  Without that, the
comparison only says that a deep pool sees more designs than 96 colonies.

    python scripts/compare_platforms.py \
        --run-dir runs/260608-AI-DBTL-v4 \
        --illumina <.../LAB_STUFFER_SUMO_CONCORDANT>

Outputs into <run-dir>/figures unless --out-dir says otherwise.
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

from nanopore3 import platforms  # noqa: E402
from nanopore3.figures import MM, save_figure, use_journal_style  # noqa: E402

# Widths chosen to fit the content rather than a column: wide enough to read at
# 6 pt, with no blank margin to spare.
NARROW = 68 * MM
COMPACT = 114 * MM

# Which barcode and culture plates hold which library.  This is knowledge about
# the experiment, not about the data: RP08 carried the stuffer library, RP05
# carried both SUMO encodings across 22 pooled culture plates.  The A+B set is
# not a fourth library - it is the same 342 designs counted as recovered if
# either encoding was.
LIBRARIES = [
    platforms.LibrarySet(
        key="stuffer",
        label="Library 1 (stuffer)",
        library_id="lab",
        plate_id="RP08",
        culture_plate_prefix="LAB_",
        encodings=("",),
        illumina_groups=("LAB_STUFFER",),
    ),
    platforms.LibrarySet(
        key="sumo_a",
        label="Library 3 (A)",
        library_id="sumo_ab",
        plate_id="RP05",
        culture_plate_prefix="SUMO_A_",
        encodings=("A",),
        illumina_groups=("SUMO_A",),
    ),
    platforms.LibrarySet(
        key="sumo_b",
        label="Library 3 (B)",
        library_id="sumo_ab",
        plate_id="RP05",
        culture_plate_prefix="SUMO_B_",
        encodings=("B",),
        illumina_groups=("SUMO_B",),
    ),
]
UNION = platforms.LibrarySet(
    key="sumo_ab",
    label="Library 3 (A+B)",
    library_id="sumo_ab",
    plate_id="RP05",
    culture_plate_prefix="SUMO_",
    encodings=("A", "B"),
    illumina_groups=("SUMO_A", "SUMO_B"),
)

# Grey, blue, orange as in the source figures, plus green for the union set.
#
#   scripts/validate_palette.py "#4D4D4D,#0072B2,#D55E00,#009E73" --pairs all
#
# gives worst normal-vision dE 17.2 (floor 15) and worst CVD dE 11.0 (target 8),
# so all four stay separable.  The grey trips the per-colour chroma policy by
# construction - it is a neutral standing in for the stuffer library.
COLOURS = {
    "stuffer": "#4D4D4D",
    "sumo_a": "#0072B2",
    "sumo_b": "#D55E00",
    "sumo_ab": "#009E73",
}
ILLUMINA_HATCH = "///"


def grouped_bars(
    axes,
    categories: list[str],
    series: list[tuple[str, str, list[float], dict]],
    *,
    width: float = 0.8,
) -> None:
    """Draw grouped bars, each series a (key, label, values, kwargs) tuple."""

    count = len(series)
    bar = width / count
    for index, (_key, _label, values, style) in enumerate(series):
        offsets = [
            position - width / 2 + bar * (index + 0.5) for position in range(len(categories))
        ]
        kwargs = {"edgecolor": "black", "linewidth": 0.5, "zorder": 2} | style
        axes.bar(offsets, values, width=bar * 0.92, **kwargs)
    axes.set_xticks(range(len(categories)))
    axes.set_xticklabels(categories)
    axes.tick_params(axis="x", length=0)
    axes.tick_params(axis="y", direction="in", length=3, width=0.6)
    for side in ("top", "right"):
        axes.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        axes.spines[side].set_linewidth(0.6)


def two_panel(figure_width: float, height: float, ratios: tuple[float, float]):
    """Two panels sharing a y scale is how the small categories stay readable.

    On one axis the Perfect bars run to 100% and everything else is a smear along
    the floor; splitting the axis spends the width on data instead of blank space.
    """

    figure, axes = plt.subplots(
        1, 2, figsize=(figure_width, height), gridspec_kw={"width_ratios": list(ratios)}
    )
    return figure, axes


def sequences_per_well_figure(histograms: dict[str, dict[int, int]], allocated: int, path: Path):
    """How many distinct sequences each allocated well yielded."""

    use_journal_style()
    figure, axes = plt.subplots(figsize=(NARROW, 1.55))
    buckets = [0, 1, 2, 3]
    series = [
        (
            library.key,
            library.label,
            [100.0 * histograms[library.key].get(b, 0) / allocated for b in buckets],
            {"color": COLOURS[library.key], "label": library.label},
        )
        for library in LIBRARIES
    ]
    grouped_bars(axes, ["0", "1", "2", "3+"], series)
    axes.set_ylim(0, 100)
    axes.set_yticks([0, 25, 50, 75, 100])
    axes.set_xlabel("Distinct sequences per allocated well")
    axes.set_ylabel("Block-allocated wells (%)")
    axes.legend(
        loc="upper right",
        frameon=False,
        fontsize=6,
        handlelength=1.0,
        handleheight=0.85,
        handletextpad=0.4,
        title=f"n = {allocated:,} wells each",
        title_fontsize=6,
        alignment="left",
        borderpad=0.1,
        labelspacing=0.35,
    )
    figure.subplots_adjust(left=0.155, right=0.995, top=0.97, bottom=0.235)
    return save_figure(figure, path)


def _nice_limit(largest: float) -> float:
    """A round ceiling just above the tallest bar, so nothing is clipped."""

    for step in (0.25, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 7.5, 10, 12.5, 20, 25):
        if largest <= step * 4:
            return step * 4
    return 100.0


def _outcome_panels(
    axes,
    values_for,
    *,
    ylabel: str,
) -> None:
    """Perfect on a full axis, the other two outcomes on a magnified one.

    The magnified ceiling comes from the data.  A fixed one silently clipped a
    bar the first time this ran, which is the sort of error that survives review
    because a clipped bar still looks like a bar.
    """

    sets = [*LIBRARIES, UNION]
    small_limit = _nice_limit(
        max(
            values_for(library, platform, name)
            for library in sets
            for platform in ("nanopore", "illumina")
            for name in platforms.CLASS_ORDER[1:]
        )
    )
    panels = (
        (axes[0], [PERFECT_LABEL], (0, 100), [0, 25, 50, 75, 100]),
        (
            axes[1],
            list(platforms.CLASS_ORDER[1:]),
            (0, small_limit),
            [small_limit * f for f in (0, 0.25, 0.5, 0.75, 1.0)],
        ),
    )
    for panel, categories, limits, ticks in panels:
        series = []
        for library in sets:
            for platform, hatch in (("nanopore", None), ("illumina", ILLUMINA_HATCH)):
                style = {"color": COLOURS[library.key]}
                if hatch:
                    style |= {
                        "hatch": hatch,
                        "color": "white",
                        "edgecolor": COLOURS[library.key],
                    }
                series.append(
                    (
                        library.key,
                        library.label,
                        [values_for(library, platform, name) for name in categories],
                        style,
                    )
                )
        grouped_bars(panel, categories, series, width=0.88)
        panel.set_ylim(*limits)
        panel.set_yticks(ticks)
        if all(float(value).is_integer() for value in ticks):
            panel.set_yticklabels([f"{int(value)}" for value in ticks])
    axes[0].set_ylabel(ylabel)
    axes[1].tick_params(axis="y", labelsize=6)


PERFECT_LABEL = platforms.CLASS_ORDER[0]


def _outcome_legend(figure, *, rows_at: tuple[float, float]) -> None:
    """Libraries on one row, platforms on the next, both under the panels."""

    sets = [*LIBRARIES, UNION]
    figure.legend(
        handles=[
            Patch(facecolor=COLOURS[s.key], edgecolor="black", linewidth=0.5, label=s.label)
            for s in sets
        ],
        loc="lower center",
        bbox_to_anchor=(0.5, rows_at[0]),
        ncol=len(sets),
        frameon=False,
        fontsize=6,
        handlelength=1.0,
        handleheight=0.85,
        handletextpad=0.4,
        columnspacing=1.0,
    )
    figure.legend(
        handles=[
            Patch(
                facecolor="#8C8C8C",
                edgecolor="black",
                linewidth=0.5,
                label="Nanopore, picked clones",
            ),
            Patch(
                facecolor="white",
                edgecolor="black",
                hatch=ILLUMINA_HATCH,
                linewidth=0.5,
                label="Illumina, block-matched depth",
            ),
        ],
        loc="lower center",
        bbox_to_anchor=(0.5, rows_at[1]),
        ncol=2,
        frameon=False,
        fontsize=6,
        handlelength=1.0,
        handleheight=0.85,
        handletextpad=0.4,
        columnspacing=1.2,
    )


def reference_recovery_figure(recovery: dict[str, dict[str, dict[str, float]]], path: Path):
    """Percentage of each designed library recovered, by platform."""

    use_journal_style()
    figure, axes = two_panel(COMPACT, 1.95, (1.0, 1.9))
    _outcome_panels(
        axes,
        lambda library, platform, name: recovery[library.key][platform][name],
        ylabel="Recovered references (%)",
    )
    figure.subplots_adjust(left=0.115, right=0.995, top=0.97, bottom=0.27, wspace=0.26)
    _outcome_legend(figure, rows_at=(0.10, 0.005))
    return save_figure(figure, path)


def population_figure(population: dict[str, dict[str, dict[str, float]]], path: Path):
    """Where the sequences themselves fall: consensuses against reads."""

    use_journal_style()
    figure, axes = two_panel(COMPACT, 1.95, (1.0, 1.9))
    _outcome_panels(
        axes,
        lambda library, platform, name: population[library.key][platform][name],
        ylabel="Sequences recovered (%)",
    )
    figure.subplots_adjust(left=0.115, right=0.995, top=0.97, bottom=0.27, wspace=0.26)
    _outcome_legend(figure, rows_at=(0.10, 0.005))
    return save_figure(figure, path)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--illumina", type=Path, required=True, help="Illumina results directory")
    parser.add_argument("--out-dir", type=Path, help="default <run-dir>/figures")
    parser.add_argument(
        "--wells-per-plate",
        type=int,
        default=95,
        help="allocated wells per culture plate; H12 was not picked in this run",
    )
    parser.add_argument("--culture-plates", type=int, default=11, help="per library set")
    parser.add_argument("--seed", type=int, default=141142, help="Illumina subsampling seed")
    args = parser.parse_args()

    out_dir = args.out_dir or args.run_dir / "figures"
    config = json.loads((args.run_dir / "run.json").read_text(encoding="utf-8"))["config"]

    universes: dict[str, set] = {}
    for library in [*LIBRARIES, UNION]:
        paths = [Path(p) for p in config["reference_libraries"][library.library_id]["fasta"]]
        universes[library.key] = platforms.load_design_universe(paths)

    clones = platforms.load_clones(args.run_dir)
    illumina = platforms.load_illumina_reads(args.illumina)
    print(
        f"loaded {sum(len(v) for v in clones.values()):,} clones "
        f"and {len(illumina):,} Illumina reads"
    )

    # Matched depth: allocated wells per block becomes reads drawn per block.
    # Allocated, not recovered - the effort spent is the colonies picked, and
    # charging nanopore's failures to Illumina as reduced depth would flatter us.
    depths: dict[tuple[str, int], int] = {}
    recovered: dict[tuple[str, int], int] = {}
    for library in LIBRARIES:
        set_clones = platforms.clones_in(clones, library)
        allocated_blocks = platforms.allocated_wells_per_block(
            set_clones, wells_per_plate=args.wells_per_plate
        )
        for block, wells in allocated_blocks.items():
            depths[(library.illumina_groups[0], block)] = wells
        for block, wells in platforms.wells_sequenced_per_block(set_clones).items():
            recovered[(library.illumina_groups[0], block)] = wells
    sample = platforms.subsample_by_block(illumina, depths, seed=args.seed)
    print(
        f"matched depth: {sum(depths.values()):,} Illumina reads "
        f"sampled across {len(depths)} blocks"
    )

    allocated = args.culture_plates * args.wells_per_plate
    histograms: dict[str, dict[int, int]] = {}
    for library in LIBRARIES:
        histogram = platforms.distinct_sequences_per_well(
            platforms.clones_in(clones, library),
            culture_plates=args.culture_plates,
            wells_per_plate=args.wells_per_plate,
        )
        histograms[library.key] = dict(histogram)
        print(f"  {library.label}: " + "  ".join(
            f"{k}:{100 * v / allocated:.1f}%" for k, v in sorted(histogram.items())
        ))

    recovery: dict[str, dict[str, dict[str, float]]] = {}
    allowed_by_set: dict[str, set] = {}
    detail_rows: list[dict[str, object]] = []
    for library in [*LIBRARIES, UNION]:
        designs = platforms.design_keys(universes[library.key], library.encodings)
        nano = platforms.nanopore_recovery(platforms.clones_in(clones, library))
        allowed = {
            (encoding, block, key)
            for encoding, block, key in universes[library.key]
            if encoding in library.encodings
        }
        allowed_by_set[library.key] = allowed
        ill = platforms.illumina_recovery(sample, allowed)
        full = platforms.illumina_recovery(illumina, allowed)
        recovery[library.key] = {
            "designs": len(designs),
            "nanopore": platforms.recovery_fractions(nano, len(designs)),
            "illumina": platforms.recovery_fractions(ill, len(designs)),
            "illumina_full_depth": platforms.recovery_fractions(full, len(designs)),
        }
        for name in platforms.CLASS_ORDER:
            detail_rows.append(
                {
                    "library_set": library.label,
                    "designs": len(designs),
                    "outcome": name,
                    "nanopore_percent": round(recovery[library.key]["nanopore"][name], 2),
                    "illumina_matched_percent": round(recovery[library.key]["illumina"][name], 2),
                    "illumina_full_depth_percent": round(
                        recovery[library.key]["illumina_full_depth"][name], 2
                    ),
                }
            )
        print(
            f"  {library.label}: {len(designs)} designs; "
            f"nanopore perfect {recovery[library.key]['nanopore']['Perfect']:.1f}%, "
            f"Illumina matched {recovery[library.key]['illumina']['Perfect']:.1f}%, "
            f"Illumina full depth {recovery[library.key]['illumina_full_depth']['Perfect']:.1f}%"
        )

    population: dict[str, dict[str, dict[str, float]]] = {}
    for library in [*LIBRARIES, UNION]:
        nano = platforms.nanopore_population(platforms.clones_in(clones, library))
        ill = platforms.illumina_population(sample, allowed_by_set[library.key])
        population[library.key] = {
            "nanopore": platforms.population_fractions(nano),
            "illumina": platforms.population_fractions(ill),
            "nanopore_counts": dict(nano),
            "illumina_counts": dict(ill),
        }
        print(
            f"  {library.label}: consensuses "
            + ", ".join(
                f"{n} {population[library.key]['nanopore'][n]:.1f}%"
                for n in platforms.CLASS_ORDER
            )
            + " | Illumina reads "
            + ", ".join(
                f"{n} {population[library.key]['illumina'][n]:.1f}%"
                for n in platforms.CLASS_ORDER
            )
        )

    read_identities = platforms.load_nanopore_read_identities(args.run_dir, [*LIBRARIES, UNION])
    consensus_identities = platforms.load_consensus_identities(args.run_dir, [*LIBRARIES, UNION])
    accuracy: dict[str, dict[str, dict[str, float]]] = {}
    for library in [*LIBRARIES, UNION]:
        accuracy[library.key] = {
            "nanopore_read": platforms.summarise_identities(read_identities[library.key]),
            "nanopore_consensus": platforms.summarise_identities(
                consensus_identities[library.key]
            ),
            "illumina_read": platforms.summarise_identities(
                platforms.illumina_read_identities(sample, allowed_by_set[library.key])
            ),
            "illumina_read_full_depth": platforms.summarise_identities(
                platforms.illumina_read_identities(illumina, allowed_by_set[library.key])
            ),
        }
        row = accuracy[library.key]
        print(
            f"  {library.label}: nanopore read {100 * row['nanopore_read']['mean']:.2f}%, "
            f"consensus {100 * row['nanopore_consensus']['mean']:.2f}%, "
            f"Illumina read {100 * row['illumina_read']['mean']:.2f}%"
        )

    out_dir.mkdir(parents=True, exist_ok=True)
    with (out_dir / "platform_reference_recovery.csv").open("w", encoding="utf-8", newline="") as h:
        writer = csv.DictWriter(h, fieldnames=list(detail_rows[0]))
        writer.writeheader()
        writer.writerows(detail_rows)
    (out_dir / "platform_comparison.json").write_text(
        json.dumps(
            {
                "run_dir": str(args.run_dir),
                "illumina_results": str(args.illumina),
                "seed": args.seed,
                "allocated_wells_per_set": allocated,
                "wells_per_plate": args.wells_per_plate,
                "matched_depth_per_block": {f"{g}:{b}": n for (g, b), n in sorted(depths.items())},
                "recovered_wells_per_block": {
                    f"{g}:{b}": n for (g, b), n in sorted(recovered.items())
                },
                "sequences_per_well": histograms,
                "reference_recovery": recovery,
                "read_accuracy": accuracy,
                "populations": population,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    paths = [
        sequences_per_well_figure(histograms, allocated, out_dir / "platform_sequences_per_well"),
        reference_recovery_figure(recovery, out_dir / "platform_reference_recovery"),
        population_figure(population, out_dir / "platform_sequence_populations"),
    ]
    print("\nwrote:")
    for path in paths:
        print(f"  {path.with_suffix('')}.{{pdf,png,svg}}")
    print(f"  {out_dir / 'platform_reference_recovery.csv'}")
    print(f"  {out_dir / 'platform_comparison.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
