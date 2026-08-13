"""Publication figures for plate occupancy, clonality, and deconvolution.

Figures are sized and styled for a print journal: single column 89 mm, double
column 183 mm, sans-serif at 6-8 pt, vector PDF plus a high-resolution PNG.

The categorical palette is validated, not chosen by eye.  ``#0072B2``,
``#D55E00`` and ``#009E73`` clear the all-pairs colour-vision-deficiency check
(worst OKLab dE 11.0 under protanopia/deuteranopia, 18.7 unsimulated); a fourth
hue drops the worst pair to 7.6, below target, so three is the cap for anything
where two marks can touch.  Re-check with ``scripts/validate_palette.py``.

Matplotlib is an optional dependency: import this module only when it is needed.
"""

from __future__ import annotations

import csv
import gzip
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

# Validated categorical slots, assigned in this fixed order and never cycled.
SERIES = ("#0072B2", "#D55E00", "#009E73")
# Neutrals carry "no category" and are exempt from the chroma floor by design.
INK, MUTED, FAINT = "#1B1F24", "#6E7781", "#D8DEE4"
EMPTY_WELL = "#F2F4F6"

MM = 1 / 25.4
SINGLE_COLUMN, DOUBLE_COLUMN = 89 * MM, 183 * MM

ROWS = "ABCDEFGH"
COLUMNS = tuple(range(1, 13))

_SEQUENTIAL = LinearSegmentedColormap.from_list(
    "nanopore3_blue", ["#E8F1F8", "#9CC6E2", "#4E97C7", "#0072B2", "#00456B"]
)


def use_print_style() -> None:
    """Apply print-journal defaults: thin recessive furniture, small sans text."""

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
            "font.size": 7,
            "axes.titlesize": 8,
            "axes.labelsize": 7,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "legend.fontsize": 6.5,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "axes.edgecolor": MUTED,
            "axes.labelcolor": INK,
            "text.color": INK,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 150,
            "savefig.dpi": 600,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.02,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


@dataclass(frozen=True, slots=True)
class PlateSummary:
    """Per-well observations for one colony PCR plate."""

    plate_id: str
    clones_per_well: Mapping[str, int]
    plates_per_well: Mapping[str, tuple[str, ...]]
    expected_clones: int | None


def _read_csv_gz(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def summarize_run(
    run_dir: Path, expected_clones: Mapping[str, int] | None = None
) -> tuple[list[PlateSummary], Counter[str], Counter[tuple[str, str]]]:
    """Derive plate occupancy and deconvolution outcome from a completed run."""

    rows = _read_csv_gz(run_dir / "stages" / "03_assignment" / "assignment_calls.csv.gz")
    assigned = [row for row in rows if row["assignment_status"].startswith("assigned")]
    clones: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    sources: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    for row in assigned:
        clones[row["plate_id"]][row["well_id"]].add(row["reference_ids"])
        if row.get("culture_plate"):
            sources[row["plate_id"]][row["well_id"]].add(row["culture_plate"])
    status = Counter(row.get("culture_plate_status", "not_configured") for row in rows)
    by_plate = Counter(
        (row["plate_id"], row.get("culture_plate_status", "not_configured"))
        for row in rows
        if row["plate_id"]
    )
    summaries = [
        PlateSummary(
            plate_id=plate_id,
            clones_per_well={well: len(genes) for well, genes in wells.items()},
            plates_per_well={
                well: tuple(sorted(found)) for well, found in sources[plate_id].items()
            },
            expected_clones=(expected_clones or {}).get(plate_id),
        )
        for plate_id, wells in sorted(clones.items())
    ]
    return summaries, status, by_plate


def _plate_grid(
    ax,
    values: Mapping[str, float],
    vmax: float,
    title: str,
    *,
    show_columns: bool = True,
    show_rows: bool = True,
) -> None:
    ax.set_xlim(0.5, 12.5)
    ax.set_ylim(0.5, 8.5)
    ax.invert_yaxis()
    ax.set_aspect("equal")
    for spine in ax.spines.values():
        spine.set_visible(False)
    for row_index, row in enumerate(ROWS, start=1):
        for column in COLUMNS:
            well = f"{row}{column}"
            value = values.get(well)
            colour = EMPTY_WELL if value is None else _SEQUENTIAL(value / vmax if vmax else 0)
            # A 2px-equivalent surface gap keeps adjacent cells legible.
            ax.add_patch(
                Rectangle(
                    (column - 0.44, row_index - 0.44), 0.88, 0.88,
                    facecolor=colour, edgecolor="white", linewidth=0.5,
                )
            )
    # Plate coordinates repeat in every panel, so they are drawn only on the
    # outer edge; repeating them 96 times per panel is noise, not information.
    shown = [c for c in COLUMNS if c % 2 == 1]
    ax.set_xticks(shown)
    ax.set_xticklabels([str(c) for c in shown] if show_columns else [])
    ax.set_yticks(range(1, 9))
    ax.set_yticklabels(list(ROWS) if show_rows else [])
    ax.tick_params(length=0, pad=1.5)
    ax.set_title(title, pad=3, loc="left", color=INK, fontsize=6.5)


def plate_occupancy_figure(summaries: Sequence[PlateSummary], path: Path) -> Path:
    """Small multiples: distinct genes recovered per well, one panel per plate."""

    use_print_style()
    count = len(summaries)
    columns = min(4, count)
    rows = (count + columns - 1) // columns
    # A 96-well plate is 12x8, so the panel aspect fixes the figure height once
    # the column count is chosen; deriving it avoids acres of dead white space.
    panel_width = DOUBLE_COLUMN / columns
    figure = plt.figure(figsize=(DOUBLE_COLUMN, rows * panel_width * (8 / 12) * 1.34 + 0.30))
    grid = figure.add_gridspec(rows, columns, hspace=0.34, wspace=0.12)

    vmax = max(
        (max(summary.clones_per_well.values(), default=0) for summary in summaries),
        default=1,
    ) or 1
    bottom_of_column: dict[int, int] = {}
    for index in range(count):
        bottom_of_column[index % columns] = index

    for index, summary in enumerate(summaries):
        ax = figure.add_subplot(grid[index // columns, index % columns])
        observed = list(summary.clones_per_well.values())
        median = sorted(observed)[len(observed) // 2] if observed else 0
        title = f"{summary.plate_id}  ·  {len(observed)} wells · median {median}"
        if summary.expected_clones is not None:
            title += f" · expected {summary.expected_clones}"
        _plate_grid(
            ax, summary.clones_per_well, vmax, title,
            show_columns=index == bottom_of_column[index % columns],
            show_rows=index % columns == 0,
        )

    mappable = plt.cm.ScalarMappable(cmap=_SEQUENTIAL)
    mappable.set_clim(0, vmax)
    if count < rows * columns:
        # Spend the empty panel on the scale rather than leaving a hole.
        spare = figure.add_subplot(grid[count // columns, count % columns])
        spare.axis("off")
        box = spare.get_position()
        bar_axes = figure.add_axes(
            (box.x0 + 0.012, box.y0 + box.height * 0.42, box.width * 0.80, 0.014)
        )
    else:
        bar_axes = figure.add_axes((0.13, 0.02, 0.30, 0.014))
    bar = figure.colorbar(mappable, cax=bar_axes, orientation="horizontal")
    bar.set_label("distinct genes per well", color=INK, labelpad=3, fontsize=6.5)
    bar.outline.set_visible(False)
    bar.ax.tick_params(length=2, width=0.6, pad=1.5, labelsize=6)
    return _save(figure, path, tight=True)


def clonality_figure(summaries: Sequence[PlateSummary], path: Path) -> Path:
    """Observed clones per well by plate, against the expected count where declared."""

    use_print_style()
    figure, ax = plt.subplots(figsize=(DOUBLE_COLUMN * 0.62, 0.34 * len(summaries) + 0.95))
    labels = []
    any_expected = False
    for position, summary in enumerate(summaries):
        observed = sorted(summary.clones_per_well.values())
        labels.append(summary.plate_id)
        if not observed:
            continue
        # Spread ties symmetrically about the row so the pile-up at each integer
        # shows its size, rather than hiding it under one opaque dot.
        offsets: list[float] = []
        for value in sorted(set(observed)):
            tied = [v for v in observed if v == value]
            span = min(0.30, 0.05 * (len(tied) - 1))
            step = 0 if len(tied) == 1 else (2 * span) / (len(tied) - 1)
            offsets.extend(-span + step * i for i in range(len(tied)))
        ordered = sorted(observed)
        ax.scatter(
            ordered, [position + o for o in offsets],
            s=3.0, linewidths=0, color=SERIES[0], alpha=0.55, zorder=2,
        )
        median = ordered[len(ordered) // 2]
        ax.plot(
            [median, median], [position - 0.36, position + 0.36],
            color=INK, linewidth=1.3, zorder=4, solid_capstyle="butt",
        )
        if summary.expected_clones is not None:
            any_expected = True
            ax.plot(
                [summary.expected_clones] * 2, [position - 0.40, position + 0.40],
                color=SERIES[1], linewidth=1.3, zorder=5,
                linestyle=(0, (1.6, 1.2)), solid_capstyle="butt",
            )
    ax.set_yticks(range(len(summaries)))
    ax.set_yticklabels(labels)
    ax.set_ylim(-0.7, len(summaries) - 0.3)
    ax.invert_yaxis()
    ax.set_xlabel("distinct genes per well")
    ax.set_xlim(left=0)
    ax.xaxis.grid(True, color=FAINT, linewidth=0.5)
    ax.set_axisbelow(True)
    handles = [plt.Line2D([], [], color=INK, lw=1.3, label="median")]
    if any_expected:
        # Only advertise a series that is actually drawn.
        handles.append(
            plt.Line2D(
                [], [], color=SERIES[1], lw=1.3, ls=(0, (1.6, 1.2)), label="expected"
            )
        )
    ax.legend(
        handles=handles, loc="lower right", bbox_to_anchor=(1.0, 1.005),
        ncol=len(handles), frameon=False, handlelength=1.6, columnspacing=1.2,
        borderpad=0.0, handletextpad=0.5,
    )
    return _save(figure, path)


def deconvolution_figure(by_plate: Counter[tuple[str, str]], path: Path) -> Path:
    """Read fate per plate: resolved to a source plate, unexpected, or unresolvable."""

    use_print_style()
    plates = sorted({plate for plate, _ in by_plate})
    order = ("resolved", "unexpected_block", "ambiguous", "unknown_block", "unknown_pcr_plate")
    # Status is not identity: only the two informative states take a hue, the
    # rest stay neutral so the eye goes to what needs attention.
    colours = {
        "resolved": SERIES[0],
        "unexpected_block": SERIES[1],
        "ambiguous": SERIES[2],
        "unknown_block": FAINT,
        "unknown_pcr_plate": MUTED,
    }
    figure, ax = plt.subplots(figsize=(SINGLE_COLUMN, 0.24 * len(plates) + 0.80))
    left = [0.0] * len(plates)
    totals = [sum(by_plate[(plate, status)] for status in order) or 1 for plate in plates]
    present = []
    for status in order:
        widths = [
            100 * by_plate[(plate, status)] / total
            for plate, total in zip(plates, totals, strict=True)
        ]
        if not any(widths):
            continue
        present.append(status)
        ax.barh(
            range(len(plates)), widths, left=left, height=0.62,
            color=colours[status], edgecolor="white", linewidth=0.7, zorder=2,
        )
        left = [a + b for a, b in zip(left, widths, strict=True)]
    ax.set_yticks(range(len(plates)))
    ax.set_yticklabels(plates)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("reads (%)")
    ax.xaxis.grid(True, color=FAINT, linewidth=0.5)
    ax.set_axisbelow(True)
    ax.legend(
        handles=[
            Patch(facecolor=colours[s], label=s.replace("_", " ")) for s in present
        ],
        loc="upper center", bbox_to_anchor=(0.5, -0.20),
        ncol=2, frameon=False, handlelength=1.1, handleheight=0.9,
        columnspacing=1.0, handletextpad=0.5, borderpad=0.0,
    )
    return _save(figure, path)


def _save(figure, path: Path, *, tight: bool = True) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    kwargs = {"bbox_inches": "tight"} if tight else {"bbox_inches": None}
    figure.savefig(path.with_suffix(".pdf"), **kwargs)
    figure.savefig(path.with_suffix(".png"), **kwargs)
    plt.close(figure)
    return path.with_suffix(".pdf")


def write_all(
    run_dir: Path,
    output_dir: Path,
    expected_clones: Mapping[str, int] | None = None,
) -> list[Path]:
    """Write every figure for a completed run and return the PDF paths."""

    summaries, _, by_plate = summarize_run(run_dir, expected_clones)
    written = [
        plate_occupancy_figure(summaries, output_dir / "fig1_plate_occupancy"),
        clonality_figure(summaries, output_dir / "fig2_clonality"),
    ]
    if any(status != "not_configured" for _, status in by_plate):
        written.append(
            deconvolution_figure(by_plate, output_dir / "fig3_deconvolution")
        )
    return written


__all__ = [
    "DOUBLE_COLUMN",
    "SERIES",
    "SINGLE_COLUMN",
    "PlateSummary",
    "clonality_figure",
    "deconvolution_figure",
    "plate_occupancy_figure",
    "summarize_run",
    "use_print_style",
    "write_all",
]
