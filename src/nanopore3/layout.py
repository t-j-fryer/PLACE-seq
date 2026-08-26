"""Read a pooling layout from a spreadsheet instead of nested YAML.

The compressed-PCR layout is *experimental design*: which assembly block was
picked into which culture plate, and which culture plates were pooled into which
colony-PCR barcode.  It is known before any read exists and must never be inferred
from the data - inferring it would use the reads to build the key that then
interprets those same reads.

It is also the most tedious part of a configuration to write by hand: the 260608
run needs sixty lines of nested YAML for what a lab records as one table.  So it
can be supplied as a CSV, in the shape the information is actually held:

    plate_barcode,culture_plate,library,block
    RP05,SUMO_A_P1,sumo_ab,A_Block_1
    RP05,SUMO_A_P6,sumo_ab,A_Block_6
    RP05,SUMO_A_P6,sumo_ab,A_Block_8
    RP08,LAB_P1,lab,Block_1

One row per (culture plate, block) pair - "this block went on this plate, and this
plate went into this barcode".  A block picked across several plates gets several
rows; so does a plate holding several blocks.  Both are normal.

``clonality`` is an optional fifth column, and applies to the barcode: give it on
any row and it must agree on every row for that barcode.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

REQUIRED_COLUMNS = ("plate_barcode", "culture_plate", "library", "block")
OPTIONAL_COLUMNS = ("clonality",)


class LayoutError(ValueError):
    """The layout table cannot be read, or contradicts itself."""


def read_layout_csv(path: Path) -> dict[str, object]:
    """Return ``pcr_plates``, ``blocks`` and ``clonality`` from a layout table.

    Ordering is preserved as first-seen rather than sorted, so a table written in
    plate order produces a configuration in plate order and a diff between two
    layouts stays readable.
    """

    if not path.is_file():
        raise LayoutError(f"compressed_pcr.layout_csv not found: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = [name.strip() for name in (reader.fieldnames or [])]
        missing = [name for name in REQUIRED_COLUMNS if name not in fields]
        if missing:
            raise LayoutError(
                f"{path}: missing column(s) {', '.join(missing)}; "
                f"expected {', '.join(REQUIRED_COLUMNS)}"
                + (f" and optionally {', '.join(OPTIONAL_COLUMNS)}" if fields else "")
            )
        unknown = [
            name for name in fields if name not in REQUIRED_COLUMNS + OPTIONAL_COLUMNS
        ]
        if unknown:
            raise LayoutError(
                f"{path}: unknown column(s) {', '.join(unknown)}. A typo in a column "
                "name would otherwise be read as a missing value on every row."
            )

        plates: dict[str, list[str]] = defaultdict(list)
        blocks: dict[str, dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))
        clonality: dict[str, str] = {}
        rows = 0
        for number, row in enumerate(reader, start=2):
            values = {name: (row.get(name) or "").strip() for name in fields}
            if not any(values.values()):  # a blank line, which spreadsheets add
                continue
            rows += 1
            for name in REQUIRED_COLUMNS:
                if not values[name]:
                    raise LayoutError(f"{path} line {number}: {name} is empty")
            barcode, culture = values["plate_barcode"], values["culture_plate"]
            library, block = values["library"], values["block"]
            if culture not in plates[barcode]:
                plates[barcode].append(culture)
            if culture not in blocks[library][block]:
                blocks[library][block].append(culture)
            mode = values.get("clonality", "")
            if mode:
                previous = clonality.setdefault(barcode, mode)
                if previous != mode:
                    raise LayoutError(
                        f"{path} line {number}: clonality for {barcode} is {mode!r} here "
                        f"but {previous!r} on an earlier row; it describes the barcode, "
                        "so it must be the same on every row for it"
                    )
    if not rows:
        raise LayoutError(f"{path}: no rows")
    return {
        "pcr_plates": {k: tuple(v) for k, v in plates.items()},
        "blocks": {
            library: {block: tuple(p) for block, p in by_block.items()}
            for library, by_block in blocks.items()
        },
        "clonality": clonality,
    }


def write_layout_csv(
    path: Path,
    pcr_plates: dict[str, tuple[str, ...]],
    blocks: dict[str, dict[str, tuple[str, ...]]],
    clonality: dict[str, str] | None = None,
) -> Path:
    """Write an existing YAML layout back out as a table.

    For migrating a configuration that already has the nested form, so nobody has
    to retype sixty lines to adopt the CSV.
    """

    plate_barcode = {
        culture: barcode
        for barcode, cultures in pcr_plates.items()
        for culture in cultures
    }
    rows = []
    for library, by_block in blocks.items():
        for block, cultures in by_block.items():
            for culture in cultures:
                barcode = plate_barcode.get(culture, "")
                row = {
                    "plate_barcode": barcode,
                    "culture_plate": culture,
                    "library": library,
                    "block": block,
                }
                if clonality:
                    row["clonality"] = clonality.get(barcode, "")
                rows.append(row)
    fields = list(REQUIRED_COLUMNS) + (["clonality"] if clonality else [])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return path


__all__ = ["LayoutError", "read_layout_csv", "write_layout_csv", "REQUIRED_COLUMNS"]
