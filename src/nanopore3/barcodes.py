"""Validated barcode registries shared by assay profiles.

The registry is deliberately a plain CSV so it can be reviewed in Excel, Git,
or a text editor.  Scientific runs select an explicit family; identifiers from
different primer contexts are never pooled merely because they share an index.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from .errors import Nanopore3Error
from .provenance import sha256_file


class BarcodeRegistryError(Nanopore3Error, ValueError):
    """A malformed registry or unsafe family selection."""


@dataclass(frozen=True, slots=True)
class BarcodePanel:
    family_id: str
    sequences: dict[str, str]
    source_path: Path
    source_sha256: str


def read_barcode_panel(path: str | Path, family_id: str) -> BarcodePanel:
    """Load exactly one barcode family from the canonical registry schema."""

    source = Path(path).expanduser().resolve(strict=True)
    if not source.is_file():
        raise BarcodeRegistryError(f"barcode registry is not a file: {source}")
    family = family_id.strip()
    if not family:
        raise BarcodeRegistryError("barcode family_id must be non-empty")
    required = {"family_id", "barcode_id", "sequence"}
    sequences: dict[str, str] = {}
    with source.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise BarcodeRegistryError(
                f"{source}: registry requires columns: {', '.join(sorted(required))}"
            )
        for line_number, row in enumerate(reader, start=2):
            if (row.get("family_id") or "").strip() != family:
                continue
            barcode_id = (row.get("barcode_id") or "").strip()
            sequence = (row.get("sequence") or "").strip().upper()
            if not barcode_id or not sequence:
                raise BarcodeRegistryError(
                    f"{source}:{line_number}: empty barcode_id or sequence"
                )
            if barcode_id in sequences:
                raise BarcodeRegistryError(
                    f"{source}:{line_number}: duplicate barcode_id {barcode_id!r} "
                    f"within family {family!r}"
                )
            sequences[barcode_id] = sequence
    if not sequences:
        raise BarcodeRegistryError(
            f"barcode family {family!r} was not found in {source}"
        )
    return BarcodePanel(family, sequences, source, sha256_file(source))


__all__ = ["BarcodePanel", "BarcodeRegistryError", "read_barcode_panel"]
