"""Streaming FASTQ exports with bounded open handles and deterministic gzip."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
from collections import OrderedDict
from pathlib import Path

from .export import safe_name


def _component(value: str) -> str:
    cleaned = safe_name(value, limit=40)
    reserved = {"CON", "PRN", "AUX", "NUL"} | {
        f"{prefix}{number}" for prefix in ("COM", "LPT") for number in range(1, 10)
    }
    if cleaned != value or value != value.upper() or cleaned.split(".")[0] in reserved:
        cleaned += "-" + hashlib.sha256(value.encode()).hexdigest()[:12]
    return cleaned


class DemuxFastqWriter:
    """Write plate/well views. Reopened files contain standard gzip members.

    Each eligible read appears in its plate file and either its well file or
    the unresolved-well file. These are overlapping views, not disjoint inputs.
    Only the parent process writes; ordering follows the original input stream.
    """

    def __init__(self, root: Path, max_open: int = 32, *, include_wells: bool = True):
        if max_open < 1:
            raise ValueError("max_open must be positive")
        self.root = root
        self.max_open = max_open
        self.include_wells = include_wells
        self.handles = OrderedDict()
        self.entries = {}

    def __enter__(self):
        self.root.mkdir(parents=True, exist_ok=False)
        return self

    def _close(self, wrappers):
        text, compressed, raw = wrappers
        try:
            text.close()
        finally:
            compressed.close()
            raw.close()

    def _append(self, relative, payload, *, scope, plate, well):
        if relative in self.handles:
            wrappers = self.handles.pop(relative)
        else:
            if len(self.handles) >= self.max_open:
                self._close(self.handles.popitem(last=False)[1])
            path = self.root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            raw = path.open("ab")
            compressed = gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0)
            wrappers = (io.TextIOWrapper(compressed, encoding="utf-8", newline=""), compressed, raw)
        self.handles[relative] = wrappers
        wrappers[0].write(payload)
        if relative not in self.entries:
            self.entries[relative] = dict(
                scope=scope, plate_id=plate, well_id=well, reads=0, file=relative
            )
        entry = self.entries[relative]
        if (entry["scope"], entry["plate_id"], entry["well_id"]) != (scope, plate, well):
            raise ValueError("FASTQ export filename collision")
        entry["reads"] += 1

    def write(self, read, *, assigned: bool):
        plate = str(read["plate_id"])
        well = str(read["well_id"]) if assigned else ""
        p = _component(plate)
        payload = (
            f"@{read['original_read_id']} read_uid={read['read_uid']} "
            f"sample_id={read['sample_id']}\n{read['sequence']}\n+\n{read['quality']}\n"
        )
        self._append(f"by_plate/{p}.fastq.gz", payload, scope="plate", plate=plate, well="")
        if not self.include_wells:
            return
        if assigned:
            relative, scope = f"by_well/{p}/{_component(well)}.fastq.gz", "well"
        else:
            relative, scope = f"unresolved_well/{p}.fastq.gz", "unresolved_well"
        self._append(relative, payload, scope=scope, plate=plate, well=well)

    def __exit__(self, exc_type, exc, tb):
        while self.handles:
            self._close(self.handles.popitem()[1])
        if exc_type is None:
            with (self.root / "index.csv").open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=("scope", "plate_id", "well_id", "reads", "file"),
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(self.entries[k] for k in sorted(self.entries))
