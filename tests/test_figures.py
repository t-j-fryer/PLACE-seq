"""Figures must render from a run's tables without a display or network."""

from __future__ import annotations

import csv
import gzip
import tempfile
import unittest
from pathlib import Path

try:
    from nanopore3 import figures
except ImportError:  # pragma: no cover - exercised only without the report extra
    figures = None

_FIELDS = [
    "read_uid", "plate_id", "well_id", "reference_ids",
    "assignment_status", "culture_plate", "culture_plate_status",
]


def _write_run(root: Path) -> Path:
    stage = root / "stages" / "03_assignment"
    stage.mkdir(parents=True)
    rows = []
    for well_index, well in enumerate(["A1", "A2", "B1", "H12"]):
        for gene in range(well_index + 1):
            rows.append(
                {
                    "read_uid": f"{well}-{gene}", "plate_id": "RP01", "well_id": well,
                    "reference_ids": f"Block_1_gene{gene}",
                    "assignment_status": "assigned_unique",
                    "culture_plate": "CP_A", "culture_plate_status": "resolved",
                }
            )
    rows.append(
        {
            "read_uid": "x", "plate_id": "RP01", "well_id": "A1",
            "reference_ids": "", "assignment_status": "no_match",
            "culture_plate": "", "culture_plate_status": "unknown_block",
        }
    )
    with gzip.open(stage / "assignment_calls.csv.gz", "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return root


@unittest.skipIf(figures is None, "matplotlib is not installed")
class FigureTests(unittest.TestCase):
    def setUp(self) -> None:
        self._dir = tempfile.TemporaryDirectory()
        self.run = _write_run(Path(self._dir.name) / "run")

    def tearDown(self) -> None:
        self._dir.cleanup()

    def test_summary_counts_distinct_genes_not_reads(self) -> None:
        summaries, status, by_plate = figures.summarize_run(self.run)
        self.assertEqual(len(summaries), 1)
        summary = summaries[0]
        # A2 has two reads of two genes; H12 has four. Reads are not clones.
        self.assertEqual(summary.clones_per_well["A1"], 1)
        self.assertEqual(summary.clones_per_well["H12"], 4)
        self.assertEqual(status["resolved"], 10)
        self.assertEqual(by_plate[("RP01", "unknown_block")], 1)

    def test_unassigned_reads_do_not_count_as_clones(self) -> None:
        summaries, _, _ = figures.summarize_run(self.run)
        self.assertEqual(summaries[0].clones_per_well["A1"], 1)

    def test_expected_clonality_is_attached_when_supplied(self) -> None:
        summaries, _, _ = figures.summarize_run(self.run, {"RP01": 4})
        self.assertEqual(summaries[0].expected_clones, 4)

    def test_every_figure_renders_as_pdf_and_png(self) -> None:
        out = Path(self._dir.name) / "figures"
        written = figures.write_all(self.run, out, {"RP01": 4})
        self.assertTrue(written)
        for path in written:
            self.assertTrue(path.exists() and path.stat().st_size > 0)
            self.assertTrue(path.with_suffix(".png").exists())

    def test_the_palette_is_the_validated_one(self) -> None:
        """Guards against a hue being edited in without re-running the check."""

        self.assertEqual(figures.SERIES, ("#0072B2", "#D55E00", "#009E73"))


if __name__ == "__main__":
    unittest.main()
