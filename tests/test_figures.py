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


@unittest.skipIf(figures is None, "matplotlib is not installed")
class CulturePlateTests(unittest.TestCase):
    """Deconvolution reporting: how much of each pool a well actually saw."""

    _FIELDS = [
        "consensus_id", "plate_id", "well_id", "reference_library_id",
        "reference_ids", "status", "culture_plate", "n_reads_available",
        "n_reads_used", "mean_depth", "min_depth", "ambiguous_bases",
        "backend", "sequence_sha256", "failure_reason",
    ]

    def setUp(self) -> None:
        self._dir = tempfile.TemporaryDirectory()
        self.run = Path(self._dir.name) / "run"
        stage = self.run / "stages" / "04_consensus"
        stage.mkdir(parents=True)
        rows = []
        # Two wells; the second sees only one of the two pooled plates.
        for well, sources in (("A1", ["CP_A", "CP_B"]), ("A2", ["CP_A"])):
            for i, source in enumerate(sources):
                rows.append({f: "" for f in self._FIELDS} | {
                    "consensus_id": f"c{well}{i}", "plate_id": "RP05",
                    "well_id": well, "reference_ids": f"g{i}",
                    "status": "consensus_pass", "culture_plate": source,
                    "n_reads_available": "10", "n_reads_used": "10",
                    "mean_depth": "10", "ambiguous_bases": "0",
                })
        # A shallow group carries a culture plate but built no sequence.
        rows.append({f: "" for f in self._FIELDS} | {
            "consensus_id": "shallow", "plate_id": "RP05", "well_id": "A2",
            "reference_ids": "g7", "status": "low_depth", "culture_plate": "CP_B",
            "n_reads_available": "2", "n_reads_used": "2", "mean_depth": "2",
            "ambiguous_bases": "0",
        })
        # An undeconvolved plate must not appear at all.
        rows.append({f: "" for f in self._FIELDS} | {
            "consensus_id": "z", "plate_id": "RP01", "well_id": "A1",
            "reference_ids": "g9", "status": "consensus_pass", "culture_plate": "",
            "n_reads_available": "9", "n_reads_used": "9", "mean_depth": "9",
            "ambiguous_bases": "0",
        })
        with gzip.open(stage / "consensus.csv.gz", "wt", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=self._FIELDS, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    def tearDown(self) -> None:
        self._dir.cleanup()

    def test_only_plates_with_deconvolution_are_summarised(self) -> None:
        summaries = figures.summarize_culture_plates(self.run, {"RP05": 2})
        self.assertEqual([s.plate_id for s in summaries], ["RP05"])

    def test_per_well_recovery_counts_distinct_source_plates(self) -> None:
        summary = figures.summarize_culture_plates(self.run, {"RP05": 2})[0]
        self.assertEqual(summary.sources_per_well, {"A1": 2, "A2": 1})
        self.assertEqual(summary.pooled, 2)

    def test_only_built_sequences_are_credited_not_shallow_groups(self) -> None:
        """A group below the depth floor produced no sequence, so it is not recovery."""

        summary = figures.summarize_culture_plates(self.run, {"RP05": 2})[0]
        self.assertEqual(summary.below_depth, 1)
        self.assertEqual(summary.consensus_per_source, {"CP_A": 2, "CP_B": 1})

    def test_consensus_sequences_are_counted_per_source_plate(self) -> None:
        summary = figures.summarize_culture_plates(self.run, {"RP05": 2})[0]
        self.assertEqual(summary.consensus_per_source, {"CP_A": 2, "CP_B": 1})

    def test_pooled_count_falls_back_to_what_was_observed(self) -> None:
        """Without the layout, the figure still draws rather than failing."""

        self.assertEqual(figures.summarize_culture_plates(self.run)[0].pooled, 2)

    def test_the_figure_renders_and_is_skipped_when_there_is_nothing(self) -> None:
        out = Path(self._dir.name) / "figures"
        path = figures.culture_plate_figure(
            figures.summarize_culture_plates(self.run, {"RP05": 2}), out / "fig4"
        )
        self.assertIsNotNone(path)
        self.assertTrue(path.exists() and path.stat().st_size > 0)
        self.assertIsNone(figures.culture_plate_figure([], out / "fig4_empty"))
