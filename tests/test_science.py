from __future__ import annotations

import unittest

from nanopore3.assignment import ReferenceIndex, assign_sequence, extract_insert
from nanopore3.consensus import ConsensusRead, build_reference_consensus
from nanopore3.config import ParallelSettings
from nanopore3.demux import call_barcode, prepare_barcode_panel, validate_barcodes
from nanopore3.qc import evaluate_consensus


class ScienceTests(unittest.TestCase):
    def test_process_backend_is_explicit_and_validated(self) -> None:
        self.assertEqual(ParallelSettings(backend="process").backend, "process")
        with self.assertRaisesRegex(ValueError, "parallel.backend"):
            ParallelSettings(backend="fork")

    def test_barcode_call_retains_margin_evidence(self) -> None:
        barcodes = validate_barcodes(
            {"P1": "AAAACCCC", "P2": "CCAATTGG"}, max_edits=1, min_margin=2
        )
        call = call_barcode(
            "AAAACCCCGATTACA",
            barcodes,
            window_size=10,
            max_edits=1,
            min_margin=2,
            search_ends=("head",),
        )
        self.assertEqual(call.status, "assigned")
        self.assertEqual(call.barcode_id, "P1")
        self.assertEqual(call.best_distance, 0)

    def test_legacy_unique_threshold_matches_optimization_sweep_policy(self) -> None:
        panel = {"A": "AAAACCCC", "B": "AAAAGGGG"}
        one_hit = call_barcode(
            "TTTAAAACCCCTTT",
            panel,
            window_size=20,
            max_edits=1,
            search_ends=("head",),
            allow_reverse_complement=False,
            decision_policy="legacy_unique_threshold",
        )
        self.assertEqual((one_hit.status, one_hit.barcode_id), ("assigned", "A"))
        multi_hit = call_barcode(
            "TTTAAAACCGGTTT",
            panel,
            window_size=20,
            max_edits=2,
            search_ends=("head",),
            allow_reverse_complement=False,
            decision_policy="legacy_unique_threshold",
        )
        self.assertEqual(multi_hit.status, "ambiguous")
        self.assertIsNone(multi_hit.barcode_id)

    def test_legacy_unique_best_matches_active_well_policy(self) -> None:
        panel = {"A": "AAAACCCC", "B": "AAAAGGGG"}
        call = call_barcode(
            "TTTAAAACCCCTTT",
            panel,
            window_size=20,
            max_edits=4,
            search_ends=("head",),
            allow_reverse_complement=False,
            decision_policy="legacy_unique_best",
        )
        self.assertEqual((call.status, call.barcode_id), ("assigned", "A"))

    def test_prepared_panel_is_equivalent_and_legacy_ties_prefer_tail(self) -> None:
        panel = {"A": "AAAACCCC"}
        read = "AAAACCCCNNNNAAAACCCC"
        expected = call_barcode(
            read,
            panel,
            window_size=8,
            max_edits=0,
            decision_policy="legacy_unique_threshold",
        )
        observed = call_barcode(
            read,
            panel,
            window_size=8,
            max_edits=0,
            decision_policy="legacy_unique_threshold",
            prepared_panel=prepare_barcode_panel(panel),
        )
        self.assertEqual(observed, expected)
        self.assertEqual(observed.matched_end, "tail")
        self.assertEqual(observed.orientation, "forward")

    def test_motif_extraction_never_silently_uses_full_read(self) -> None:
        found = extract_insert("AAAAGATTACAACGTCCGGAATTT", "GATTACA", "CCGGAAT")
        self.assertEqual(found.status, "found")
        self.assertEqual(found.sequence, "ACGT")
        missing = extract_insert("AAAACCCC", "GATTACA", "CCGGAAT")
        self.assertEqual(missing.status, "motif_missing")

    def test_assignment_reports_alias_set_and_rejects_missing_motif(self) -> None:
        references = {
            "a": "ACGTACGTACGT",
            "alias": "ACGTACGTACGT",
            "b": "TTTTGGGGCCCC",
        }
        index = ReferenceIndex(references, k=5, max_kmer_owners=2)
        call = assign_sequence(
            "GATTACAACGTACGTACGTCCGGAAT",
            index,
            left_motif="GATTACA",
            right_motif="CCGGAAT",
            motif_max_edits=0,
            min_kmer_score=1,
            min_identity=0.95,
            min_query_coverage=0.95,
            min_reference_coverage=0.95,
        )
        self.assertEqual(call.status, "assigned_alias_set")
        self.assertEqual(call.reference_ids, ("a", "alias"))
        missing = assign_sequence(
            "ACGTACGTACGT",
            index,
            left_motif="GATTACA",
            right_motif="CCGGAAT",
        )
        self.assertEqual(missing.status, "motif_missing")

    def test_consensus_is_input_order_invariant_and_qc_is_tristate(self) -> None:
        reads = [ConsensusRead(f"r{i}", "ACGTACGT") for i in range(6)]
        first = build_reference_consensus("ACGTACGT", reads, group_id="g", max_reads=4)
        second = build_reference_consensus(
            "ACGTACGT", list(reversed(reads)), group_id="g", max_reads=4
        )
        self.assertEqual(first, second)
        qc = evaluate_consensus(first.sequence, "ACGTACGT")
        self.assertEqual(qc.overall, "pass")
        self.assertEqual(qc.reading_frame, "not_evaluable")


if __name__ == "__main__":
    unittest.main()
