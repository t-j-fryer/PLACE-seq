from __future__ import annotations

import unittest

from nanopore3.assignment import ReferenceIndex, assign_sequence, extract_insert
from nanopore3.consensus import ConsensusRead, build_reference_consensus
from nanopore3.demux import call_barcode, validate_barcodes
from nanopore3.qc import evaluate_consensus


class ScienceTests(unittest.TestCase):
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

