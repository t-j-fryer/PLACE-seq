"""Reading-frame and internal-stop QC driven by user-declared constant regions."""

from __future__ import annotations

import unittest

from nanopore3.config import ConfigError, QcSettings
from nanopore3.qc import evaluate_consensus, first_internal_stop, translate

# A miniature construct in the same shape as the 20260506 assay: the start codon
# sits in the upstream constant, the consensus carries a constant linker plus the
# variable insert, and the downstream constant ends with a tag and a stop.
UPSTREAM = "ATG"
LINKER = "CAGCTT"
INSERT = "GTCGTTGGTGGAGTGCTGATTTTCAAA"
DOWNSTREAM = "AGTGGATCCCATCACCACCACCATCACTAA"


class TranslationTests(unittest.TestCase):
    def test_translate_reads_frame_zero_and_marks_stops(self) -> None:
        self.assertEqual(translate("ATGCAGCTTTAA"), "MQL*")

    def test_translate_ignores_a_trailing_partial_codon(self) -> None:
        self.assertEqual(translate("ATGCAGCTTTAAGG"), "MQL*")

    def test_ambiguous_codons_translate_to_x_rather_than_a_guess(self) -> None:
        self.assertEqual(translate("ATGNNNCTT"), "MXL")

    def test_first_internal_stop_ignores_the_terminal_stop(self) -> None:
        self.assertIsNone(first_internal_stop("MQL*"))
        self.assertEqual(first_internal_stop("MQ*L*"), 2)
        self.assertIsNone(first_internal_stop(""))


class CodingQcTests(unittest.TestCase):
    def _evaluate(self, consensus: str, reference: str, **kwargs):
        return evaluate_consensus(
            consensus,
            reference,
            upstream_constant=UPSTREAM,
            downstream_constant=DOWNSTREAM,
            **kwargs,
        )

    def test_coding_checks_are_not_evaluable_without_constant_regions(self) -> None:
        result = evaluate_consensus(LINKER + INSERT, LINKER + INSERT)
        self.assertEqual(result.reading_frame, "not_evaluable")
        self.assertEqual(result.internal_stops, "not_evaluable")
        self.assertEqual(result.protein_length, None)

    def test_an_in_frame_insert_passes_both_coding_checks(self) -> None:
        consensus = LINKER + INSERT
        result = self._evaluate(consensus, consensus)
        self.assertEqual(result.reading_frame, "pass")
        self.assertEqual(result.internal_stops, "pass")
        self.assertEqual(result.overall, "pass")
        # M + Q L + 9 insert residues + S G S H H H H H H, stop excluded.
        expected = len(UPSTREAM + consensus + DOWNSTREAM) // 3 - 1
        self.assertEqual(result.protein_length, expected)
        self.assertIsNone(result.internal_stop_codon)

    def test_a_single_base_insertion_is_reported_as_a_frameshift(self) -> None:
        consensus = LINKER + INSERT + "G"
        result = self._evaluate(consensus, consensus, length_tolerance=10)
        self.assertEqual(result.reading_frame, "fail")
        self.assertEqual(result.overall, "fail")
        self.assertIn("past a codon boundary", result.reason)
        # Stops are not reported from a frameshifted ORF, where they are noise.
        self.assertEqual(result.internal_stops, "not_evaluable")

    def test_a_premature_stop_inside_the_insert_is_caught(self) -> None:
        consensus = LINKER + "GTCGTT" + "TAA" + "GGAGTGCTGATTTTCAAA"
        result = self._evaluate(consensus, consensus)
        self.assertEqual(result.reading_frame, "pass")
        self.assertEqual(result.internal_stops, "fail")
        self.assertEqual(result.overall, "fail")
        # ATG=1, CAG=2, CTT=3, GTC=4, GTT=5, TAA=6 -> index 5, reported as 6.
        self.assertEqual(result.internal_stop_codon, 5)
        self.assertIn("residue 6", result.reason)

    def test_the_terminal_stop_of_the_construct_is_not_an_internal_stop(self) -> None:
        consensus = LINKER + INSERT
        result = self._evaluate(consensus, consensus)
        self.assertEqual(result.internal_stops, "pass")

    def test_coding_failure_is_independent_of_alignment_identity(self) -> None:
        """A consensus can match its reference perfectly and still be frameshifted."""

        consensus = LINKER + INSERT + "G"
        result = self._evaluate(consensus, consensus)
        self.assertEqual(result.full_amplicon, "pass")
        self.assertEqual(result.expected_length, "pass")
        self.assertEqual(result.reading_frame, "fail")
        self.assertEqual(result.overall, "fail")


class CodingConfigTests(unittest.TestCase):
    def test_both_constants_are_required_together(self) -> None:
        with self.assertRaisesRegex(ConfigError, "must be supplied together"):
            QcSettings(upstream_constant="ATG")
        with self.assertRaisesRegex(ConfigError, "must be supplied together"):
            QcSettings(downstream_constant="TAA")

    def test_upstream_constant_must_begin_at_a_start_codon(self) -> None:
        with self.assertRaisesRegex(ConfigError, "start codon"):
            QcSettings(upstream_constant="CAGCTT", downstream_constant="TAA")
        self.assertEqual(
            QcSettings(upstream_constant="ATG", downstream_constant="TAA").upstream_constant,
            "ATG",
        )

    def test_constants_must_be_dna(self) -> None:
        with self.assertRaises(ConfigError):
            QcSettings(upstream_constant="ATGZZZ", downstream_constant="TAA")

    def test_coding_checks_default_to_disabled(self) -> None:
        settings = QcSettings()
        self.assertIsNone(settings.upstream_constant)
        self.assertIsNone(settings.downstream_constant)


if __name__ == "__main__":
    unittest.main()
