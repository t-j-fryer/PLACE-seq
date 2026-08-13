"""Compressed-PCR deconvolution: gene identity as an extra demultiplexing key."""

from __future__ import annotations

import re
import unittest

from nanopore3.config import CompressedPcrSettings, ConfigError
from nanopore3.deconvolution import (
    CompressedPcrPlan,
    DeconvolutionError,
    block_from_reference_id,
)

# RP01 pools two culture plates; RP02 holds the second half of block 3, which is
# split across CP_C and CP_D and therefore must not share one PCR plate.
PCR_PLATES = {"RP01": ("CP_A", "CP_B", "CP_C"), "RP02": ("CP_D",)}
BLOCKS = {"1": ("CP_A",), "2": ("CP_B",), "3": ("CP_C", "CP_D")}


class BlockFromReferenceIdTests(unittest.TestCase):
    def setUp(self) -> None:
        self.pattern = re.compile(r"^Block_(\d+)_")

    def test_the_block_prefix_is_extracted(self) -> None:
        self.assertEqual(
            block_from_reference_id("Block_11_binder_dLK39_458_model", self.pattern),
            "11",
        )

    def test_an_identifier_without_a_block_returns_none(self) -> None:
        self.assertIsNone(block_from_reference_id("dTF141_186_c5", self.pattern))


class PlanValidationTests(unittest.TestCase):
    def test_a_valid_layout_is_accepted(self) -> None:
        plan = CompressedPcrPlan(PCR_PLATES, BLOCKS)
        self.assertEqual(plan.culture_plates, ("CP_A", "CP_B", "CP_C", "CP_D"))
        self.assertEqual(plan.summary()["pcr_plates"], 2)

    def test_a_block_split_within_one_pcr_plate_is_rejected(self) -> None:
        """Both halves of a split block in one PCR plate can never be separated."""

        with self.assertRaises(DeconvolutionError) as caught:
            CompressedPcrPlan({"RP01": ("CP_C", "CP_D")}, {"3": ("CP_C", "CP_D")})
        message = str(caught.exception)
        self.assertIn("unresolvable", message)
        self.assertIn("CP_C", message)
        self.assertIn("different colony PCR plates", message)

    def test_a_culture_plate_never_pooled_anywhere_is_rejected(self) -> None:
        with self.assertRaisesRegex(DeconvolutionError, "never pooled"):
            CompressedPcrPlan({"RP01": ("CP_A",)}, {"1": ("CP_A",), "2": ("CP_ZZ",)})

    def test_empty_layouts_are_rejected(self) -> None:
        with self.assertRaises(DeconvolutionError):
            CompressedPcrPlan({}, BLOCKS)
        with self.assertRaises(DeconvolutionError):
            CompressedPcrPlan(PCR_PLATES, {})


class ResolutionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.plan = CompressedPcrPlan(PCR_PLATES, BLOCKS)

    def test_a_gene_block_names_one_pooled_culture_plate(self) -> None:
        result = self.plan.resolve("RP01", "1")
        self.assertEqual(result.status, "resolved")
        self.assertEqual(result.culture_plate, "CP_A")
        self.assertEqual(result.block, "1")

    def test_two_blocks_in_one_pcr_plate_resolve_to_different_plates(self) -> None:
        """This is the whole point: one PCR well, two source plates, separated."""

        self.assertEqual(self.plan.resolve("RP01", "1").culture_plate, "CP_A")
        self.assertEqual(self.plan.resolve("RP01", "2").culture_plate, "CP_B")

    def test_a_split_block_is_separated_by_the_reverse_barcode(self) -> None:
        # Block 3 spans CP_C and CP_D, which sit in different PCR plates, so the
        # plate barcode decides which half a read came from.
        self.assertEqual(self.plan.resolve("RP01", "3").culture_plate, "CP_C")
        self.assertEqual(self.plan.resolve("RP02", "3").culture_plate, "CP_D")

    def test_a_block_absent_from_this_pcr_plate_is_flagged(self) -> None:
        result = self.plan.resolve("RP02", "1")
        self.assertEqual(result.status, "unexpected_block")
        self.assertIsNone(result.culture_plate)
        self.assertIn("not pooled into", result.reason)

    def test_an_undescribed_pcr_plate_is_flagged(self) -> None:
        self.assertEqual(self.plan.resolve("RP99", "1").status, "unknown_pcr_plate")

    def test_an_unknown_or_missing_block_is_flagged(self) -> None:
        self.assertEqual(self.plan.resolve("RP01", "99").status, "unknown_block")
        self.assertEqual(self.plan.resolve("RP01", None).status, "unknown_block")

    def test_ambiguity_is_reported_rather_than_guessed(self) -> None:
        # Constructed directly, bypassing the layout check, to prove that an
        # ambiguous case is never silently attributed to one plate.
        plan = CompressedPcrPlan.__new__(CompressedPcrPlan)
        plan.pcr_plates = {"RP01": ("CP_C", "CP_D")}
        plan.blocks = {"3": ("CP_C", "CP_D")}
        result = plan.resolve("RP01", "3")
        self.assertEqual(result.status, "ambiguous")
        self.assertIsNone(result.culture_plate)


class CompressedPcrSettingsTests(unittest.TestCase):
    def test_disabled_by_default_and_needs_no_maps(self) -> None:
        settings = CompressedPcrSettings()
        self.assertFalse(settings.enabled)

    def test_enabling_requires_both_maps(self) -> None:
        with self.assertRaisesRegex(ConfigError, "pcr_plates is required"):
            CompressedPcrSettings(enabled=True, blocks=BLOCKS)
        with self.assertRaisesRegex(ConfigError, "blocks is required"):
            CompressedPcrSettings(enabled=True, pcr_plates=PCR_PLATES)

    def test_an_unresolvable_layout_fails_at_configuration_time(self) -> None:
        """Catch a bad pooling design before a run, not as ambiguous reads after."""

        with self.assertRaisesRegex(ConfigError, "unresolvable"):
            CompressedPcrSettings(
                enabled=True,
                pcr_plates={"RP01": ("CP_C", "CP_D")},
                blocks={"3": ("CP_C", "CP_D")},
            )

    def test_a_malformed_block_pattern_is_rejected(self) -> None:
        with self.assertRaisesRegex(ConfigError, "not a regex"):
            CompressedPcrSettings(
                enabled=True, pcr_plates=PCR_PLATES, blocks=BLOCKS, block_pattern="Block_(",
            )


if __name__ == "__main__":
    unittest.main()
