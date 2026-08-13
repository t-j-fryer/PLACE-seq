"""Recover the source culture plate when several were pooled for colony PCR.

In compressed-PCR mode a colony-PCR plate is loaded from more than one culture
plate.  The forward and reverse primer barcodes therefore identify a *PCR* well,
which may hold colonies drawn from several culture plates at the same position.
The assigned gene supplies the missing coordinate: a gene belongs to exactly one
assembly block, and a block was picked into known culture plates, so

``gene -> block -> culture plate``

resolves the well as long as at most one of that block's culture plates was
pooled into this particular PCR plate.  When a block is split across culture
plates, those plates must go to different PCR plates so the reverse barcode
separates them; otherwise the read is genuinely ambiguous and is reported as
such rather than being attributed to an arbitrary plate.
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Literal

DeconvolutionStatus = Literal[
    "resolved",
    "ambiguous",
    "unexpected_block",
    "unknown_block",
    "unknown_pcr_plate",
    "not_configured",
]


class DeconvolutionError(ValueError):
    """A compressed-PCR layout cannot resolve the wells it describes."""


@dataclass(frozen=True, slots=True)
class Deconvolution:
    """Where one read's colony physically came from."""

    culture_plate: str | None
    block: str | None
    status: DeconvolutionStatus
    reason: str


def block_from_reference_id(reference_id: str, pattern: re.Pattern[str]) -> str | None:
    """Extract an assembly block from a reference identifier, if it encodes one."""

    match = pattern.search(reference_id)
    if match is None:
        return None
    return match.group(1) if match.groups() else match.group(0)


class CompressedPcrPlan:
    """Validated pooling layout: which culture plates each PCR plate was loaded from."""

    def __init__(
        self,
        pcr_plates: Mapping[str, tuple[str, ...]],
        blocks: Mapping[str, tuple[str, ...]],
    ) -> None:
        if not pcr_plates:
            raise DeconvolutionError("compressed_pcr.pcr_plates must not be empty")
        if not blocks:
            raise DeconvolutionError("compressed_pcr.blocks must not be empty")
        self.pcr_plates = {
            plate: tuple(sorted(set(sources))) for plate, sources in pcr_plates.items()
        }
        self.blocks = {
            block: tuple(sorted(set(plates))) for block, plates in blocks.items()
        }

        loaded = {source for sources in self.pcr_plates.values() for source in sources}
        orphaned = sorted(
            {plate for plates in self.blocks.values() for plate in plates} - loaded
        )
        if orphaned:
            raise DeconvolutionError(
                "culture plate(s) never pooled into any colony PCR plate, so reads "
                f"from them can never be resolved: {', '.join(orphaned)}"
            )

        # A block split across two culture plates that share one PCR plate is
        # unresolvable by construction. That is a layout error worth catching
        # before a run rather than discovering as ambiguous reads afterwards.
        collisions = []
        for block, plates in sorted(self.blocks.items()):
            for pcr_plate, sources in sorted(self.pcr_plates.items()):
                shared = sorted(set(plates) & set(sources))
                if len(shared) > 1:
                    collisions.append(
                        f"block {block} occupies {', '.join(shared)} which were both "
                        f"pooled into {pcr_plate}"
                    )
        if collisions:
            raise DeconvolutionError(
                "compressed PCR layout is unresolvable: "
                + "; ".join(collisions)
                + ". Split these culture plates across different colony PCR plates "
                "so the reverse barcode separates them."
            )

    def resolve(self, pcr_plate: str, block: str | None) -> Deconvolution:
        """Resolve one read's source culture plate from its PCR plate and block."""

        if block is None:
            return Deconvolution(
                None, None, "unknown_block",
                "the assigned reference does not identify an assembly block",
            )
        sources = self.pcr_plates.get(pcr_plate)
        if sources is None:
            return Deconvolution(
                None, block, "unknown_pcr_plate",
                f"plate barcode {pcr_plate!r} is not described by compressed_pcr.pcr_plates",
            )
        plates = self.blocks.get(block)
        if plates is None:
            return Deconvolution(
                None, block, "unknown_block",
                f"block {block!r} is not described by compressed_pcr.blocks",
            )
        candidates = sorted(set(plates) & set(sources))
        if len(candidates) == 1:
            return Deconvolution(
                candidates[0], block, "resolved",
                "gene block identifies exactly one pooled culture plate",
            )
        if not candidates:
            # The block was never loaded into this PCR plate, so the colony is
            # misplaced, cross-contaminating, or the layout is misdescribed.
            return Deconvolution(
                None, block, "unexpected_block",
                f"block {block!r} was not pooled into {pcr_plate!r}",
            )
        return Deconvolution(
            None, block, "ambiguous",
            f"block {block!r} maps to several culture plates pooled into "
            f"{pcr_plate!r}: {', '.join(candidates)}",
        )

    @property
    def culture_plates(self) -> tuple[str, ...]:
        return tuple(sorted({p for plates in self.blocks.values() for p in plates}))

    def summary(self) -> dict[str, int]:
        return {
            "pcr_plates": len(self.pcr_plates),
            "culture_plates": len(self.culture_plates),
            "blocks": len(self.blocks),
        }


NOT_CONFIGURED = Deconvolution(
    None, None, "not_configured", "compressed PCR deconvolution is not enabled"
)


__all__ = [
    "NOT_CONFIGURED",
    "CompressedPcrPlan",
    "Deconvolution",
    "DeconvolutionError",
    "DeconvolutionStatus",
    "block_from_reference_id",
]
