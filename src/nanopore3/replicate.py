"""Concordance between a dedicated plate and the same culture plate in a pool.

A compressed colony PCR pools several culture plates into one barcoded plate, so
one well of that barcode holds one clone per pooled culture plate.  Where a
culture plate was *also* sequenced on a barcode of its own, the two datasets are
independent measurements of the same physical wells: the dedicated barcode reads
one clone per well directly, the pooled barcode reads the same clone alongside
many others and relies on the block-to-culture-plate map to separate them.

Agreement between the two therefore measures the whole chain end to end -
demultiplexing, assignment, consensus and deconvolution - against data rather
than against itself.  Disagreement is informative in both directions: a well
recovered only from the pool is extra yield, and a well recovered only from the
dedicated barcode is a deconvolution miss.

Chimeric clones participate.  A chimera is real sequence in a real well, and its
parent signature identifies it as well as a design name does, so a chimera seen
in both datasets is a match like any other.
"""

from __future__ import annotations

import csv
import gzip
import json
import re
from dataclasses import dataclass
from pathlib import Path

# A call must carry a sequence to be comparable.  `mixed_variants` does - it is a
# polyclonal well whose consensus is reported with ambiguity codes - so it counts.
BUILT_STATUSES = frozenset({"consensus_pass", "mixed_variants"})

MATCHED_EXACT = "matched_exact"
MATCHED_DIFFERENT_SEQUENCE = "matched_different_sequence"
DIFFERENT_DESIGN = "different_design"
POOL_ONLY = "pool_only"
DEDICATED_ONLY = "dedicated_only"

# Presentation order for stacked bars: agreement first, then disagreement.
CATEGORY_ORDER = (
    MATCHED_EXACT,
    MATCHED_DIFFERENT_SEQUENCE,
    DIFFERENT_DESIGN,
    POOL_ONLY,
    DEDICATED_ONLY,
)

CATEGORY_LABELS = {
    MATCHED_EXACT: "Matched, exact seq.",
    MATCHED_DIFFERENT_SEQUENCE: "Matched, different seq.",
    DIFFERENT_DESIGN: "Different clone",
    POOL_ONLY: "Pooled barcode only",
    DEDICATED_ONLY: "Dedicated barcode only",
}


@dataclass(frozen=True)
class Call:
    """One sequenced clone in one well."""

    identity: str
    kind: str  # "consensus" or "chimera"
    status: str
    reads: int
    sequence: str

    @property
    def sequence_available(self) -> bool:
        return bool(self.sequence)


@dataclass(frozen=True)
class WellOutcome:
    """The comparison for a single well."""

    well_id: str
    category: str
    dedicated: tuple[str, ...]
    pooled: tuple[str, ...]
    shared: tuple[str, ...]
    exact: tuple[str, ...]
    differing: tuple[str, ...]
    edit_distances: tuple[tuple[str, int], ...]


def _normalise_well(well_id: str) -> str:
    """`A1` and `A01` are the same well; compare on a single form."""

    match = re.fullmatch(r"([A-Ha-h])0*(\d{1,2})", well_id.strip())
    if not match:
        return well_id.strip().upper()
    return f"{match.group(1).upper()}{int(match.group(2)):02d}"


def read_fasta(path: Path) -> dict[str, str]:
    """Read a FASTA keyed by the first whitespace-delimited header token."""

    sequences: dict[str, str] = {}
    if not path.exists():
        return sequences
    name: str | None = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(chunks)
                name = line[1:].split()[0] if len(line) > 1 else ""
                chunks = []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        sequences[name] = "".join(chunks)
    return sequences


def load_calls(
    run_dir: Path,
    plate_id: str,
    *,
    culture_plate: str | None = None,
    include_chimeras: bool = True,
) -> dict[str, dict[str, Call]]:
    """Collect every sequenced clone for one barcode, keyed by well then identity.

    `culture_plate` restricts a pooled barcode to the clones the deconvolution
    assigned to that source plate; leave it None for a dedicated barcode.
    """

    stages = run_dir / "stages"
    sequences = read_fasta(stages / "04_consensus" / "consensus.fasta")
    sequences.update(read_fasta(stages / "04b_chimera" / "clones" / "scaffolds.fasta"))

    wells: dict[str, dict[str, Call]] = {}

    consensus_csv = stages / "04_consensus" / "consensus.csv.gz"
    with gzip.open(consensus_csv, "rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["plate_id"] != plate_id or row["status"] not in BUILT_STATUSES:
                continue
            if culture_plate is not None and row["culture_plate"] != culture_plate:
                continue
            well = _normalise_well(row["well_id"])
            wells.setdefault(well, {})[row["reference_ids"]] = Call(
                identity=row["reference_ids"],
                kind="consensus",
                status=row["status"],
                reads=int(row["n_reads_used"] or 0),
                sequence=sequences.get(row["consensus_id"], ""),
            )

    clones_csv = stages / "04b_chimera" / "clones" / "clones.csv"
    if include_chimeras and clones_csv.exists():
        with clones_csv.open("r", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                if row["plate_id"] != plate_id:
                    continue
                if culture_plate is not None and row["culture_plate"] != culture_plate:
                    continue
                well = _normalise_well(row["well_id"])
                wells.setdefault(well, {})[row["parents"]] = Call(
                    identity=row["parents"],
                    kind="chimera",
                    status="chimera",
                    reads=int(row["reads"] or 0),
                    sequence=sequences.get(row["chimera_id"], ""),
                )

    return wells


def _edit_distance(left: str, right: str) -> int:
    """Edit distance, for reporting how far apart two non-identical calls are."""

    import edlib

    return int(edlib.align(left, right, mode="NW", task="distance")["editDistance"])


def compare_wells(
    dedicated: dict[str, dict[str, Call]],
    pooled: dict[str, dict[str, Call]],
) -> list[WellOutcome]:
    """Classify every well seen in either dataset.

    A well matches if the two datasets name at least one clone in common.  It is
    exact only if *every* shared clone has a byte-identical sequence, so a
    polyclonal well cannot be called exact on the strength of one of its clones.
    """

    outcomes: list[WellOutcome] = []
    for well in sorted(set(dedicated) | set(pooled)):
        left = dedicated.get(well, {})
        right = pooled.get(well, {})
        shared = sorted(set(left) & set(right))

        exact: list[str] = []
        differing: list[str] = []
        distances: list[tuple[str, int]] = []
        for identity in shared:
            a, b = left[identity].sequence, right[identity].sequence
            if a and b and a == b:
                exact.append(identity)
            else:
                differing.append(identity)
                if a and b:
                    distances.append((identity, _edit_distance(a, b)))

        if shared:
            category = MATCHED_EXACT if not differing else MATCHED_DIFFERENT_SEQUENCE
        elif left and right:
            category = DIFFERENT_DESIGN
        elif left:
            category = DEDICATED_ONLY
        else:
            category = POOL_ONLY

        outcomes.append(
            WellOutcome(
                well_id=well,
                category=category,
                dedicated=tuple(sorted(left)),
                pooled=tuple(sorted(right)),
                shared=tuple(shared),
                exact=tuple(exact),
                differing=tuple(differing),
                edit_distances=tuple(distances),
            )
        )
    return outcomes


def summarise(outcomes: list[WellOutcome]) -> dict[str, object]:
    """Well-level and clone-level counts for one plate pair."""

    counts = {name: 0 for name in CATEGORY_ORDER}
    for outcome in outcomes:
        counts[outcome.category] += 1

    dedicated_wells = sum(1 for o in outcomes if o.dedicated)
    matched = counts[MATCHED_EXACT] + counts[MATCHED_DIFFERENT_SEQUENCE]
    shared_clones = sum(len(o.shared) for o in outcomes)
    exact_clones = sum(len(o.exact) for o in outcomes)

    return {
        "wells_total": len(outcomes),
        "wells_dedicated": dedicated_wells,
        "wells_pooled": sum(1 for o in outcomes if o.pooled),
        "counts": counts,
        "wells_matched": matched,
        # Concordance is measured over wells the dedicated barcode actually saw:
        # a well only the pool recovered is extra yield, not a disagreement.
        "matched_fraction": matched / dedicated_wells if dedicated_wells else 0.0,
        "clones_shared": shared_clones,
        "clones_exact": exact_clones,
        "exact_fraction": exact_clones / shared_clones if shared_clones else 0.0,
        "chimeras_matched": sum(
            1
            for o in outcomes
            for identity in o.shared
            if " >> " in identity
        ),
    }


def verify_pairing(run_dir: Path, plate_id: str, culture_plate: str) -> dict[str, object]:
    """Check from the data that a dedicated barcode really is that culture plate.

    The pairing is an experimental claim, so it is tested rather than trusted:
    every design on the dedicated barcode carries a block name, and the run's own
    block map says which culture plate that block came from.
    """

    config = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))["config"]
    compressed = config.get("compressed_pcr") or {}
    pattern = re.compile(compressed.get("block_pattern") or r"^((?:[AB]_)?Block_\d+)_")
    blocks: dict[str, dict[str, list[str]]] = compressed.get("blocks") or {}

    calls = load_calls(run_dir, plate_id, include_chimeras=False)
    resolved: dict[str, int] = {}
    unresolved = 0
    for well_calls in calls.values():
        for identity in well_calls:
            match = pattern.match(identity)
            if not match:
                unresolved += 1
                continue
            block = match.group(1)
            sources = {
                plate
                for library_blocks in blocks.values()
                for name, plate_list in library_blocks.items()
                if name == block
                for plate in plate_list
            }
            if not sources:
                unresolved += 1
                continue
            for plate in sources:
                resolved[plate] = resolved.get(plate, 0) + 1

    total = sum(resolved.values())
    consistent = resolved.get(culture_plate, 0)
    return {
        "plate_id": plate_id,
        "expected_culture_plate": culture_plate,
        "clones_examined": total + unresolved,
        "clones_consistent": consistent,
        "clones_unresolved": unresolved,
        "sources_seen": dict(sorted(resolved.items())),
        "confirmed": total > 0 and consistent == total and unresolved == 0,
    }
