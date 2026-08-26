"""Compare clone-picked nanopore sequencing against pooled Illumina sequencing.

The two platforms answer the same question by different means.  Nanopore reads
one picked colony per well, so a design is recovered when some well held it.
Illumina reads the pool, so a design is recovered when some read carried it.
Comparing recovery between them is only fair at matched sampling effort: picking
N colonies from a block is the same investment as taking N reads of that block,
so the Illumina data is subsampled per block to the number of wells the nanopore
run actually sequenced of that block.

Reference names agree exactly across the two datasets once the block prefix is
parsed - the Illumina names separate the encoding with ``|`` where ours use
``_`` - so designs join on (encoding, block, design key) with nothing left over.

Three quantities are computed here:

* distinct sequences per allocated well - how often a well held one clone,
  none, or several;
* recovered references - what fraction of a designed library was seen at all,
  split by how good the best observation was;
* read accuracy - per-read identity, which is where the platforms differ most
  and where the nanopore consensus earns its keep.
"""

from __future__ import annotations

import csv
import gzip
import random
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

# ``A|Block_9_x`` (Illumina) and ``A_Block_9_x`` (ours) name the same design.
REFERENCE_NAME = re.compile(r"^(?:([AB])[|_])?Block_(\d+)_")

# A well never yields a sequence without a consensus, so ``low_depth`` is absence
# of data rather than a poor result and is excluded from every count here.
NO_SEQUENCE_GRADE = "low_depth"

PERFECT, SCREENABLE, OTHER = "Perfect", "Screenable", "Other"
CLASS_ORDER = (PERFECT, SCREENABLE, OTHER)
_CLASS_RANK = {name: index for index, name in enumerate(CLASS_ORDER)}

# Our consensus grades and the Illumina read categories describe the same three
# outcomes: right, usable, and not usable.  ``screenable`` and
# ``synonymous``/``nonsynonymous`` both mean full length and in frame with
# substitutions, which is a clone you would still put in an assay.
GRADE_CLASS = {
    "perfect": PERFECT,
    "screenable": SCREENABLE,
    # A damage-signature mixture is a real clone carrying a real, heritable base -
    # it screens like any other non-exact clone. An unexplained mixture does not.
    "mixed_damage": SCREENABLE,
    "mixed_variants": OTHER,
    "mismatched": OTHER,
    "truncated": OTHER,
    "premature_stop": OTHER,
    "frameshift": OTHER,
}
CATEGORY_CLASS = {
    "perfect_match": PERFECT,
    "synonymous": SCREENABLE,
    "nonsynonymous": SCREENABLE,
    "premature_stop": OTHER,
    "frameshift": OTHER,
}


@dataclass(frozen=True)
class LibrarySet:
    """One comparable library: where its clones are, and what its designs are.

    ``encodings`` holds more than one entry for a set that unions encodings of
    the same designs, where a design counts as recovered if either encoding was.
    """

    key: str
    label: str
    library_id: str
    plate_id: str
    culture_plate_prefix: str
    encodings: tuple[str, ...]
    illumina_groups: tuple[str, ...]


@dataclass(frozen=True)
class Clone:
    """One sequenced clone from one well.

    ``insert_class`` is set only for a full-length run, where the consensus spans
    the whole amplicon but the outcome over the insert alone is also recorded.
    """

    culture_plate: str
    well_id: str
    encoding: str
    block: int
    design_key: str
    grade: str
    kind: str
    sequence_id: str
    insert_class: str | None = None


@dataclass(frozen=True)
class IlluminaRead:
    """One classified Illumina read."""

    group: str
    encoding: str
    block: int
    expected_block: int
    design_key: str
    category: str
    edit_distance: int
    reference_length: int


def parse_reference_name(name: str) -> tuple[str, int, str] | None:
    """``A_Block_9_dTF090_x`` -> ``("A", 9, "dTF090_x")``; None if unparseable."""

    match = REFERENCE_NAME.match(name)
    if not match:
        return None
    return match.group(1) or "", int(match.group(2)), name[match.end() :]


def load_design_universe(fasta_paths: list[Path]) -> set[tuple[str, int, str]]:
    """Every designed reference, as (encoding, block, design key)."""

    universe: set[tuple[str, int, str]] = set()
    for path in fasta_paths:
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                if not line.startswith(">"):
                    continue
                parsed = parse_reference_name(line[1:].split()[0])
                if parsed:
                    universe.add(parsed)
    return universe


def design_keys(universe: set[tuple[str, int, str]], encodings: tuple[str, ...]) -> set[str]:
    """The distinct designs in one or more encodings of a library."""

    return {key for encoding, _, key in universe if encoding in encodings}


def load_insert_outcomes(run_dir: Path) -> dict[str, str]:
    """Per-consensus outcome over the insert region alone, if the run records it.

    A full-length run scores the whole amplicon, but Illumina only ever sees the
    insert, so comparing the two on the amplicon holds nanopore to a standard
    three times longer.  These columns make the like-for-like comparison possible:
    perfect when the insert matches exactly, screenable when it does not but the
    reading frame and stop checks still pass, otherwise other.
    """

    path = run_dir / "stages" / "05_qc" / "qc.csv.gz"
    outcomes: dict[str, str] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if not row.get("insert_edit_distance"):
                continue
            if int(row["insert_edit_distance"]) == 0:
                outcome = PERFECT
            elif row.get("reading_frame") == "pass" and row.get("internal_stops") == "pass":
                outcome = SCREENABLE
            else:
                outcome = OTHER
            outcomes[row["consensus_id"]] = outcome
    return outcomes


def load_clones(run_dir: Path, *, scope: str = "amplicon") -> dict[str, list[Clone]]:
    """Every sequenced clone in the run, keyed by plate barcode.

    Grades come from the QC index, which already carries chimeras alongside
    designed consensuses; sequence identity comes from the consensus table so
    that "distinct sequences" means distinct sequence content, not distinct name.

    ``scope="insert"`` additionally records the insert-only outcome, for
    comparison against a platform that reads the insert alone.
    """

    stages = run_dir / "stages"
    insert_outcomes = load_insert_outcomes(run_dir) if scope == "insert" else {}
    digests: dict[str, str] = {}
    with gzip.open(stages / "04_consensus" / "consensus.csv.gz", "rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            digests[row["consensus_id"]] = row["sequence_sha256"]

    clones: dict[str, list[Clone]] = defaultdict(list)
    index = stages / "05_qc" / "consensus_by_plate" / "index.csv"
    with index.open("r", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["grade"] == NO_SEQUENCE_GRADE:
                continue
            # A chimera is named by its parents; its block is the first parent's.
            parsed = parse_reference_name(row["design"].split(" >> ")[0])
            if not parsed:
                continue
            encoding, block, key = parsed
            clones[row["plate_id"]].append(
                Clone(
                    culture_plate=row["culture_plate"],
                    well_id=row["well_id"],
                    encoding=encoding,
                    block=block,
                    design_key=key,
                    grade=row["grade"],
                    kind=row["kind"],
                    sequence_id=digests.get(row["consensus_id"]) or row["consensus_id"],
                    insert_class=insert_outcomes.get(row["consensus_id"]),
                )
            )
    return clones


def clones_in(clones: dict[str, list[Clone]], library: LibrarySet) -> list[Clone]:
    """The clones belonging to one library set."""

    return [
        clone
        for clone in clones.get(library.plate_id, ())
        if clone.culture_plate.startswith(library.culture_plate_prefix)
        and clone.encoding in library.encodings
    ]


def load_illumina_reads(results_dir: Path) -> list[IlluminaRead]:
    """Read the Illumina per-read classification, joined to reference lengths."""

    lengths: dict[str, int] = {}
    with (results_dir / "reference_counts.tsv").open("r", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            lengths[row["reference_name"]] = int(row["reference_length"])

    reads: list[IlluminaRead] = []
    path = results_dir / "read_classification.tsv.gz"
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            parsed = parse_reference_name(row["reference_name"])
            if not parsed:
                continue
            reads.append(
                IlluminaRead(
                    group=row["expected_group"],
                    encoding=parsed[0],
                    block=parsed[1],
                    expected_block=int(row["expected_block"]),
                    design_key=parsed[2],
                    category=row["category"],
                    edit_distance=int(row["edit_distance"]),
                    reference_length=lengths.get(row["reference_name"], 0),
                )
            )
    return reads


def wells_sequenced_per_block(clones: list[Clone]) -> Counter[int]:
    """How many wells yielded a clone of each block.

    This is the sampling effort the nanopore run actually spent on a block, and
    therefore the number of Illumina reads that matches it.  Wells that yielded
    nothing cannot be attributed to a block and so are not counted - the match
    is to clones sequenced, not to colonies picked.
    """

    per_block: dict[int, set[tuple[str, str]]] = defaultdict(set)
    for clone in clones:
        per_block[clone.block].add((clone.culture_plate, clone.well_id))
    return Counter({block: len(wells) for block, wells in per_block.items()})


def allocated_wells_per_block(
    clones: list[Clone],
    *,
    wells_per_plate: int,
) -> Counter[int]:
    """Allocated wells per block: every well picked, not only those that worked.

    A culture plate holds a fixed number of allocated wells but may carry several
    blocks, and which well got which block is recorded nowhere - only the clones
    that came back say.  So each plate's allocated wells are apportioned among
    its blocks in the ratio the recovered clones show, by largest remainder, and
    every plate contributes exactly `wells_per_plate`.

    This counts the colonies actually picked, which is the sampling effort to
    match.  Using recovered wells instead would credit nanopore's own failures to
    Illumina as reduced depth.
    """

    by_plate: dict[str, dict[int, set[tuple[str, str]]]] = defaultdict(lambda: defaultdict(set))
    for clone in clones:
        by_plate[clone.culture_plate][clone.block].add((clone.culture_plate, clone.well_id))

    allocated: Counter[int] = Counter()
    for blocks in by_plate.values():
        observed = {block: len(wells) for block, wells in blocks.items()}
        total = sum(observed.values())
        if not total:
            continue
        # Largest remainder, so a plate's wells are neither lost nor invented.
        exact = {block: wells_per_plate * count / total for block, count in observed.items()}
        floors = {block: int(value) for block, value in exact.items()}
        remainder = wells_per_plate - sum(floors.values())
        order = sorted(exact, key=lambda block: (-(exact[block] - floors[block]), block))
        for block in order[:remainder]:
            floors[block] += 1
        for block, count in floors.items():
            allocated[block] += count
    return allocated


def subsample_by_block(
    reads: list[IlluminaRead],
    depths: dict[tuple[str, int], int],
    *,
    seed: int,
) -> list[IlluminaRead]:
    """Take `depths[(group, block)]` reads from each block, without replacement.

    Deterministic for a given seed: blocks are drawn in sorted order so the
    sample does not depend on dictionary or file ordering.
    """

    pools: dict[tuple[str, int], list[IlluminaRead]] = defaultdict(list)
    for read in reads:
        # Keyed on the well's expected block, which is the sampling unit; the
        # reference a read matched may sit in another block, and that is a
        # finding rather than a reason to re-bin it.
        pools[(read.group, read.expected_block)].append(read)

    rng = random.Random(seed)
    sample: list[IlluminaRead] = []
    for key in sorted(depths):
        pool = pools.get(key, [])
        wanted = depths[key]
        sample.extend(pool if wanted >= len(pool) else rng.sample(pool, wanted))
    return sample


def distinct_sequences_per_well(
    clones: list[Clone],
    *,
    culture_plates: int,
    wells_per_plate: int,
    maximum: int = 3,
) -> Counter[int]:
    """Histogram of distinct sequences per allocated well, zeros included.

    Wells that yielded nothing are the point of the figure, so the denominator is
    the allocated well count rather than the wells that happen to appear in the
    data.  Counts above `maximum` are folded into `maximum`.
    """

    per_well: dict[tuple[str, str], set[str]] = defaultdict(set)
    for clone in clones:
        per_well[(clone.culture_plate, clone.well_id)].add(clone.sequence_id)

    histogram: Counter[int] = Counter()
    for sequences in per_well.values():
        histogram[min(len(sequences), maximum)] += 1

    allocated = culture_plates * wells_per_plate
    histogram[0] = allocated - sum(histogram.values())
    return histogram


def nanopore_population(clones: list[Clone]) -> Counter[str]:
    """How the consensus sequences themselves fall across the three outcomes.

    Counts sequences, not designs: a design recovered in forty wells contributes
    forty observations, which is the quantity comparable to a read population.
    """

    return Counter(
        clone.insert_class or GRADE_CLASS[clone.grade]
        for clone in clones
        if clone.kind != "chimera"
        and (clone.insert_class or clone.grade in GRADE_CLASS)
    )


def illumina_population(
    reads: list[IlluminaRead],
    universe: set[tuple[str, int, str]],
) -> Counter[str]:
    """How the reads themselves fall across the three outcomes."""

    return Counter(
        CATEGORY_CLASS[read.category]
        for read in reads
        if (read.encoding, read.block, read.design_key) in universe
        and read.category in CATEGORY_CLASS
    )


def population_fractions(population: Counter[str]) -> dict[str, float]:
    """Percentage of observations at each outcome."""

    total = sum(population.values())
    if not total:
        return {name: 0.0 for name in CLASS_ORDER}
    return {name: 100.0 * population.get(name, 0) / total for name in CLASS_ORDER}


def best_class(classes: list[str]) -> str | None:
    """The best outcome among several observations of one design."""

    ranked = [c for c in classes if c in _CLASS_RANK]
    return min(ranked, key=lambda name: _CLASS_RANK[name]) if ranked else None


def nanopore_recovery(clones: list[Clone]) -> dict[str, str]:
    """Best outcome per design across every well that yielded it.

    Uses the insert-only outcome where a clone carries one, so a full-length run
    is not held to a longer standard than the platform it is compared against.
    """

    observed: dict[str, list[str]] = defaultdict(list)
    for clone in clones:
        # A chimera is not an observation of either parent design.
        if clone.kind == "chimera":
            continue
        outcome = clone.insert_class or GRADE_CLASS.get(clone.grade)
        if outcome:
            observed[clone.design_key].append(outcome)
    return {key: best_class(values) for key, values in observed.items() if best_class(values)}


def illumina_recovery(
    reads: list[IlluminaRead],
    universe: set[tuple[str, int, str]],
) -> dict[str, str]:
    """Best outcome per design across the reads matching a library's references.

    Membership is decided by the reference a read actually matched, not by the
    well it came from: a read from a SUMO_A well that matches a B reference is
    evidence about B.  Filtering on the well's group instead let references from
    outside a library count towards it and pushed recovery above 100%.
    """

    observed: dict[str, list[str]] = defaultdict(list)
    for read in reads:
        if (read.encoding, read.block, read.design_key) not in universe:
            continue
        outcome = CATEGORY_CLASS.get(read.category)
        if outcome:
            observed[read.design_key].append(outcome)
    return {key: best_class(values) for key, values in observed.items() if best_class(values)}


def recovery_fractions(recovery: dict[str, str], universe: int) -> dict[str, float]:
    """Percentage of a designed library recovered at each outcome."""

    counts = Counter(recovery.values())
    return {name: 100.0 * counts.get(name, 0) / universe for name in CLASS_ORDER}


def load_nanopore_read_identities(
    run_dir: Path,
    libraries: list[LibrarySet],
) -> dict[str, list[float]]:
    """Per-read alignment identity for the reads behind each library's clones."""

    identities: dict[str, list[float]] = {library.key: [] for library in libraries}
    path = run_dir / "stages" / "03_assignment" / "assignment_calls.csv.gz"
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["assignment_status"] != "assigned_unique" or not row["best_identity"]:
                continue
            for library in libraries:
                if row["plate_id"] != library.plate_id:
                    continue
                if not row["culture_plate"].startswith(library.culture_plate_prefix):
                    continue
                parsed = parse_reference_name(row["reference_ids"])
                if parsed and parsed[0] in library.encodings:
                    identities[library.key].append(float(row["best_identity"]))
    return identities


def illumina_read_identities(
    reads: list[IlluminaRead],
    universe: set[tuple[str, int, str]],
) -> list[float]:
    """Per-read identity over the reference, comparable to the nanopore figure."""

    return [
        1.0 - read.edit_distance / read.reference_length
        for read in reads
        if (read.encoding, read.block, read.design_key) in universe and read.reference_length
    ]


def load_consensus_identities(run_dir: Path, libraries: list[LibrarySet]) -> dict[str, list[float]]:
    """Per-consensus identity to the design, for context against read accuracy."""

    stages = run_dir / "stages"
    qc: dict[str, str] = {}
    with gzip.open(stages / "05_qc" / "qc.csv.gz", "rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["alignment_identity"]:
                qc[row["consensus_id"]] = row["alignment_identity"]

    identities: dict[str, list[float]] = {library.key: [] for library in libraries}
    index = stages / "05_qc" / "consensus_by_plate" / "index.csv"
    with index.open("r", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["grade"] == NO_SEQUENCE_GRADE or row["kind"] == "chimera":
                continue
            value = qc.get(row["consensus_id"])
            if not value:
                continue
            parsed = parse_reference_name(row["design"])
            for library in libraries:
                if row["plate_id"] != library.plate_id:
                    continue
                if not row["culture_plate"].startswith(library.culture_plate_prefix):
                    continue
                if parsed and parsed[0] in library.encodings:
                    identities[library.key].append(float(value))
    return identities


def summarise_identities(values: list[float]) -> dict[str, float]:
    """Mean, median and exact-match rate for a set of identities."""

    if not values:
        return {"n": 0, "mean": 0.0, "median": 0.0, "exact_fraction": 0.0}
    ordered = sorted(values)
    middle = len(ordered) // 2
    median = (
        ordered[middle]
        if len(ordered) % 2
        else (ordered[middle - 1] + ordered[middle]) / 2
    )
    return {
        "n": len(ordered),
        "mean": sum(ordered) / len(ordered),
        "median": median,
        "exact_fraction": sum(1 for value in ordered if value >= 1.0) / len(ordered),
    }
