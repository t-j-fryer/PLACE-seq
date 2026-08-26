"""Browsable, graded consensus output organised by plate and well.

One FASTA per consensus, named with its design and a single-word grade, so a
screening decision can be made from the file listing alone without opening
anything.

Files are filed under ``<plate barcode>/<culture plate>/<well>/`` wherever
compressed-PCR deconvolution resolved a source plate, and under
``<plate barcode>/<well>/`` where it did not.  Pooling several culture plates into
one colony-PCR plate is exactly what makes a bare well coordinate ambiguous, so
once the source plate is known the tree should say so.

Chimeric clones appear alongside the reference-guided consensuses: they are
sequences that are present in the well.  They are graded against the spliced
parent scaffold and carry ``chimera`` in the file name.

The grade collapses the QC criteria into one ordered vocabulary.  Precedence runs
from "no usable data" through "the construct is broken" to "the construct is
fine", so the reported grade is always the most actionable problem rather than
the first one encountered:

``low_depth`` → ``mixed_variants`` → ``frameshift`` → ``premature_stop`` →
``truncated`` → ``mismatched`` → ``perfect`` / ``screenable`` / ``mixed_damage``

``mixed_variants`` is *not* "the well is polyclonal".  A well holding several
designs is the normal case here and simply yields several files.  This grade is
per design: the reads assigned to one design disagree with each other beyond the
support threshold, so that design is not a single clean clone.

``mixed_damage`` is the subset of those whose contested positions all carry the
G:C → T:A signature of 8-oxoguanine.  That lesion is *pre-mutagenic* rather than a
mutation - it pairs with C or with A depending on the replication event - so a
single transformed molecule carrying one produces both sequences inside one
colony, without two plasmids ever meeting.  Measured on this data, that class is
the one shared between independently grown cultures, so it is real heritable
sequence rather than a sequencing artefact, and the damage happened before the
clone existed: during synthesis, amplification or assembly.  Such a clone is
**usable, with a caveat**, which is why it grades above ``mixed_variants`` - and
the file name and header say where the contested position sits, because one in the
vector backbone matters far less than one in the reading frame.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import shutil
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping
from pathlib import Path

# Ordered best to worst for reporting; see the module docstring for precedence.
GRADES = (
    "perfect",
    "screenable",
    "mixed_damage",
    "mixed_variants",
    "mismatched",
    "truncated",
    "premature_stop",
    "frameshift",
    "low_depth",
)

GRADE_DESCRIPTIONS = {
    "perfect": "exact match to the designed reference, in frame, no internal stop",
    "screenable": "full length and in frame with no internal stop, but carries substitutions",
    "mixed_damage": (
        "two alleles at one or more positions, all of them G:C>T:A - the "
        "8-oxoguanine signature of DNA damage that occurred before transformation, "
        "so one molecule carrying the lesion produced both sequences in one colony; "
        "a usable clone, with the caveat named in the file header"
    ),
    "mixed_variants": (
        "reads for this one design disagree beyond the support threshold, so it is "
        "not a single clean clone; a well holding several designs is normal and "
        "simply yields several files"
    ),
    "mismatched": "identity or coverage below the QC floor",
    "truncated": "length outside the configured tolerance",
    "premature_stop": "a stop codon occurs before the end of the reading frame",
    "frameshift": "the assembled reading frame is not a whole number of codons",
    "low_depth": "fewer contributing reads than the configured minimum",
}

_UNSAFE = re.compile(r"[^A-Za-z0-9._-]+")
# Keep paths well inside the Windows MAX_PATH budget; long design names are
# truncated with a short digest so distinct designs never collide.
_MAX_DESIGN_CHARS = 56


def grade_consensus(consensus: Mapping[str, str], qc: Mapping[str, str]) -> str:
    """Reduce one consensus and its QC row to a single actionable grade."""

    if consensus.get("status") == "low_depth" or qc.get("overall") == "not_evaluable":
        return "low_depth"
    if int(consensus.get("ambiguous_bases") or 0) > 0:
        # A sample problem rather than a sequence one, and it explains any
        # identity or frame failure that follows from it. Damage-signature
        # mixtures are separated out because they are real, usable clones: see
        # the module docstring.
        if qc.get("mixed_signature") == "oxidative":
            return "mixed_damage"
        return "mixed_variants"
    if qc.get("reading_frame") == "fail":
        return "frameshift"
    if qc.get("internal_stops") == "fail":
        return "premature_stop"
    if qc.get("expected_length") == "fail":
        return "truncated"
    if qc.get("full_amplicon") == "fail":
        return "mismatched"
    distance = qc.get("alignment_edit_distance")
    if distance not in (None, "") and int(distance) == 0:
        return "perfect"
    return "screenable"


def safe_name(value: str, *, limit: int = _MAX_DESIGN_CHARS) -> str:
    """Return a filesystem-safe fragment, shortened stably when overlong."""

    cleaned = _UNSAFE.sub("_", value).strip("._-") or "unnamed"
    if len(cleaned) <= limit:
        return cleaned
    digest = hashlib.sha256(value.encode("utf-8")).hexdigest()[:8]
    return f"{cleaned[: limit - 9]}~{digest}"


def normalized_well(well: str) -> str:
    """Zero-pad a well so a directory listing sorts A01 before A10."""

    match = re.fullmatch(r"([A-Za-z]+)\s*0*(\d+)", well.strip())
    if not match:
        return safe_name(well)
    return f"{match.group(1).upper()}{int(match.group(2)):02d}"


def _read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    current: str | None = None
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith(">"):
            current = line[1:].split()[0]
            sequences[current] = ""
        elif current is not None:
            sequences[current] += line.strip()
    return sequences


def write_consensus_tree(
    consensus_rows: Iterable[Mapping[str, str]],
    qc_by_id: Mapping[str, Mapping[str, str]],
    sequences: Mapping[str, str],
    root: Path,
) -> dict[str, object]:
    """Write one graded FASTA per consensus under ``<plate>/<well>/``.

    A consensus with no sequence (too few reads to build one) is still recorded
    in the index, so a missing file never has to be interpreted as an oversight.

    The tree is built beside the target and swapped in, so a rerun cannot leave
    files from a previous build behind: a stale FASTA is indistinguishable from a
    current one once written, and would be read as a real result.  An existing
    directory is only replaced when it carries this function's ``index.csv``,
    so an unrelated directory is never deleted.
    """

    if root.exists():
        if not (root / "index.csv").is_file():
            raise FileExistsError(
                f"{root} exists and was not written by this exporter "
                "(no index.csv); refusing to replace it"
            )
    staging = root.with_name(root.name + ".partial")
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True)
    written_root, root = root, staging
    index_rows: list[dict[str, object]] = []
    per_plate: dict[str, Counter[str]] = defaultdict(Counter)
    totals: Counter[str] = Counter()
    used_paths: set[Path] = set()

    def order(item: Mapping[str, str]) -> tuple[str, str, str, str]:
        return (
            item["plate_id"],
            item.get("culture_plate", ""),
            normalized_well(item["well_id"]),
            item["reference_ids"],
        )

    for row in sorted(consensus_rows, key=order):
        consensus_id = row["consensus_id"]
        qc = qc_by_id.get(consensus_id, {})
        grade = grade_consensus(row, qc)
        totals[grade] += 1
        per_plate[row["plate_id"]][grade] += 1

        aliases = [alias for alias in row["reference_ids"].split("|") if alias]
        design = aliases[0] if aliases else "unassigned"
        if len(aliases) > 1:
            design = f"{design}+{len(aliases) - 1}"
        plate = safe_name(row["plate_id"] or "unknown_plate")
        well = normalized_well(row["well_id"] or "unknown_well")
        # A single culture plate: a "|"-joined set means the read could not be
        # attributed to one, so it stays at well level rather than being filed
        # under a plate it may not have come from.
        source = row.get("culture_plate", "")
        culture = safe_name(source) if source and "|" not in source else ""
        sequence = sequences.get(consensus_id, "")

        relative: str | None = None
        if sequence:
            directory = root / plate / culture / well if culture else root / plate / well
            directory.mkdir(parents=True, exist_ok=True)
            marker = "chimera__" if row.get("status") == "chimera" else ""
            prefix = f"{plate}_{culture}_{well}" if culture else f"{plate}_{well}"
            stem = f"{prefix}__{marker}{safe_name(design)}__{grade}"
            path = directory / f"{stem}.fasta"
            if path in used_paths:  # distinct designs that shorten alike
                path = directory / f"{stem}__{consensus_id[-8:]}.fasta"
            used_paths.add(path)
            header = (
                f">{consensus_id} grade={grade} "
                + ("kind=chimera " if row.get("status") == "chimera" else "")
                + f"design={row['reference_ids']} "
                f"library={row['reference_library_id']} plate={row['plate_id']} "
                f"well={row['well_id']}"
            )
            if row.get("culture_plate"):
                header += f" culture_plate={row['culture_plate']}"
            header += (
                f" reads_used={row['n_reads_used']} mean_depth={row['mean_depth']}"
                f" ambiguous_bases={row['ambiguous_bases']}"
            )
            if qc.get("alignment_identity"):
                header += f" identity={qc['alignment_identity']}"
            if qc.get("mixed_detail"):
                # Where the contested position is, and what it does to the
                # protein: a mixture in the vector backbone matters far less than
                # one in the reading frame, and the grade alone cannot say which.
                header += (
                    f" mixed_signature={qc.get('mixed_signature', '')}"
                    f" mixed_at={qc['mixed_detail']}"
                    f" mixed_in_reading_frame={qc.get('mixed_in_reading_frame', '')}"
                    f" mixed_worst_effect={qc.get('mixed_worst_effect', '')}"
                    f" designed_allele_fraction={qc.get('designed_allele_fraction', '')}"
                )
            # 60-column wrapping keeps the files readable in any viewer.
            wrapped = "\n".join(sequence[i : i + 60] for i in range(0, len(sequence), 60))
            path.write_text(f"{header}\n{wrapped}\n", encoding="ascii")
            relative = str(path.relative_to(root))

        index_rows.append(
            {
                "grade": grade,
                "kind": "chimera" if row.get("status") == "chimera" else "consensus",
                "plate_id": row["plate_id"],
                "well_id": row["well_id"],
                "design": row["reference_ids"],
                "reference_library_id": row["reference_library_id"],
                "culture_plate": row.get("culture_plate", ""),
                "consensus_id": consensus_id,
                "reads_used": row["n_reads_used"],
                "mean_depth": row["mean_depth"],
                "ambiguous_bases": row["ambiguous_bases"],
                # Per-region accuracy, where a full-length library reports it.
                # "Is the error in the part I designed, or in the vector?" is the
                # first question a screener asks, and a whole-amplicon identity
                # cannot answer it.
                "insert_identity": qc.get("insert_identity", ""),
                "insert_edit_distance": qc.get("insert_edit_distance", ""),
                "flank_5p_edit_distance": qc.get("flank_5p_edit_distance", ""),
                "flank_3p_edit_distance": qc.get("flank_3p_edit_distance", ""),
                "mixed_signature": qc.get("mixed_signature", ""),
                "mixed_detail": qc.get("mixed_detail", ""),
                "mixed_in_reading_frame": qc.get("mixed_in_reading_frame", ""),
                "mixed_worst_effect": qc.get("mixed_worst_effect", ""),
                "designed_allele_fraction": qc.get("designed_allele_fraction", ""),
                "identity": qc.get("alignment_identity", ""),
                "edit_distance": qc.get("alignment_edit_distance", ""),
                "protein_length": qc.get("protein_length", ""),
                "qc_reason": qc.get("reason", row.get("failure_reason", "")),
                "file": relative or "",
            }
        )

    fields = list(index_rows[0]) if index_rows else ["grade", "plate_id", "well_id"]
    with (root / "index.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(index_rows)

    summary = {
        "grades": {grade: totals[grade] for grade in GRADES if totals[grade]},
        "by_plate": {
            plate: {grade: counts[grade] for grade in GRADES if counts[grade]}
            for plate, counts in sorted(per_plate.items())
        },
        "descriptions": {g: GRADE_DESCRIPTIONS[g] for g in GRADES if totals[g]},
        "files_written": sum(1 for row in index_rows if row["file"]),
        "consensuses": len(index_rows),
    }
    (root / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if written_root.exists():
        shutil.rmtree(written_root)
    root.rename(written_root)
    return summary


def write_tree_from_stages(
    consensus_dir: Path, qc_csv: Path, root: Path
) -> dict[str, object]:
    """Build the graded tree from a completed run's consensus and QC tables."""

    import gzip

    path = consensus_dir / "consensus.csv.gz"
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        consensus_rows = list(csv.DictReader(handle))
    with gzip.open(qc_csv, "rt", encoding="utf-8", newline="") as handle:
        qc_by_id = {row["consensus_id"]: row for row in csv.DictReader(handle)}
    sequences = _read_fasta(consensus_dir / "consensus.fasta")
    return write_consensus_tree(consensus_rows, qc_by_id, sequences, root)


__all__ = [
    "GRADES",
    "GRADE_DESCRIPTIONS",
    "grade_consensus",
    "normalized_well",
    "safe_name",
    "write_consensus_tree",
    "write_tree_from_stages",
]
