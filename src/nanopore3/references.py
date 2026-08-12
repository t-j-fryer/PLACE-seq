"""Strict reference FASTA parsing and duplicate-sequence alias detection."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .io import open_text_auto
from .provenance import canonical_digest, sha256_bytes, sha256_file


REFERENCE_DNA_IUPAC = frozenset("ACGTRYSWKMBDHVN")


class FastaFormatError(ValueError):
    """A FASTA structure, alphabet, or identifier validation failure."""


@dataclass(frozen=True, slots=True)
class ReferenceRecord:
    """One reference sequence with its original description and source."""

    id: str
    description: str
    sequence: str
    sequence_sha256: str
    source_path: Path

    def as_fasta(self, *, line_width: int = 80) -> str:
        """Return a normalized FASTA representation."""

        if line_width < 1:
            raise ValueError("line_width must be >= 1")
        lines = [f">{self.description}"]
        lines.extend(
            self.sequence[index : index + line_width]
            for index in range(0, len(self.sequence), line_width)
        )
        return "\n".join(lines) + "\n"


@dataclass(frozen=True, slots=True)
class ReferenceBundle:
    """Validated references and equivalence groups for identical sequences.

    ``alias_groups`` contains only groups with two or more distinct identifiers.
    Assignment code should report such a group rather than breaking a biological
    tie by filename or lexical ordering.
    """

    records: tuple[ReferenceRecord, ...]
    alias_groups: tuple[tuple[str, ...], ...]
    source_sha256: tuple[tuple[str, str], ...]
    digest: str

    def get(self, reference_id: str) -> ReferenceRecord:
        """Return one reference by ID, raising ``KeyError`` when absent."""

        for record in self.records:
            if record.id == reference_id:
                return record
        raise KeyError(reference_id)

    def aliases_for(self, reference_id: str) -> tuple[str, ...]:
        """Return the identical-sequence equivalence group for an ID."""

        for group in self.alias_groups:
            if reference_id in group:
                return group
        # Validate the ID even when it has no aliases.
        self.get(reference_id)
        return (reference_id,)

    @property
    def ids(self) -> tuple[str, ...]:
        """Reference identifiers in deterministic input order."""

        return tuple(record.id for record in self.records)


def _coerce_paths(paths: str | Path | Iterable[str | Path]) -> tuple[Path, ...]:
    if isinstance(paths, (str, Path)):
        values = (paths,)
    else:
        values = tuple(paths)
    if not values:
        raise FastaFormatError("At least one reference FASTA is required")
    result: list[Path] = []
    for value in values:
        path = Path(value).expanduser().resolve(strict=True)
        if not path.is_file():
            raise FastaFormatError(f"Reference path is not a regular file: {path}")
        result.append(path)
    return tuple(result)


def _records_from_path(path: Path) -> list[ReferenceRecord]:
    records: list[ReferenceRecord] = []
    description: str | None = None
    sequence_parts: list[str] = []
    header_line = 0

    def finish_record() -> None:
        nonlocal description, sequence_parts
        if description is None:
            return
        sequence = "".join(sequence_parts).upper()
        reference_id = description.split(maxsplit=1)[0]
        if not sequence:
            raise FastaFormatError(
                f"{path}:{header_line}: reference {reference_id!r} has no sequence"
            )
        invalid = sorted(set(sequence) - REFERENCE_DNA_IUPAC)
        if invalid:
            raise FastaFormatError(
                f"{path}:{header_line}: reference {reference_id!r} contains "
                f"unsupported symbol(s): {''.join(invalid)}"
            )
        records.append(
            ReferenceRecord(
                id=reference_id,
                description=description,
                sequence=sequence,
                sequence_sha256=sha256_bytes(sequence.encode("ascii")),
                source_path=path,
            )
        )
        description = None
        sequence_parts = []

    try:
        with open_text_auto(path) as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.rstrip("\r\n")
                if not line:
                    continue
                if line.startswith(">"):
                    finish_record()
                    description = line[1:]
                    header_line = line_number
                    if not description or not description.split(maxsplit=1)[0]:
                        raise FastaFormatError(f"{path}:{line_number}: empty FASTA header")
                    if any(character.isspace() for character in description.split(maxsplit=1)[0]):
                        raise FastaFormatError(
                            f"{path}:{line_number}: invalid whitespace in reference ID"
                        )
                    continue
                if description is None:
                    raise FastaFormatError(
                        f"{path}:{line_number}: sequence data occurs before the first header"
                    )
                if any(character.isspace() for character in line):
                    raise FastaFormatError(
                        f"{path}:{line_number}: sequence line contains whitespace"
                    )
                sequence_parts.append(line)
        finish_record()
    except UnicodeDecodeError as exc:
        raise FastaFormatError(f"{path}: FASTA is not ASCII text: {exc}") from exc
    if not records:
        raise FastaFormatError(f"Reference FASTA contains no records: {path}")
    return records


def read_fasta(paths: str | Path | Iterable[str | Path]) -> ReferenceBundle:
    """Read one or more plain/gzip FASTAs into a validated reference bundle.

    Identifiers must be globally unique across all inputs.  Exact duplicate
    sequences are preserved but exposed as alias groups, preventing arbitrary
    assignment to one of several experimentally indistinguishable names.
    """

    source_paths = _coerce_paths(paths)
    records: list[ReferenceRecord] = []
    seen_ids: dict[str, Path] = {}
    source_digests: list[tuple[str, str]] = []
    for path in source_paths:
        source_digests.append((str(path), sha256_file(path)))
        for record in _records_from_path(path):
            previous = seen_ids.get(record.id)
            if previous is not None:
                raise FastaFormatError(
                    f"Duplicate reference ID {record.id!r} in {path}; first seen in {previous}"
                )
            seen_ids[record.id] = path
            records.append(record)

    sequence_owners: dict[str, list[str]] = {}
    for record in records:
        # Key by full sequence, rather than trusting a digest to establish equality.
        sequence_owners.setdefault(record.sequence, []).append(record.id)
    alias_groups = tuple(
        sorted(
            (tuple(sorted(owners)) for owners in sequence_owners.values() if len(owners) > 1),
            key=lambda owners: owners,
        )
    )
    bundle_digest = canonical_digest(
        {
            "schema_version": 1,
            "references": [
                {
                    "id": record.id,
                    "sequence_sha256": record.sequence_sha256,
                    "length": len(record.sequence),
                }
                for record in sorted(records, key=lambda item: item.id)
            ],
            "alias_groups": alias_groups,
        }
    )
    return ReferenceBundle(
        records=tuple(records),
        alias_groups=alias_groups,
        source_sha256=tuple(source_digests),
        digest=bundle_digest,
    )


__all__ = [
    "FastaFormatError",
    "ReferenceBundle",
    "ReferenceRecord",
    "read_fasta",
]
