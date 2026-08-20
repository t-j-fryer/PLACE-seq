"""Strict reference FASTA parsing and duplicate-sequence alias detection."""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Callable, Iterable, Mapping

from .io import open_text_auto
from .provenance import canonical_digest, sha256_bytes, sha256_file


REFERENCE_DNA_IUPAC = frozenset("ACGTRYSWKMBDHVN")


class FastaFormatError(ValueError):
    """A FASTA structure, alphabet, or identifier validation failure."""


ALIAS_SEPARATOR = "|"


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


@dataclass(frozen=True, slots=True)
class ReferenceLibraryCollection:
    """A deterministic collection of independently validated reference sets.

    Identifiers must be unique within a library, while the same biological ID
    may intentionally occur in two libraries that are never searched together.
    ``libraries`` is stored as sorted pairs to make iteration and provenance
    independent of input mapping order.
    """

    libraries: tuple[tuple[str, ReferenceBundle], ...]
    digest: str

    @property
    def ids(self) -> tuple[str, ...]:
        """Library identifiers in deterministic lexical order."""

        return tuple(library_id for library_id, _ in self.libraries)

    def get(self, library_id: str) -> ReferenceBundle:
        """Return one library bundle, raising ``KeyError`` when absent."""

        for candidate, bundle in self.libraries:
            if candidate == library_id:
                return bundle
        raise KeyError(library_id)


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
    if len(set(result)) != len(result):
        raise FastaFormatError(
            "A reference FASTA path must not occur more than once in one library"
        )
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


def read_fasta(
    paths: str | Path | Iterable[str | Path],
    *,
    transform: Callable[[str, str], str] | None = None,
) -> ReferenceBundle:
    """Read one or more plain/gzip FASTAs into a validated reference bundle.

    Identifiers must be globally unique across all inputs.  Exact duplicate
    sequences are preserved but exposed as alias groups, preventing arbitrary
    assignment to one of several experimentally indistinguishable names.

    ``transform`` rewrites each sequence, given its ID, as it is read - used to
    join constant flanks onto insert-only references.  It runs before alias
    detection and digesting, so the bundle describes the sequences that are
    actually searched rather than the file on disk.  Flanking cannot merge two
    distinct inserts or split two identical ones, so alias groups are unaffected.
    """

    source_paths = _coerce_paths(paths)
    records: list[ReferenceRecord] = []
    seen_ids: dict[str, Path] = {}
    source_digests: list[tuple[str, str]] = []
    for path in source_paths:
        source_digests.append((str(path), sha256_file(path)))
        for record in _records_from_path(path):
            if transform is not None:
                sequence = transform(record.id, record.sequence)
                record = replace(
                    record,
                    sequence=sequence,
                    sequence_sha256=sha256_bytes(sequence.encode("ascii")),
                )
            # "|" joins alias groups in every downstream table, and those tables
            # are split back on it. A reference whose own ID contains "|" would
            # survive assignment and only fail when the group was re-parsed, so
            # it is rejected here rather than several stages later.
            if ALIAS_SEPARATOR in record.id:
                raise FastaFormatError(
                    f"Reference ID {record.id!r} in {path} contains "
                    f"{ALIAS_SEPARATOR!r}, which is reserved as the alias-group "
                    "separator in the pipeline's output tables"
                )
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


def read_reference_libraries(
    libraries: Mapping[str, str | Path | Iterable[str | Path]],
    *,
    transforms: Mapping[str, Callable[[str, str], str]] | None = None,
) -> ReferenceLibraryCollection:
    """Read named reference libraries without merging their namespaces.

    Each value accepts the same path forms as :func:`read_fasta`.  Validation
    (including duplicate paths and duplicate record IDs) is performed within
    each library.  Library identifiers and output ordering are normalized for
    stable manifests and fingerprints.
    """

    if not isinstance(libraries, Mapping) or not libraries:
        raise FastaFormatError(
            "At least one named reference library is required"
        )
    bundles: list[tuple[str, ReferenceBundle]] = []
    for raw_library_id in sorted(libraries, key=str):
        if not isinstance(raw_library_id, str) or not raw_library_id.strip():
            raise FastaFormatError(
                "Reference library identifiers must be non-empty strings"
            )
        library_id = raw_library_id.strip()
        if library_id != raw_library_id:
            raise FastaFormatError(
                f"Reference library identifier {raw_library_id!r} has surrounding whitespace"
            )
        bundles.append(
            (
                library_id,
                read_fasta(
                    libraries[raw_library_id],
                    transform=(transforms or {}).get(library_id),
                ),
            )
        )
    collection_digest = canonical_digest(
        {
            "schema_version": 1,
            "reference_libraries": [
                {"id": library_id, "digest": bundle.digest}
                for library_id, bundle in bundles
            ],
        }
    )
    return ReferenceLibraryCollection(
        libraries=tuple(bundles),
        digest=collection_digest,
    )


__all__ = [
    "FastaFormatError",
    "ReferenceBundle",
    "ReferenceLibraryCollection",
    "ReferenceRecord",
    "read_fasta",
    "read_reference_libraries",
]
