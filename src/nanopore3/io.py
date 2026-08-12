"""Strict, streaming FASTQ input/output utilities.

The reader intentionally supports the portable four-line FASTQ representation,
not wrapped sequences.  It fails at the first malformed record with a path and
record number rather than silently truncating a run.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass
import gzip
import hashlib
import io
from pathlib import Path
from typing import TextIO

from .provenance import atomic_write_bytes, sha256_file


FASTQ_DNA_IUPAC = frozenset("ACGTRYSWKMBDHVN.-")


class FastqFormatError(ValueError):
    """A FASTQ structural or sequence validation failure."""

    def __init__(self, path: Path, record_index: int, message: str) -> None:
        self.path = path
        self.record_index = record_index
        self.detail = message
        super().__init__(f"{path}: FASTQ record {record_index + 1}: {message}")


@dataclass(frozen=True, slots=True)
class FastqRecord:
    """One validated FASTQ record with a path-independent stable identity."""

    read_uid: str
    name: str
    description: str
    sequence: str
    quality: str
    record_index: int
    source_sha256: str

    @property
    def mean_quality(self) -> float:
        """Mean Phred+33 quality, or 0 for an empty record."""

        if not self.quality:
            return 0.0
        return sum(ord(symbol) - 33 for symbol in self.quality) / len(self.quality)

    def as_fastq(self) -> str:
        """Serialize using the original header description and a normalized plus line."""

        return f"@{self.description}\n{self.sequence}\n+\n{self.quality}\n"


def stable_read_uid(source_sha256: str, record_index: int) -> str:
    """Return a stable UID from a source artifact digest and zero-based record index."""

    if (
        len(source_sha256) != 64
        or any(character not in "0123456789abcdefABCDEF" for character in source_sha256)
    ):
        raise ValueError("source_sha256 must be a 64-character hexadecimal SHA-256")
    if record_index < 0:
        raise ValueError("record_index must be >= 0")
    material = f"{source_sha256.lower()}\0{record_index}".encode("ascii")
    return hashlib.sha256(material).hexdigest()


def _is_gzip(path: Path) -> bool:
    with path.open("rb") as handle:
        magic = handle.read(2)
    if path.name.lower().endswith(".gz") and magic != b"\x1f\x8b":
        raise OSError(f"File has a .gz suffix but is not gzip data: {path}")
    return magic == b"\x1f\x8b"


def open_text_auto(path: str | Path) -> TextIO:
    """Open plain or gzip-compressed ASCII sequence data for reading.

    Compression is detected from gzip magic bytes; a misleading ``.gz`` suffix
    is rejected.  ``newline=None`` handles LF and CRLF consistently.
    """

    source = Path(path).expanduser().resolve(strict=True)
    if not source.is_file():
        raise OSError(f"Not a regular file: {source}")
    if _is_gzip(source):
        return gzip.open(source, mode="rt", encoding="ascii", newline=None)
    return source.open(mode="r", encoding="ascii", newline=None)


def _without_newline(line: str) -> str:
    return line[:-1] if line.endswith("\n") else line


def iter_fastq(
    path: str | Path,
    *,
    source_sha256: str | None = None,
    validate_sequence: bool = True,
) -> Iterator[FastqRecord]:
    """Stream strict four-line FASTQ records from a plain or gzip file.

    If ``source_sha256`` is omitted the physical input is hashed before parsing.
    Passing a digest previously obtained during ingest avoids a second hashing
    pass.  Stable UIDs depend on the physical input digest and record position,
    so duplicate or blank FASTQ identifiers remain unambiguous.
    """

    source = Path(path).expanduser().resolve(strict=True)
    digest = source_sha256.lower() if source_sha256 is not None else sha256_file(source)
    # Validate a supplied digest before opening the potentially large file.
    stable_read_uid(digest, 0)
    try:
        with open_text_auto(source) as handle:
            record_index = 0
            while True:
                header_line = handle.readline()
                if header_line == "":
                    break
                sequence_line = handle.readline()
                plus_line = handle.readline()
                quality_line = handle.readline()
                if "" in (sequence_line, plus_line, quality_line):
                    raise FastqFormatError(
                        source, record_index, "truncated four-line record"
                    )
                header = _without_newline(header_line)
                sequence = _without_newline(sequence_line).upper()
                plus = _without_newline(plus_line)
                quality = _without_newline(quality_line)
                if not header.startswith("@"):
                    raise FastqFormatError(source, record_index, "header must start with '@'")
                description = header[1:]
                if not description or not description.split(maxsplit=1)[0]:
                    raise FastqFormatError(source, record_index, "read identifier is empty")
                name = description.split(maxsplit=1)[0]
                if not sequence:
                    raise FastqFormatError(source, record_index, "sequence is empty")
                if any(character.isspace() for character in sequence):
                    raise FastqFormatError(
                        source, record_index, "sequence contains whitespace"
                    )
                if validate_sequence:
                    invalid = sorted(set(sequence) - FASTQ_DNA_IUPAC)
                    if invalid:
                        raise FastqFormatError(
                            source,
                            record_index,
                            "sequence contains unsupported symbol(s): "
                            + "".join(invalid),
                        )
                if not plus.startswith("+"):
                    raise FastqFormatError(source, record_index, "third line must start with '+'")
                repeated_header = plus[1:]
                if repeated_header:
                    repeated_name = repeated_header.split(maxsplit=1)[0]
                    if repeated_name != name:
                        raise FastqFormatError(
                            source,
                            record_index,
                            "identifier after '+' does not match the '@' identifier",
                        )
                if len(quality) != len(sequence):
                    raise FastqFormatError(
                        source,
                        record_index,
                        f"quality length {len(quality)} does not match sequence length "
                        f"{len(sequence)}",
                    )
                invalid_quality = [
                    character for character in quality if not 33 <= ord(character) <= 126
                ]
                if invalid_quality:
                    raise FastqFormatError(
                        source,
                        record_index,
                        "quality contains symbols outside printable Phred+33 range",
                    )
                yield FastqRecord(
                    read_uid=stable_read_uid(digest, record_index),
                    name=name,
                    description=description,
                    sequence=sequence,
                    quality=quality,
                    record_index=record_index,
                    source_sha256=digest,
                )
                record_index += 1
    except (gzip.BadGzipFile, EOFError, UnicodeDecodeError) as exc:
        raise FastqFormatError(source, locals().get("record_index", 0), str(exc)) from exc


def _render_fastq(records: Iterable[FastqRecord]) -> bytes:
    buffer = io.StringIO(newline="\n")
    for record in records:
        if len(record.sequence) != len(record.quality):
            raise ValueError(
                f"Cannot write {record.name!r}: sequence and quality lengths differ"
            )
        buffer.write(record.as_fastq())
    return buffer.getvalue().encode("ascii")


def write_fastq(
    records: Iterable[FastqRecord],
    path: str | Path,
    *,
    compress: bool | None = None,
) -> Path:
    """Atomically write FASTQ, with deterministic gzip output when requested.

    Gzip defaults from a ``.gz`` suffix.  The gzip header has no filename or
    timestamp, allowing byte-identical output on repeated runs.
    """

    destination = Path(path)
    use_gzip = destination.name.lower().endswith(".gz") if compress is None else compress
    raw = _render_fastq(records)
    if use_gzip:
        output = io.BytesIO()
        with gzip.GzipFile(fileobj=output, mode="wb", filename="", mtime=0) as handle:
            handle.write(raw)
        raw = output.getvalue()
    atomic_write_bytes(destination, raw)
    return destination


__all__ = [
    "FastqFormatError",
    "FastqRecord",
    "iter_fastq",
    "open_text_auto",
    "stable_read_uid",
    "write_fastq",
]
