"""Temporary disk-backed grouping with stable Python tuple ordering.

SQLite holds the dataset and output accumulators; callers materialize only one
group. Its cache is fixed at 2 MiB and sort temporaries stay on disk. Files are
removed before the enclosing immutable stage is published, including on errors.
"""

from __future__ import annotations

import json
import sqlite3
from collections.abc import Iterable, Iterator
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

from .errors import Nanopore3Error


class GroupMemoryError(Nanopore3Error):
    """A group exceeds the configured retained-memory estimate budget."""


class GroupStorageError(Nanopore3Error):
    """Temporary grouping storage cannot be read or written."""


def check_group_memory(estimated_bytes: int, limit: int, key: object) -> None:
    if estimated_bytes > limit:
        raise GroupMemoryError(
            f"group {key!r} exceeds the retained-memory estimate budget "
            f"({estimated_bytes / 1048576:.1f} MiB > {limit / 1048576:.1f} MiB). "
            "No reads were silently dropped. Increase parallel.group_memory_mb "
            "on a machine with sufficient memory, or split the input. "
            "Cgroup memory limits also constrain this budget."
        )


def _tuple(value: Any) -> Any:
    return tuple(_tuple(item) for item in value) if isinstance(value, list) else value


def _compare(left: str, right: str) -> int:
    a, b = _tuple(json.loads(left)), _tuple(json.loads(right))
    return (a > b) - (a < b)


class DiskStore:
    def __init__(self, directory: Path) -> None:
        directory.mkdir(parents=True, exist_ok=True)
        self._temporary = TemporaryDirectory(prefix=".groups-", dir=directory)
        self._cursors: set[sqlite3.Cursor] = set()
        try:
            self.db = sqlite3.connect(str(Path(self._temporary.name) / "groups.sqlite"))
            self.db.create_collation("PYKEY", _compare)
            self.db.executescript("""
                PRAGMA journal_mode=OFF;
                PRAGMA synchronous=OFF;
                PRAGMA cache_size=-2048;
                PRAGMA temp_store=FILE;
                PRAGMA mmap_size=0;
                CREATE TABLE records (
                    ordinal INTEGER PRIMARY KEY, bucket TEXT, key TEXT, payload TEXT
                );
                CREATE INDEX groups_idx ON records(bucket, key, ordinal);
            """)
        except BaseException as exc:
            if hasattr(self, "db"):
                self.db.close()
            self._temporary.cleanup()
            if isinstance(exc, sqlite3.Error):
                raise GroupStorageError(f"cannot initialize grouping storage: {exc}") from exc
            raise

    def __enter__(self) -> DiskStore:
        return self

    def __exit__(self, *exc: object) -> None:
        # Suspended generators can retain live cursors after a memory guard
        # raises. Close them before unlinking the database, including on Windows.
        try:
            for cursor in self._cursors:
                cursor.close()
            self._cursors.clear()
            self.db.close()
        finally:
            self._temporary.cleanup()
        if len(exc) > 1 and isinstance(exc[1], sqlite3.Error):
            raise GroupStorageError(
                f"temporary grouping storage failed: {exc[1]}; check free disk space "
                "and filesystem access in the output directory"
            ) from exc[1]

    def _query(self, sql: str, parameters: tuple) -> Iterator[tuple]:
        cursor = self.db.execute(sql, parameters)
        self._cursors.add(cursor)
        try:
            # Explicit iteration avoids delegating generator.close() to an
            # already-closed SQLite cursor after the store's context has exited.
            for row in cursor:  # noqa: UP028 -- cursor.close() delegation is unsafe here
                yield row
        finally:
            if cursor in self._cursors:
                self._cursors.remove(cursor)
                cursor.close()

    def append(self, bucket: str, payload: Any, *, key: Any = ()) -> None:
        self.db.execute(
            "INSERT INTO records(bucket, key, payload) VALUES (?, ?, ?)",
            (bucket, json.dumps(key), json.dumps(payload)),
        )

    def extend(self, bucket: str, values: Iterable[Any]) -> None:
        self.db.executemany(
            "INSERT INTO records(bucket, key, payload) VALUES (?, '[]', ?)",
            ((bucket, json.dumps(value)) for value in values),
        )

    def contains(self, bucket: str, key: Any) -> bool:
        return (
            self.db.execute(
                "SELECT 1 FROM records WHERE bucket=? AND key=? LIMIT 1", (bucket, json.dumps(key))
            ).fetchone()
            is not None
        )

    def count(self, bucket: str) -> int:
        return self.db.execute("SELECT COUNT(*) FROM records WHERE bucket=?", (bucket,)).fetchone()[
            0
        ]

    def values(self, bucket: str, *, key: Any = ()) -> Iterator[Any]:
        cursor = self._query(
            "SELECT payload FROM records WHERE bucket=? AND key=? ORDER BY ordinal",
            (bucket, json.dumps(key)),
        )
        for (payload,) in cursor:
            yield json.loads(payload)

    def groups(self, bucket: str) -> Iterator[tuple[Any, Iterator[Any]]]:
        # Sort distinct group keys, preserving the old Python tuple order even
        # for aliases containing prefixes, quotes, non-ASCII or backslashes.
        cursor = self._query(
            "SELECT DISTINCT key FROM records WHERE bucket=? ORDER BY key COLLATE PYKEY", (bucket,)
        )
        for (encoded,) in cursor:
            key = _tuple(json.loads(encoded))
            yield key, self.values(bucket, key=key)

    def sorted_values(self, bucket: str) -> Iterator[Any]:
        for (payload,) in self._query(
            "SELECT payload FROM records WHERE bucket=? ORDER BY payload COLLATE PYKEY", (bucket,)
        ):
            yield json.loads(payload)
