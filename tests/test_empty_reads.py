"""Zero-length basecalls must be rejected loudly or counted explicitly."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from nanopore3.config import ConfigError, LibrarySettings
from nanopore3.io import FastqFormatError, iter_fastq

_FASTQ = """@read_one runid=x
ACGTACGTAC
+
IIIIIIIIII
@read_empty runid=x

+

@read_three runid=x
GGGGCCCCAA
+
IIIIIIIIII
"""


class EmptyReadTests(unittest.TestCase):
    def setUp(self) -> None:
        self._directory = tempfile.TemporaryDirectory()
        self.path = Path(self._directory.name) / "reads.fastq"
        self.path.write_text(_FASTQ, encoding="ascii")

    def tearDown(self) -> None:
        self._directory.cleanup()

    def test_empty_sequence_is_rejected_by_default(self) -> None:
        with self.assertRaises(FastqFormatError) as caught:
            list(iter_fastq(self.path))
        self.assertIn("allow_empty_reads", str(caught.exception))

    def test_empty_sequence_is_retained_when_allowed(self) -> None:
        records = list(iter_fastq(self.path, allow_empty_sequence=True))
        self.assertEqual(
            [record.name for record in records],
            ["read_one", "read_empty", "read_three"],
        )
        self.assertEqual(records[1].sequence, "")
        self.assertEqual(records[1].quality, "")
        self.assertEqual(records[1].mean_quality, 0.0)
        # Reads keep their positions, so a skipped record cannot shift the
        # identity of the reads after it.
        self.assertEqual([record.record_index for record in records], [0, 1, 2])
        self.assertEqual(len({record.read_uid for record in records}), 3)

    def test_a_truncated_final_record_is_still_an_error(self) -> None:
        self.path.write_text("@only_a_header\nACGT\n+\n", encoding="ascii")
        with self.assertRaises(FastqFormatError):
            list(iter_fastq(self.path, allow_empty_sequence=True))

    def test_allow_empty_reads_defaults_to_false_and_must_be_boolean(self) -> None:
        self.assertFalse(LibrarySettings().allow_empty_reads)
        from nanopore3.config import _boolean

        self.assertTrue(_boolean(True, "library.allow_empty_reads"))
        with self.assertRaisesRegex(ConfigError, "must be true or false"):
            _boolean("yes", "library.allow_empty_reads")


if __name__ == "__main__":
    unittest.main()
