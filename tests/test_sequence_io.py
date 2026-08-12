from __future__ import annotations

import gzip
from pathlib import Path
import tempfile
import unittest

from nanopore3.io import FastqFormatError, iter_fastq
from nanopore3.references import read_fasta
from nanopore3.sequence import canonical_kmer, normalize_sequence, reverse_complement


class SequenceIoTests(unittest.TestCase):
    def test_iupac_normalization_and_reverse_complement(self) -> None:
        self.assertEqual(normalize_sequence(" acgu\n"), "ACGT")
        self.assertEqual(reverse_complement("ARYN"), "NRYT")
        sequence = "ACGTRYSWKMBDHVN"
        self.assertEqual(reverse_complement(reverse_complement(sequence)), sequence)
        self.assertEqual(canonical_kmer("GTT"), "AAC")

    def test_strict_fastq_rejects_truncation(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.fastq"
            path.write_text("@r\nACGT\n+\n", encoding="ascii")
            with self.assertRaisesRegex(FastqFormatError, "truncated"):
                list(iter_fastq(path))

    def test_fastq_gzip_and_stable_ids(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reads.fastq.gz"
            with gzip.open(path, "wt", encoding="ascii") as handle:
                handle.write("@same description\nACGT\n+\nIIII\n@same\nTGCA\n+\nIIII\n")
            records = list(iter_fastq(path))
            self.assertEqual(len(records), 2)
            self.assertEqual(records[0].name, records[1].name)
            self.assertNotEqual(records[0].read_uid, records[1].read_uid)
            self.assertEqual(
                [record.read_uid for record in iter_fastq(path)],
                [record.read_uid for record in records],
            )

    def test_reference_alias_groups(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "refs.fa"
            path.write_text(">a\nACGT\n>b\nACGT\n>c\nTGCA\n", encoding="ascii")
            bundle = read_fasta(path)
            self.assertEqual(bundle.alias_groups, (("a", "b"),))
            self.assertEqual(bundle.aliases_for("a"), ("a", "b"))


if __name__ == "__main__":
    unittest.main()

