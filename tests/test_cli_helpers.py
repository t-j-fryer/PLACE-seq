"""The CLI helpers that make a long run legible and a trial run cheap."""

from __future__ import annotations

import gzip
import logging
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from nanopore3.cli import _configure_progress, _subsample
from nanopore3.config import ConfigError, ParallelSettings
from nanopore3.pipeline import PipelineError, _duration

RECORD = "@read_{n}\nACGTACGTAC\n+\n##########\n"


class SubsampleTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def fastq(self, count: int, *, gzipped: bool = False) -> Path:
        text = "".join(RECORD.format(n=n) for n in range(count))
        path = self.root / ("reads.fastq.gz" if gzipped else "reads.fastq")
        if gzipped:
            path.write_bytes(gzip.compress(text.encode("ascii")))
        else:
            path.write_text(text, encoding="ascii")
        return path

    def test_the_first_n_reads_are_taken_in_file_order(self) -> None:
        out = self.root / "sub.fastq"
        self.assertEqual(_subsample(self.fastq(50), out, 5), 5)
        headers = [
            line
            for line in out.read_text(encoding="ascii").splitlines()
            if line.startswith("@read_")
        ]
        self.assertEqual(headers, [f"@read_{n}" for n in range(5)])

    def test_asking_for_more_reads_than_exist_takes_what_there_is(self) -> None:
        out = self.root / "sub.fastq"
        self.assertEqual(_subsample(self.fastq(3), out, 100), 3)

    def test_gzip_is_handled_on_both_sides(self) -> None:
        out = self.root / "sub.fastq.gz"
        self.assertEqual(_subsample(self.fastq(20, gzipped=True), out, 4), 4)
        self.assertEqual(
            gzip.decompress(out.read_bytes()).decode("ascii").count("@read_"), 4
        )

    def test_a_file_that_is_not_a_fastq_is_refused(self) -> None:
        bad = self.root / "notes.txt"
        bad.write_text("this is not a fastq\n", encoding="ascii")
        with self.assertRaises(PipelineError) as caught:
            _subsample(bad, self.root / "out.fastq", 1)
        self.assertIn("not a FASTQ", str(caught.exception))

    def test_zero_reads_is_refused_rather_than_writing_nothing(self) -> None:
        with self.assertRaises(PipelineError):
            _subsample(self.fastq(5), self.root / "out.fastq", 0)

    def test_the_destination_directory_is_created(self) -> None:
        out = self.root / "nested" / "deeper" / "sub.fastq"
        _subsample(self.fastq(5), out, 2)
        self.assertTrue(out.is_file())


class ProgressTests(unittest.TestCase):
    def tearDown(self) -> None:
        logger = logging.getLogger("nanopore3")
        logger.handlers.clear()
        logger.propagate = True

    def test_progress_is_scoped_to_this_package(self) -> None:
        # Raising the root logger to INFO turns on every dependency, and
        # matplotlib's font machinery alone buries what this was meant to show.
        _configure_progress(quiet=False)
        self.assertEqual(logging.getLogger("nanopore3").level, logging.INFO)
        self.assertNotEqual(logging.getLogger().level, logging.INFO)

    def test_quiet_leaves_only_warnings(self) -> None:
        _configure_progress(quiet=True)
        self.assertEqual(logging.getLogger("nanopore3").level, logging.WARNING)

    def test_calling_it_twice_does_not_double_the_output(self) -> None:
        _configure_progress(quiet=False)
        _configure_progress(quiet=False)
        handlers = logging.getLogger("nanopore3").handlers
        self.assertEqual(sum(isinstance(h, logging.StreamHandler) for h in handlers), 1)


class DurationTests(unittest.TestCase):
    def test_it_reads_at_every_scale(self) -> None:
        self.assertEqual(_duration(9.4), "9s")
        self.assertEqual(_duration(75), "1m15s")
        self.assertEqual(_duration(3 * 3600 + 25 * 60), "3h25m")


class AutoJobsTests(unittest.TestCase):
    def test_zero_jobs_means_detect(self) -> None:
        # One configuration should run sensibly on a workstation and on a two-core
        # hosted notebook without being edited.
        self.assertEqual(ParallelSettings(jobs=0).jobs, 0)

    def test_negative_jobs_is_still_refused(self) -> None:
        with self.assertRaises(ConfigError):
            ParallelSettings(jobs=-2)


if __name__ == "__main__":
    unittest.main()
