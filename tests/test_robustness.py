"""Failures a first-time user will actually hit, and what they are told.

Each of these was accepted silently, or reported in a way that pointed somewhere
other than the mistake. A configuration error caught at preflight costs seconds; the
same error caught at stage five costs the run, and often does not name its cause.
"""

from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import yaml

from nanopore3.config import ConfigError, expand_variables, load_config
from nanopore3.pipeline import PipelineError, validate_inputs


class VariableMistakeTests(unittest.TestCase):
    def test_a_variable_that_is_set_but_empty_is_refused(self) -> None:
        # `export SEQ_DATA=` would otherwise give "/reads.fastq" and a
        # file-not-found that says nothing about the cause.
        with mock.patch.dict(os.environ, {"SEQ_DATA": ""}):
            with self.assertRaises(ConfigError) as caught:
                expand_variables("${SEQ_DATA}/reads.fastq", "inputs[0].path")
        self.assertIn("set but empty", str(caught.exception))
        self.assertIn("SEQ_DATA", str(caught.exception))

    def test_whitespace_only_counts_as_empty(self) -> None:
        with mock.patch.dict(os.environ, {"SEQ_DATA": "   "}):
            with self.assertRaises(ConfigError):
                expand_variables("${SEQ_DATA}/reads.fastq", "p")

    def test_a_missing_brace_is_refused_rather_than_left_literal(self) -> None:
        for text in ("$SEQ_DATA/reads.fastq", "${SEQ_DATA/reads.fastq"):
            with mock.patch.dict(os.environ, {"SEQ_DATA": "/data"}):
                with self.assertRaises(ConfigError) as caught:
                    expand_variables(text, "p")
            self.assertIn("braces", str(caught.exception), text)

    def test_a_dollar_that_is_not_a_variable_attempt_is_left_alone(self) -> None:
        # A directory may legitimately contain one; only something that looks like
        # an attempted reference is a mistake.
        for text in ("/data/cost$5/r.fastq", "/data/dir$", "price is $ here"):
            self.assertEqual(expand_variables(text, "p"), text)

    def test_a_variable_whose_value_contains_a_dollar_is_not_flagged(self) -> None:
        # The check runs on the config text, not on the expansion: the value is the
        # user's data.
        with mock.patch.dict(os.environ, {"ODD": "/mnt/we$rd"}):
            self.assertEqual(expand_variables("${ODD}/r.fastq", "p"), "/mnt/we$rd/r.fastq")


class PresetRegressionTests(unittest.TestCase):
    def test_adding_a_preset_does_not_break_a_config_without_run_name(self) -> None:
        # run_name defaults to the config's filename. A line meant to preserve
        # provenance instead injected run_name=None, so adding `preset:` to a
        # working configuration made it fail.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "reads.fastq").write_text("@a\nACGT\n+\n####\n", encoding="ascii")
            (root / "r.fasta").write_text(">r1\n" + "ACGT" * 8 + "\n", encoding="ascii")
            body = {
                "schema_version": 1,
                "preset": "ont-r10-amplicon",
                "output_root": str(root / "runs"),
                "inputs": [{"path": str(root / "reads.fastq"), "sample_id": "s"}],
                "references": {"fasta": str(root / "r.fasta")},
                "library": {
                    "name": "x",
                    "forward_motif": "ACGTACGTACGTACGTACGT",
                    "reverse_motif": "TGCATGCATGCATGCATGCA",
                },
            }
            path = root / "my_experiment.yaml"
            path.write_text(yaml.safe_dump(body), encoding="utf-8")
            self.assertEqual(load_config(path).run_name, "my_experiment")


class RoutingTests(unittest.TestCase):
    """Several libraries and no map means every read is unroutable."""

    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        (self.root / "reads.fastq").write_text("@a\nACGT\n+\n####\n", encoding="ascii")
        (self.root / "a.fasta").write_text(">a1\n" + "ACGT" * 8 + "\n", encoding="ascii")

    def config(self, libraries: int, mapped: bool) -> Path:
        names = ["libA", "libB", "libC"][:libraries]
        body = {
            "schema_version": 1,
            "run_name": "t",
            "output_root": str(self.root / "runs"),
            "inputs": [{"path": str(self.root / "reads.fastq"), "sample_id": "s"}],
            "reference_libraries": {
                name: {"fasta": str(self.root / "a.fasta")} for name in names
            },
            "library": {
                "name": "x",
                "forward_motif": "ACGTACGTACGTACGTACGT",
                "reverse_motif": "TGCATGCATGCATGCATGCA",
            },
        }
        if mapped:
            body["plate_reference_map"] = {
                f"BC{i + 1:02d}": name for i, name in enumerate(names)
            }
        path = self.root / "c.yaml"
        path.write_text(yaml.safe_dump(body), encoding="utf-8")
        return path

    def test_several_libraries_with_no_map_is_refused_at_preflight(self) -> None:
        with self.assertRaises(PipelineError) as caught:
            validate_inputs(load_config(self.config(2, mapped=False)), scan_fastq=False)
        message = str(caught.exception)
        self.assertIn("plate_reference_map", message)
        # And the message shows what to write, not just what is wrong.
        self.assertIn("BC01:", message)

    def test_one_library_needs_no_map(self) -> None:
        validate_inputs(load_config(self.config(1, mapped=False)), scan_fastq=False)

    def test_a_mapped_multi_library_config_is_fine(self) -> None:
        validate_inputs(load_config(self.config(3, mapped=True)), scan_fastq=False)


if __name__ == "__main__":
    unittest.main()
