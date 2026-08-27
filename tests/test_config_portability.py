"""What makes a configuration shareable: named paths and named tuning.

A configuration describes an experiment. Where that experiment's data sits, and
which tuning a platform needs, are not part of the experiment - and writing either
into the file is what makes it unshareable, or worse, what puts somebody's home
directory in a public repository.
"""

from __future__ import annotations

import os
import unittest
from pathlib import Path
from unittest import mock

from nanopore3.config import (
    ConfigError,
    _apply_preset,
    _deep_merge,
    available_presets,
    expand_variables,
)

REPO = Path(__file__).resolve().parents[1]


class VariableExpansionTests(unittest.TestCase):
    def test_a_named_variable_comes_from_the_environment(self) -> None:
        with mock.patch.dict(os.environ, {"DATA_HOME": "/mnt/run7"}):
            self.assertEqual(
                expand_variables("${DATA_HOME}/reads.fastq", "inputs[0].path"),
                "/mnt/run7/reads.fastq",
            )

    def test_several_variables_in_one_path(self) -> None:
        with mock.patch.dict(os.environ, {"A": "/one", "B": "two"}):
            self.assertEqual(expand_variables("${A}/${B}/x.fa", "p"), "/one/two/x.fa")

    def test_an_unset_variable_is_an_error_not_an_empty_string(self) -> None:
        # Expanding to "" would give "/reads.fastq" and a file-not-found three
        # lines later that says nothing about the cause.
        with mock.patch.dict(os.environ, {}, clear=True):
            with self.assertRaises(ConfigError) as caught:
                expand_variables("${NOT_SET_ANYWHERE}/reads.fastq", "inputs[0].path")
        message = str(caught.exception)
        self.assertIn("NOT_SET_ANYWHERE", message)
        self.assertIn("inputs[0].path", message)
        self.assertIn("export", message)

    def test_every_unset_variable_is_named_at_once(self) -> None:
        with mock.patch.dict(os.environ, {}, clear=True):
            with self.assertRaises(ConfigError) as caught:
                expand_variables("${ONE}/${TWO}", "p")
        self.assertIn("ONE", str(caught.exception))
        self.assertIn("TWO", str(caught.exception))

    def test_a_path_without_variables_is_untouched(self) -> None:
        self.assertEqual(expand_variables("/plain/path.fa", "p"), "/plain/path.fa")
        self.assertEqual(expand_variables("cost $5", "p"), "cost $5")

    def test_no_tracked_config_contains_a_machine_specific_path(self) -> None:
        # The reason this module exists. A tracked config with an absolute home
        # directory in it cannot be shared and should never reach a public repo.
        offenders = []
        for config in sorted((REPO / "configs").rglob("*.yaml")):
            if "local" in config.parts:
                continue
            text = config.read_text(encoding="utf-8")
            for marker in ("/Users/", "/home/", "C:\\\\Users", "CloudStorage", "/Volumes/"):
                if marker in text:
                    offenders.append(f"{config.relative_to(REPO)}: {marker}")
        self.assertEqual(offenders, [], f"machine-specific paths found: {offenders}")


class PresetTests(unittest.TestCase):
    def test_the_shipped_presets_are_discoverable(self) -> None:
        self.assertIn("ont-r10-amplicon", available_presets())

    def test_a_preset_supplies_values_the_config_omits(self) -> None:
        merged = _apply_preset({"preset": "ont-r10-amplicon", "run_name": "x"})
        self.assertEqual(merged["consensus"]["maximum_reads"], 300)
        self.assertEqual(merged["run_name"], "x")
        self.assertNotIn("preset", merged)

    def test_the_config_wins_wherever_both_mention_a_key(self) -> None:
        merged = _apply_preset(
            {"preset": "ont-r10-amplicon", "consensus": {"maximum_reads": 50}}
        )
        self.assertEqual(merged["consensus"]["maximum_reads"], 50)
        # and the siblings it did not mention are still inherited
        self.assertEqual(merged["consensus"]["minimum_depth"], 6)

    def test_an_unknown_preset_lists_the_real_ones(self) -> None:
        with self.assertRaises(ConfigError) as caught:
            _apply_preset({"preset": "ont-r9"})
        self.assertIn("ont-r10-amplicon", str(caught.exception))

    def test_no_preset_leaves_the_config_alone(self) -> None:
        original = {"run_name": "x", "consensus": {"minimum_depth": 2}}
        self.assertEqual(_apply_preset(original), original)

    def test_a_list_replaces_rather_than_extends(self) -> None:
        # Extending would make the result depend on which file was written first,
        # and leave no way to remove an inherited item.
        merged = _deep_merge({"a": [1, 2]}, {"a": [3]})
        self.assertEqual(merged["a"], [3])

    def test_merging_reaches_nested_mappings(self) -> None:
        merged = _deep_merge({"a": {"b": 1, "c": 2}}, {"a": {"c": 9}})
        self.assertEqual(merged["a"], {"b": 1, "c": 9})


class ShortConfigEquivalenceTests(unittest.TestCase):
    """The short form of the shipped run must mean the same as the long form."""

    LONG = REPO / "configs" / "runs" / "260608_full_length.yaml"
    SHORT = REPO / "configs" / "runs" / "260608_full_length_short.yaml"

    @unittest.skipUnless(LONG.is_file() and SHORT.is_file(), "shipped configs absent")
    def test_the_two_forms_resolve_to_the_same_settings(self) -> None:
        from nanopore3 import pipeline
        from nanopore3.config import load_config

        with mock.patch.dict(os.environ, {"NANOPORE_RAW": "/placeholder"}):
            long_config = load_config(self.LONG)
            short_config = load_config(self.SHORT)
        long_dict, short_dict = long_config.as_dict(), short_config.as_dict()
        # Everything but how the flanks are *declared* must already agree.
        for key in set(long_dict) | set(short_dict):
            if key == "reference_libraries":
                continue
            self.assertEqual(long_dict.get(key), short_dict.get(key), key)

        try:
            long_flanks = pipeline.resolve_flanks(long_config)
            short_flanks = pipeline.resolve_flanks(short_config)
        except Exception as exc:  # references are not tracked
            self.skipTest(f"shipped references not present: {exc}")
        self.assertEqual(set(long_flanks), set(short_flanks))
        for library, flanks in long_flanks.items():
            other = short_flanks[library]
            for field in ("upstream", "downstream", "left_anchor", "right_anchor"):
                self.assertEqual(
                    getattr(flanks, field), getattr(other, field), f"{library}.{field}"
                )


if __name__ == "__main__":
    unittest.main()
