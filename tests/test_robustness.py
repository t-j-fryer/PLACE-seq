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


class DoubleFlankGuardTests(unittest.TestCase):
    """Flanking a sequence that already has the flanks must be refused.

    Joining the constant regions on a second time is silent - the reference simply
    grows - and it produced a real defect: flanks derived *from* a set of assembled
    references were then joined *onto* those same references. The guard belongs at
    the point of work, so every caller inherits it.
    """

    UP = "TTAGGCATCCGAATTGCACGATTGCCATAGGTTCAGGATCCAATTGCACGTAGCTTAAGG"
    DOWN = "CCATTGACCTAGGTTCAGCATGCATTGACCGGTTAAGGCCATTAGCTAGGATCCTTAAGC"
    INSERT = "ACGTTGCAAGGTCCATTAGCAA"

    def flanks(self):
        from nanopore3.flanks import from_sequences

        return from_sequences(self.UP, self.DOWN, anchor_length=20)

    def test_an_insert_is_flanked_normally(self) -> None:
        built = self.flanks().flank(self.INSERT)
        self.assertEqual(len(built), len(self.INSERT) + 2 * (60 - 20))

    def test_a_reference_this_already_built_is_refused(self) -> None:
        from nanopore3.flanks import FlankError

        flanks = self.flanks()
        with self.assertRaises(FlankError) as caught:
            flanks.flank(flanks.flank(self.INSERT))
        self.assertIn("already carries", str(caught.exception))

    def test_an_assembled_construct_is_refused_and_points_at_derive(self) -> None:
        from nanopore3.flanks import FlankError

        with self.assertRaises(FlankError) as caught:
            self.flanks().flank(self.UP + self.INSERT + self.DOWN)
        self.assertIn("flanks.derive", str(caught.exception))

    def test_the_derive_route_is_unaffected(self) -> None:
        from nanopore3.flanks import from_assembled

        inserts = (self.INSERT, "TTGACCAGTTCAGGATCCAACC", "GGCATTACCGTAAGGTTCCATT")
        assembled = [self.UP + i + self.DOWN for i in inserts]
        derived = from_assembled(assembled, anchor_length=20, minimum_constant=40)
        # trim, not flank: 20 nt of primer anchor off each end
        self.assertEqual(len(derived.transform(assembled[0])), len(assembled[0]) - 40)


class NothingRecoveredTests(unittest.TestCase):
    """A run that recovers nothing must say why, not raise from a figure."""

    def test_figures_report_a_reason_instead_of_dividing_by_zero(self) -> None:
        # An empty input, or read-length limits that match nothing, left
        # plate_occupancy_figure computing rows from zero columns.
        try:
            from nanopore3.figures import NothingToPlot, plate_occupancy_figure
        except ImportError:
            self.skipTest("matplotlib not installed")
        with self.assertRaises(NothingToPlot):
            plate_occupancy_figure([], Path("unused"))

    def test_a_stage_with_no_successes_is_reported_loudly(self) -> None:
        import logging

        from nanopore3.pipeline import _report_stage_counts

        logger = logging.getLogger("nanopore3")
        with self.assertLogs(logger, level="WARNING") as captured:
            _report_stage_counts("02_demux", {"out_of_length": 2}, "assigned")
        joined = "\n".join(captured.output)
        self.assertIn("no assigned reads", joined)
        self.assertIn("read-length", joined)

    def test_a_stage_with_successes_is_reported_quietly(self) -> None:
        import logging

        from nanopore3.pipeline import _report_stage_counts

        logger = logging.getLogger("nanopore3")
        with self.assertLogs(logger, level="INFO") as captured:
            _report_stage_counts("02_demux", {"assigned": 7, "no_match": 1}, "assigned")
        self.assertFalse([m for m in captured.output if m.startswith("WARNING")])

    def test_an_empty_stage_says_nothing_reached_it(self) -> None:
        import logging

        from nanopore3.pipeline import _report_stage_counts

        with self.assertLogs(logging.getLogger("nanopore3"), level="WARNING") as caught:
            _report_stage_counts("03_assignment", {}, "assigned_unique")
        self.assertIn("no reads reached", "\n".join(caught.output))


class ResolvedJobsTests(unittest.TestCase):
    """`jobs: 0` must survive preflight, not only the stage that uses it.

    It was resolved in two places: the run used the detected CPU count, and
    preflight reported the plan from the raw value and raised on 0. The feature
    passed its unit test and failed before the run began.
    """

    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        (self.root / "reads.fastq").write_text("@a\nACGT\n+\n####\n", encoding="ascii")
        (self.root / "r.fasta").write_text(">r1\n" + "ACGT" * 8 + "\n", encoding="ascii")

    def config(self, jobs: int, backend: str) -> Path:
        body = {
            "schema_version": 1,
            "run_name": "t",
            "output_root": str(self.root / "runs"),
            "inputs": [{"path": str(self.root / "reads.fastq"), "sample_id": "s"}],
            "references": {"fasta": str(self.root / "r.fasta")},
            "library": {
                "name": "x",
                "forward_motif": "ACGTACGTACGTACGTACGT",
                "reverse_motif": "TGCATGCATGCATGCATGCA",
            },
            "parallel": {"jobs": jobs, "backend": backend},
        }
        path = self.root / f"c_{jobs}_{backend}.yaml"
        path.write_text(yaml.safe_dump(body), encoding="utf-8")
        return path

    def test_preflight_accepts_zero_jobs(self) -> None:
        # The whole point: this is what failed, while the unit test passed.
        report = validate_inputs(load_config(self.config(0, "process")), scan_fastq=False)
        self.assertGreaterEqual(report["resources"]["jobs"], 1)

    def test_zero_jobs_resolves_to_the_machine_for_every_parallel_backend(self) -> None:
        from nanopore3.pipeline import resolved_jobs
        from nanopore3.runtime import effective_cpu_count

        for backend in ("auto", "thread", "process"):
            self.assertEqual(
                resolved_jobs(load_config(self.config(0, backend))),
                effective_cpu_count(),
                backend,
            )

    def test_serial_is_one_worker_whatever_jobs_says(self) -> None:
        from nanopore3.pipeline import resolved_jobs

        for jobs in (0, 1, 8):
            self.assertEqual(resolved_jobs(load_config(self.config(jobs, "serial"))), 1)

    def test_an_explicit_count_is_respected(self) -> None:
        from nanopore3.pipeline import resolved_jobs

        self.assertEqual(resolved_jobs(load_config(self.config(3, "process"))), 3)


class UnknownIsNotZeroTests(unittest.TestCase):
    """A missing measurement must not be reported as a measurement of zero.

    `designed_allele_fraction` answers "is the designed sequence still in this
    well". A run recorded before the read fractions existed has the position but
    not its shares, and 0.0 there says "the design is gone" - the opposite
    conclusion, drawn from absent data.
    """

    REFERENCE = "ATG" + "AAACCCGGGTTTTGGCATGAC" + "TAA"

    def test_no_recorded_shares_reports_unknown(self) -> None:
        from nanopore3.qc import classify_mixed_positions, designed_allele_fraction

        positions = classify_mixed_positions(((9, "G", "GT"),), self.REFERENCE)
        self.assertIsNone(positions[0].designed_fraction)
        self.assertIsNone(designed_allele_fraction(positions))

    def test_the_design_being_absent_still_reports_zero(self) -> None:
        from nanopore3.qc import classify_mixed_positions

        # Neither allele is the designed base: measured, and it is gone.
        positions = classify_mixed_positions(
            ((9, "G", "AT", (0.55, 0.45)),), self.REFERENCE
        )
        self.assertEqual(positions[0].designed_fraction, 0.0)

    def test_one_unknown_position_makes_the_clone_unknown(self) -> None:
        from nanopore3.qc import classify_mixed_positions, designed_allele_fraction

        positions = classify_mixed_positions(
            ((9, "G", "GT", (0.62, 0.38)), (12, "T", "CT")), self.REFERENCE
        )
        # Not 0.62: a confident number from partial data is worse than none.
        self.assertIsNone(designed_allele_fraction(positions))

    def test_an_old_consensus_column_still_parses(self) -> None:
        from nanopore3.pipeline import parse_mixed_alleles

        self.assertEqual(parse_mixed_alleles("12:G>GT"), ((12, "G", "GT", ()),))


class ErrorPresentationTests(unittest.TestCase):
    """Every error meant for a person must reach them as a sentence.

    The command line caught a hand-written tuple of exception types, so anything
    not on the list arrived as a traceback. Two did: an unknown `rerun --from`
    stage, and choosing a consensus backend whose native tools are absent. Both had
    good messages; both were presented as a crash.
    """

    def test_every_user_facing_error_shares_one_base(self) -> None:
        from nanopore3.barcodes import BarcodeRegistryError
        from nanopore3.config import ConfigError
        from nanopore3.consensus import (
            ConsensusBackendError,
            ConsensusBackendUnavailable,
        )
        from nanopore3.deconvolution import DeconvolutionError
        from nanopore3.demux import BarcodeValidationError
        from nanopore3.errors import Nanopore3Error
        from nanopore3.flanks import FlankError
        from nanopore3.fragments import FragmentError
        from nanopore3.io import FastqFormatError
        from nanopore3.layout import LayoutError
        from nanopore3.pipeline import PipelineError
        from nanopore3.provenance import ProvenanceError
        from nanopore3.references import FastaFormatError
        from nanopore3.rerun import RerunError
        from nanopore3.sequence import SequenceValidationError

        for error in (
            BarcodeRegistryError, BarcodeValidationError, ConfigError,
            ConsensusBackendError, ConsensusBackendUnavailable, DeconvolutionError,
            FastaFormatError, FastqFormatError, FlankError, FragmentError,
            LayoutError, PipelineError, ProvenanceError, RerunError,
            SequenceValidationError,
        ):
            self.assertTrue(
                issubclass(error, Nanopore3Error), f"{error.__name__} is not caught"
            )

    def test_the_natural_bases_are_preserved(self) -> None:
        # Existing `except ValueError` handlers, in this package and in callers,
        # must keep working.
        from nanopore3.config import ConfigError
        from nanopore3.flanks import FlankError
        from nanopore3.rerun import RerunError

        for error in (ConfigError, FlankError, RerunError):
            self.assertTrue(issubclass(error, ValueError), error.__name__)

    def test_the_cli_catches_the_base_not_a_list(self) -> None:
        import inspect

        from nanopore3 import cli

        source = inspect.getsource(cli.main)
        self.assertIn("Nanopore3Error", source)
        # The old tuple is what let two error classes through.
        self.assertNotIn("except (ConfigError, PipelineError, OSError, ValueError)", source)

    def test_an_unknown_rerun_stage_lists_the_real_ones(self) -> None:
        from nanopore3.errors import Nanopore3Error
        from nanopore3.rerun import stages_before

        with self.assertRaises(Nanopore3Error) as caught:
            stages_before("04_consensuss")
        self.assertIn("04_consensus", str(caught.exception))


class PlateRoutingTests(unittest.TestCase):
    """An explicit plate map says which barcodes belong to the analysis.

    The fallback to a sole reference library used to apply even when a map was
    given, so a configuration naming one barcode silently analysed every barcode on
    the flow cell. On a real run that was eight barcodes and 42% of the reported
    reads.
    """

    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        (self.root / "reads.fastq").write_text("@a\nACGT\n+\n####\n", encoding="ascii")
        (self.root / "r.fasta").write_text(">r1\n" + "ACGT" * 8 + "\n", encoding="ascii")

    def config(self, plate_map: dict[str, str] | None) -> object:
        body = {
            "schema_version": 1,
            "run_name": "t",
            "output_root": str(self.root / "runs"),
            "inputs": [{"path": str(self.root / "reads.fastq"), "sample_id": "s"}],
            "reference_libraries": {"only": {"fasta": str(self.root / "r.fasta")}},
            "library": {
                "name": "x",
                "forward_motif": "ACGTACGTACGTACGTACGT",
                "reverse_motif": "TGCATGCATGCATGCATGCA",
            },
        }
        if plate_map is not None:
            body["plate_reference_map"] = plate_map
        path = self.root / "c.yaml"
        path.write_text(yaml.safe_dump(body), encoding="utf-8")
        return load_config(path)

    def test_a_barcode_absent_from_an_explicit_map_is_not_analysed(self) -> None:
        config = self.config({"RP04": "only"})
        self.assertEqual(config.reference_library_id_for_plate("RP04"), "only")
        with self.assertRaises(KeyError) as caught:
            config.reference_library_id_for_plate("RP08")
        self.assertIn("not part of this analysis", str(caught.exception))

    def test_no_map_still_falls_back_to_a_sole_library(self) -> None:
        # A simple one-library run should not have to write a map at all.
        config = self.config(None)
        self.assertEqual(config.reference_library_id_for_plate("anything"), "only")

    def test_the_fallback_does_not_depend_on_library_count(self) -> None:
        # The old rule keyed off "exactly one library", which is why the defect
        # appeared only in single-library configurations.
        config = self.config({"RP04": "only"})
        self.assertEqual(len(config.reference_sets), 1)
        with self.assertRaises(KeyError):
            config.reference_library_id_for_plate("RP99")
