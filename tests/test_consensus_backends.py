from __future__ import annotations

import subprocess
import unittest
from pathlib import Path
from unittest.mock import patch

from nanopore3.config import ConfigError, ConsensusSettings
from nanopore3.consensus import (
    ConsensusBackendError,
    ConsensusBackendUnavailable,
    ConsensusRead,
    build_reference_consensus,
)


class ConsensusBackendTests(unittest.TestCase):
    def test_backend_names_are_strictly_validated(self) -> None:
        self.assertEqual(ConsensusSettings().backend, "portable")
        self.assertEqual(ConsensusSettings(backend="mafft_spoa").backend, "mafft_spoa")
        with self.assertRaisesRegex(ConfigError, "portable or mafft_spoa"):
            ConsensusSettings(backend="spoa")

    @patch("nanopore3.consensus.shutil.which")
    def test_mafft_spoa_requires_both_optional_tools(self, which) -> None:
        which.side_effect = lambda name: "/tools/mafft" if name == "mafft" else None
        reads = [ConsensusRead(f"read-{index}", "ACGT") for index in range(3)]
        with self.assertRaisesRegex(
            ConsensusBackendUnavailable, "requires both MAFFT and SPOA.*missing: spoa"
        ):
            build_reference_consensus(
                "ACGT", reads, group_id="group", backend="mafft_spoa"
            )

    @patch("nanopore3.consensus.subprocess.run")
    @patch("nanopore3.consensus.shutil.which")
    def test_mafft_spoa_uses_safe_commands_and_stable_input(self, which, run) -> None:
        which.side_effect = lambda name: f"/tools/{name}"
        observed_inputs: list[str] = []

        def fake_run(command, **kwargs):
            self.assertIsInstance(command, list)
            self.assertFalse(kwargs["shell"])
            if Path(command[0]).name == "mafft":
                observed_inputs.append(Path(command[-1]).read_text(encoding="ascii"))
                return subprocess.CompletedProcess(
                    command, 0, stdout=">read_000001\nAC-GT\n", stderr=""
                )
            self.assertEqual(command[1:3], ["--algorithm", "msa"])
            self.assertEqual(
                Path(command[-1]).read_text(encoding="ascii"),
                ">read_000001\nAC-GT\n",
            )
            return subprocess.CompletedProcess(
                command, 0, stdout=">Consensus\nAC-GT\n", stderr=""
            )

        run.side_effect = fake_run
        reads = [ConsensusRead(f"read-{index}", "ACGT") for index in range(5)]
        first = build_reference_consensus(
            "ACGT",
            reads,
            group_id="group",
            max_reads=3,
            backend="mafft_spoa",
            threads=2,
        )
        second = build_reference_consensus(
            "ACGT",
            list(reversed(reads)),
            group_id="group",
            max_reads=3,
            backend="mafft_spoa",
            threads=2,
        )

        self.assertEqual(first.sequence, second.sequence)
        self.assertEqual(first.contributor_ids, second.contributor_ids)
        self.assertEqual(first.status, second.status)
        self.assertEqual(first.sequence, "ACGT")
        self.assertEqual(first.backend, "mafft_spoa")
        self.assertEqual(observed_inputs[0], observed_inputs[1])
        mafft_command = run.call_args_list[0].args[0]
        self.assertEqual(mafft_command[1:4], ["--quiet", "--thread", "2"])
        unlimited = build_reference_consensus(
            "ACGT", reads, group_id="group", max_reads=0, backend="mafft_spoa", threads=2
        )
        self.assertEqual(unlimited.n_reads_used, len(reads))
        self.assertEqual(set(unlimited.contributor_ids), {read.read_uid for read in reads})
        self.assertEqual(observed_inputs[-1].count(">"), len(reads))

    @patch("nanopore3.consensus.subprocess.run")
    @patch("nanopore3.consensus.shutil.which", return_value="/tools/tool")
    def test_external_tool_error_includes_actionable_stderr(self, _which, run) -> None:
        run.return_value = subprocess.CompletedProcess(
            ["mafft"], 2, stdout="", stderr="bad input FASTA"
        )
        reads = [ConsensusRead(f"read-{index}", "ACGT") for index in range(3)]
        with self.assertRaisesRegex(
            ConsensusBackendError, "MAFFT failed with exit code 2: bad input FASTA"
        ):
            build_reference_consensus(
                "ACGT", reads, group_id="group", backend="mafft_spoa"
            )


if __name__ == "__main__":
    unittest.main()
