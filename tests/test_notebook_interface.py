"""Execute the shipped notebook's local workflow without installation or Drive."""

from __future__ import annotations

import ast
import gzip
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
import zipfile
from contextlib import nullcontext
from pathlib import Path
from types import ModuleType
from unittest.mock import Mock, patch

NOTEBOOK = Path(__file__).resolve().parents[1] / "notebooks/Nanopore3_Colab.ipynb"


class NotebookInterfaceTests(unittest.TestCase):
    def test_installer_uses_the_local_checkout_and_active_interpreter(self):
        notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
        source = "".join(notebook["cells"][2]["source"])
        repo = NOTEBOOK.parents[1]
        with (
            patch("pathlib.Path.cwd", return_value=repo / "notebooks"),
            patch("subprocess.run") as run,
        ):
            exec(source, {})
        self.assertEqual(
            run.call_args.args[0],
            [sys.executable, "-m", "pip", "install", "-e", f"{repo}[report,mcp]"],
        )
        self.assertTrue(run.call_args.kwargs["check"])
        with patch("subprocess.run") as run:
            exec(
                source.replace("PACKAGE_SPEC = None", "PACKAGE_SPEC = 'custom.whl[report,mcp]'"), {}
            )
        self.assertEqual(run.call_args.args[0][-1], "custom.whl[report,mcp]")

    def test_drive_mount_cell_is_independent_of_setup(self):
        notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
        google = ModuleType("google")
        colab = ModuleType("google.colab")
        google.colab = colab
        colab.drive = Mock()
        with patch.dict(sys.modules, {"google": google, "google.colab": colab}):
            exec("".join(notebook["cells"][12]["source"]), {})
        colab.drive.mount.assert_called_once_with("/content/drive")

    def test_every_code_cell_is_valid_python_including_top_level_await(self):
        notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
        for index, cell in enumerate(notebook["cells"]):
            if cell["cell_type"] == "code":
                compile(
                    "".join(cell["source"]),
                    f"cell-{index}",
                    "exec",
                    flags=ast.PyCF_ALLOW_TOP_LEVEL_AWAIT,
                )
                self.assertEqual(cell["outputs"], [])
                self.assertIsNone(cell["execution_count"])

    def test_local_notebook_runs_fresh_examples_and_propagates_cli_failure(self):
        notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
        with tempfile.TemporaryDirectory(prefix="nanopore notebook ") as directory:
            root = Path(directory).resolve()
            # pandas is used only for display cells; the core-only CI environment
            # still executes setup, doctor, CLI example and input staging cells.
            optional = (
                {} if importlib.util.find_spec("pandas") else {"pandas": ModuleType("pandas")}
            )
            with (
                patch.dict(os.environ, {"NANOPORE3_WORKSPACE": str(root)}),
                patch.dict(sys.modules, optional) if optional else nullcontext(),
            ):
                namespace = {}

                def execute(index, **replacements):
                    source = "".join(notebook["cells"][index]["source"])
                    for old, new in replacements.items():
                        source = source.replace(old, new)
                    exec(compile(source, f"cell-{index}", "exec"), namespace)

                execute(4)
                execute(5)
                runtime = namespace["runtime"]
                self.assertEqual(runtime["parallel_defaults"]["jobs"], 0)
                self.assertGreater(runtime["group_memory_limit_bytes"], 0)
                self.assertLessEqual(runtime["available_cpus"], runtime["cpu_count"])
                execute(8)
                first_run = namespace["demo_run"]
                self.assertTrue((first_run / "stages/06_report/_SUCCESS").exists())
                execute(8)
                self.assertNotEqual(first_run, namespace["demo_run"])
                self.assertTrue((first_run / "stages/06_report/_SUCCESS").exists())
                fixtures = namespace["example"] / "fixtures"
                payload = (fixtures / "reads.fastq").read_bytes()
                source = fixtures / "my_run.fastq.gz"
                source.write_bytes(gzip.compress(payload))
                os.environ["SEQ_DATA"] = str(fixtures)
                execute(16, **{'"my_run.fastq"': '"my_run.fastq.gz"'})
                subset = namespace["local"]
                self.assertEqual(gzip.decompress(subset.read_bytes()), payload)
                execute(
                    16,
                    **{
                        '"my_run.fastq"': '"my_run.fastq.gz"',
                        "SUBSAMPLE = 100_000": "SUBSAMPLE = None",
                    },
                )
                self.assertNotEqual(namespace["local"], subset)
                self.assertEqual(namespace["local"].read_bytes(), source.read_bytes())
                with self.assertRaises(subprocess.CalledProcessError):
                    namespace["np3"](
                        "subsample", "--input", source, "--output", source, "--reads", 1
                    )
                self.assertEqual(gzip.decompress(source.read_bytes()), payload)
                execute(18)
                config = namespace["body"]
                self.assertEqual(config["parallel"]["jobs"], 0)
                self.assertEqual(config["parallel"]["backend"], "process")
                self.assertEqual(config["parallel"]["group_memory_mb"], 512)
                self.assertEqual(config["inputs"][0]["path"], str(namespace["local"]))
                self.assertEqual(config["output_root"], str(root / "runs"))
                namespace["run"] = namespace["demo_run"]
                namespace["config"] = namespace["example"] / "configs/example.yaml"
                execute(
                    25,
                    **{"Path.home()": "WORKSPACE", "SAVE_FULL_RUN = False": "SAVE_FULL_RUN = True"},
                )
                archive = namespace["OUT"]
                self.assertTrue((archive / "report/report.html").is_file())
                self.assertEqual(
                    (archive / "run.json").read_bytes(),
                    (namespace["run"] / "run.json").read_bytes(),
                )
                with zipfile.ZipFile(archive / "full-run.zip") as handle:
                    self.assertIn(
                        f"{namespace['run'].name}/stages/04_consensus/manifest.json",
                        handle.namelist(),
                    )


if __name__ == "__main__":
    unittest.main()
