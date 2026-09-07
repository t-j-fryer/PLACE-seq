"""Exercise the optional server over real MCP stdio and HTTP sessions."""

from __future__ import annotations

import asyncio
import gzip
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch

HAS_MCP = all(importlib.util.find_spec(name) is not None for name in ("mcp", "psutil"))


@unittest.skipUnless(HAS_MCP, "install nanopore3[mcp] to test the optional server")
class McpSessionTests(unittest.IsolatedAsyncioTestCase):
    async def asyncSetUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="nanopore3 mcp ")
        self.root = Path(self.temp.name).resolve()

    async def asyncTearDown(self):
        self.temp.cleanup()

    async def call(self, session, name, arguments=None):
        result = await session.call_tool(name, arguments or {})
        self.assertFalse(result.isError, result.content)
        return result.structuredContent

    async def completed(self, session, job):
        for _ in range(300):
            status = await self.call(session, "job_status", {"job_id": job["job_id"]})
            if status["state"] in {"succeeded", "failed", "cancelled"}:
                return status
            await asyncio.sleep(0.05)
        self.fail("job did not finish in 15 seconds")

    async def test_stdio_complete_workflow_errors_resume_and_rerun(self):
        from mcp import ClientSession, StdioServerParameters
        from mcp.client.stdio import stdio_client

        params = StdioServerParameters(
            command=sys.executable,
            args=["-m", "nanopore3.mcp_server", "--workspace", str(self.root)],
            env=dict(os.environ),
        )
        async with stdio_client(params) as (read, write), ClientSession(read, write) as session:
            await session.initialize()
            listed = await session.list_tools()
            self.assertTrue({"start_run", "job_status", "read_artifact"}.issubset(
                {tool.name for tool in listed.tools}
            ))
            self.assertTrue(all(t.outputSchema for t in listed.tools))
            resource = await session.read_resource("nanopore3://guide")
            self.assertIn("job_status", resource.contents[0].text)
            info = await self.call(session, "workspace_info")
            self.assertEqual(info["runtime"]["parallel_defaults"]["jobs"], 0)
            self.assertGreater(info["runtime"]["group_memory_limit_bytes"], 0)
            self.assertIn("sha256", info["runtime"]["analysis_implementation"])
            prompt = await session.get_prompt("analyse_run", {"config_path": "example.yaml"})
            self.assertTrue(prompt.messages)
            job = await self.call(session, "init_example", {"directory": "example project"})
            self.assertEqual((await self.completed(session, job))["state"], "succeeded")
            config = "example project/configs/example.yaml"
            job = await self.call(session, "start_subsample", {
                "input_path": "example project/fixtures/reads.fastq",
                "output_path": "subset.fastq.gz", "reads": 2,
            })
            self.assertEqual((await self.completed(session, job))["state"], "succeeded")
            before = (self.root / "subset.fastq.gz").read_bytes()
            self.assertEqual(len(gzip.decompress(before).splitlines()), 8)
            for output in ("subset.fastq.gz", "example project/fixtures/reads.fastq", "../escape"):
                result = await session.call_tool("start_subsample", {
                    "input_path": "example project/fixtures/reads.fastq",
                    "output_path": output, "reads": 2,
                })
                self.assertTrue(result.isError)
            self.assertEqual((self.root / "subset.fastq.gz").read_bytes(), before)
            job = await self.call(session, "start_validation", {"config_path": config})
            validated = await self.completed(session, job)
            self.assertEqual(validated["state"], "succeeded", validated)
            self.assertIn("references", validated["result"])
            job = await self.call(session, "start_run", {"config_path": config, "run_id": "demo"})
            status = await self.completed(session, job)
            self.assertEqual(status["state"], "succeeded", status)
            summary = await self.call(session, "run_summary", {"run_dir": "runs/demo"})
            self.assertEqual(len(summary["stages"]), 6)
            self.assertIn("qc", summary["outputs"])
            self.assertFalse(summary["checksums_verified"])
            self.assertTrue(summary["implementation_matches_current"])
            self.assertEqual(
                summary["analysis_implementation"], info["runtime"]["analysis_implementation"]
            )
            preview = await self.call(
                session,
                "read_artifact",
                {
                    "path": summary["outputs"]["qc"],
                    "max_chars": 30,
                },
            )
            self.assertEqual(len(preview["text"]), 30)
            self.assertTrue(preview["truncated"])
            job = await self.call(session, "start_run", {
                "config_path": config, "run_id": "demo", "resume": True,
            })
            self.assertEqual((await self.completed(session, job))["state"], "succeeded")
            metadata_path = self.root / "runs/demo/run.json"
            original_metadata = metadata_path.read_bytes()
            metadata = json.loads(original_metadata)
            metadata.pop("analysis_implementation")
            metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
            legacy = await self.call(session, "run_summary", {"run_dir": "runs/demo"})
            self.assertFalse(legacy["implementation_matches_current"])
            job = await self.call(session, "start_run", {
                "config_path": config, "run_id": "demo", "resume": True,
            })
            refused = await self.completed(session, job)
            self.assertEqual(refused["state"], "failed")
            self.assertIn("missing implementation identity", refused["stderr_tail"])
            metadata_path.write_bytes(original_metadata)
            job = await self.call(session, "start_rerun", {
                "config_path": config, "from_run": "runs/demo", "run_id": "rerun",
            })
            self.assertEqual((await self.completed(session, job))["state"], "succeeded")
            await self.call(session, "save_config", {"path": "invalid.yaml", "yaml_text": "x: 1"})
            job = await self.call(session, "start_validation", {"config_path": "invalid.yaml"})
            failed = await self.completed(session, job)
            self.assertEqual(failed["state"], "failed")
            self.assertIn("error:", failed["stderr_tail"])
            for name, args in (
                ("read_artifact", {"path": "../outside"}),
                ("read_artifact", {"path": ".git/config"}),
                ("start_run", {"config_path": config, "run_id": "../escaped"}),
                ("start_run", {"config_path": config, "run_id": "CON"}),
                ("start_run", {"config_path": config, "resume": True}),
                ("start_run", {"config_path": config, "run_id": "demo"}),
                ("start_run", {"config_path": config, "output_root": "../outside"}),
                ("save_config", {"path": "invalid.yaml", "yaml_text": "x: 2"}),
                ("read_artifact", {"path": config, "max_chars": 1_000_000}),
            ):
                result = await session.call_tool(name, args)
                self.assertTrue(result.isError, (name, args))
            history = await self.call(session, "list_jobs")
            self.assertGreaterEqual(len(history["jobs"]), 6)

        # A new server can inspect completed jobs without pretending to own a PID.
        async with stdio_client(params) as (read, write), ClientSession(read, write) as session:
            await session.initialize()
            stored = await self.call(session, "job_status", {"job_id": job["job_id"]})
            self.assertEqual(stored["state"], "failed")

    async def test_streamable_http_protocol(self):
        import httpx
        from mcp import ClientSession
        from mcp.client.streamable_http import streamable_http_client

        from nanopore3.mcp_server import create_server

        server = create_server(self.root)
        app = server.streamable_http_app()

        async with app.router.lifespan_context(app), httpx.AsyncClient(
            transport=httpx.ASGITransport(app=app),
        ) as http, streamable_http_client(
            "http://127.0.0.1:8000/mcp", http_client=http,
        ) as (read, write, _), ClientSession(read, write) as session:
            await session.initialize()
            result = await self.call(session, "list_files")
            self.assertEqual(result["entries"], [])
            job = await self.call(session, "init_example", {"directory": "http-example"})
            self.assertEqual((await self.completed(session, job))["state"], "succeeded")
            job = await self.call(session, "start_validation", {
                "config_path": "http-example/configs/example.yaml",
            })
            self.assertEqual((await self.completed(session, job))["state"], "succeeded")
            job = await self.call(session, "start_run", {
                "config_path": "http-example/configs/example.yaml",
            })
        stored = json.loads((self.root / ".nanopore3-mcp/jobs" / job["job_id"] / "job.json")
                            .read_text(encoding="utf-8"))
        self.assertIn(stored["state"], {"cancelled", "succeeded"}, stored)

    async def test_stdio_disconnect_finishes_job_lifecycle(self):
        from mcp import ClientSession, StdioServerParameters
        from mcp.client.stdio import stdio_client

        params = StdioServerParameters(
            command=sys.executable,
            args=["-m", "nanopore3.mcp_server", "--workspace", str(self.root)],
            env=dict(os.environ),
        )
        async with stdio_client(params) as (read, write), ClientSession(read, write) as session:
            await session.initialize()
            job = await self.call(session, "init_example", {"directory": "example"})
            await self.completed(session, job)
            job = await self.call(session, "start_run", {
                "config_path": "example/configs/example.yaml",
            })
        path = self.root / ".nanopore3-mcp/jobs" / job["job_id"] / "job.json"
        stored = json.loads(path.read_text(encoding="utf-8"))
        self.assertIn(stored["state"], {"cancelled", "succeeded"}, stored)

    async def test_bounded_gzip_and_symlink_escape(self):
        from nanopore3.mcp_server import create_server

        server = create_server(self.root)
        with gzip.open(self.root / "large.csv.gz", "wt") as handle:
            handle.write("a" * 2_000_000)
        _, result = await server.call_tool("read_artifact", {
            "path": "large.csv.gz", "offset": 100, "max_chars": 200,
        })
        self.assertEqual(result["text"], "a" * 200)
        self.assertEqual(result["next_offset"], 300)
        with tempfile.TemporaryDirectory() as outside:
            try:
                (self.root / "escape").symlink_to(outside, target_is_directory=True)
            except OSError:
                self.skipTest("symlinks need extra privileges on this Windows installation")
            with self.assertRaises(Exception) as caught:
                await server.call_tool("save_config", {
                    "path": "escape/stolen.yaml", "yaml_text": "x: 1",
                })
            self.assertIn("workspace", str(caught.exception))
            self.assertFalse((Path(outside) / "stolen.yaml").exists())


@unittest.skipUnless(HAS_MCP, "install nanopore3[mcp] to test job lifecycle")
class JobLifecycleTests(unittest.TestCase):
    def test_duplicate_destination_is_reserved_with_multiple_slots(self):
        from nanopore3.mcp_jobs import JobManager, Workspace

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            manager = JobManager(Workspace(root), max_jobs=2)
            popen = subprocess.Popen

            def launch(_args, **kwargs):
                return popen([sys.executable, "-c", "import time; time.sleep(60)"], **kwargs)

            try:
                with patch("nanopore3.mcp_jobs.subprocess.Popen", side_effect=launch):
                    manager.start(["run"], run_dir=root / "runs/demo")
                with self.assertRaisesRegex(ValueError, "already owns"):
                    manager.start(["run"], run_dir=root / "runs/demo")
            finally:
                manager.close()

    def test_capacity_cancellation_descendants_shutdown_and_restart(self):
        import psutil

        from nanopore3.mcp_jobs import JobManager, Workspace

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            manager = JobManager(Workspace(root))
            popen = subprocess.Popen
            code = (
                "import subprocess,sys,time; "
                "child=subprocess.Popen([sys.executable,'-c','import time; time.sleep(60)']); "
                "print(child.pid,flush=True); time.sleep(60)"
            )

            def launch(_args, **kwargs):
                return popen([sys.executable, "-c", code], **kwargs)

            try:
                with patch("nanopore3.mcp_jobs.subprocess.Popen", side_effect=launch):
                    job = manager.start(["run"])
                for _ in range(100):
                    status = manager.status(job["job_id"])
                    if status["stdout_tail"].strip():
                        break
                    time.sleep(0.02)
                child_pid = int(status["stdout_tail"].strip())
                with self.assertRaisesRegex(ValueError, "capacity"):
                    manager.start(["doctor"])
                other = JobManager(Workspace(root))
                self.assertEqual(other.status(job["job_id"])["state"], "unknown")
                other.close()
                manager.close()
                for _ in range(100):
                    status = manager.status(job["job_id"])
                    if status["state"] == "cancelled":
                        break
                    time.sleep(0.02)
                self.assertEqual(status["state"], "cancelled", status)
                try:
                    child = psutil.Process(child_pid)
                    self.assertEqual(child.status(), psutil.STATUS_ZOMBIE)
                except psutil.NoSuchProcess:
                    pass
                with self.assertRaisesRegex(ValueError, "shutting down"):
                    manager.start(["doctor"])
                saved = json.loads(next((root / ".nanopore3-mcp/jobs").glob("*/job.json"))
                                   .read_text(encoding="utf-8"))
                self.assertIsNotNone(saved["completed_utc"])
            finally:
                manager.close()


if __name__ == "__main__":
    unittest.main()
