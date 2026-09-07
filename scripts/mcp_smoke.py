#!/usr/bin/env python3
"""Connect to PLACE-seq via MCP; optionally run a new packaged synthetic example."""

from __future__ import annotations

import argparse
import asyncio
import json
import os
import sys
from pathlib import Path
from uuid import uuid4

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
from mcp.client.streamable_http import streamable_http_client


async def call(session: ClientSession, name: str, arguments: dict) -> dict:
    result = await session.call_tool(name, arguments)
    if result.isError:
        raise RuntimeError(str(result.content))
    return result.structuredContent or json.loads(result.content[0].text)


async def finish(session: ClientSession, job: dict) -> dict:
    for _ in range(120):
        result = await call(session, "job_status", {"job_id": job["job_id"]})
        if result["state"] == "succeeded":
            return result
        if result["state"] in {"failed", "cancelled", "unknown"}:
            raise RuntimeError(json.dumps(result, indent=2))
        await asyncio.sleep(0.5)
    raise TimeoutError(f"Job {job['job_id']} still running; inspect job_status")


async def exercise(read, write, example: bool) -> None:
    async with ClientSession(read, write) as session:
        await session.initialize()
        print("Tools:", ", ".join(tool.name for tool in (await session.list_tools()).tools))
        await session.read_resource("nanopore3://guide")
        info = await call(session, "workspace_info", {})
        print("Allocated CPUs:", info["runtime"]["available_cpus"])
        print("Group estimate budget (bytes):", info["runtime"]["group_memory_limit_bytes"])
        print("Source identity:", info["runtime"]["analysis_implementation"]["sha256"])
        if example:
            name = "smoke-" + uuid4().hex[:12]
            await finish(session, await call(session, "init_example", {"directory": name}))
            await finish(session, await call(session, "start_subsample", {
                "input_path": f"{name}/fixtures/reads.fastq",
                "output_path": f"{name}/subset.fastq.gz", "reads": 2,
            }))
            config = f"{name}/configs/example.yaml"
            await finish(session, await call(session, "start_validation", {"config_path": config}))
            job = await finish(session, await call(session, "start_run", {"config_path": config}))
            summary = await call(session, "run_summary", {"run_dir": job["run_dir"]})
            if summary["implementation_matches_current"] is not True:
                raise RuntimeError("Fresh run source identity does not match the server")
            print(json.dumps(summary, indent=2))


async def run(args) -> None:
    if args.url:
        async with streamable_http_client(args.url) as (read, write, _):
            await exercise(read, write, args.example)
    else:
        params = StdioServerParameters(
            command=sys.executable,
            args=["-m", "nanopore3.mcp_server", "--workspace", str(args.workspace.resolve())],
            env=dict(os.environ),
        )
        async with stdio_client(params) as (read, write):
            await exercise(read, write, args.example)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    target = parser.add_mutually_exclusive_group(required=True)
    target.add_argument("--workspace", type=Path, help="existing directory; launch a stdio server")
    target.add_argument("--url", help="existing Streamable HTTP endpoint, including /mcp")
    parser.add_argument("--example", action="store_true", help="create and run synthetic inputs")
    asyncio.run(run(parser.parse_args()))


if __name__ == "__main__":
    main()
