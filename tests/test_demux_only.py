"""Demux-only exports preserve reads, stop early and need no reference files."""

from __future__ import annotations

import csv
import gzip
import json
import subprocess
import sys
from dataclasses import replace
from unittest.mock import patch

import pytest
import yaml

from nanopore3.config import ConfigError, load_config
from nanopore3.demux_export import DemuxFastqWriter
from nanopore3.pipeline import run_pipeline, validate_inputs
from nanopore3.provenance import StageValidationError, validate_stage_directory


@pytest.fixture
def config_file(tmp_path):
    reads = [
        ("forward", "AAAACCCC" + "ACGTACGT" + "GATTACAGG"),
        ("unknown_well", "AAAACCCC" + "TATATATA" + "GATTACAGG"),
        ("unknown_plate", "TATATATA" + "ACGTACGT" + "GATTACAGG"),
        ("too_short", "AAAACCCC"),
    ]
    from nanopore3.demux import reverse_complement

    sequence = reads[0][1]
    reads.append(("reverse", reverse_complement(sequence)))
    quality = "".join(chr(40 + i) for i in range(len(sequence)))
    with (tmp_path / "reads.fastq").open("w") as handle:
        for name, seq in reads:
            qual = quality[::-1] if name == "reverse" else "I" * len(seq)
            handle.write(f"@{name}\n{seq}\n+\n{qual}\n")
    body = dict(
        schema_version=1,
        workflow="demux_only",
        output_root="runs",
        inputs=[dict(path="reads.fastq", sample_id="sample")],
        library=dict(minimum_read_length=20),
        barcodes=dict(
            plate=dict(
                sequences=dict(P1="AAAACCCC"),
                max_edits=0,
                search_window=30,
                search_ends=["head"],
                allow_reverse_complement=True,
            ),
            well=dict(
                sequences=dict(A1="ACGTACGT"),
                max_edits=0,
                search_window=30,
                allow_reverse_complement=False,
            ),
        ),
        parallel=dict(jobs=1, chunk_reads=1),
    )
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(body))
    return path


def read_fastq(path):
    with gzip.open(path, "rt") as handle:
        lines = handle.read().splitlines()
    return [(lines[i], lines[i + 1], lines[i + 3]) for i in range(0, len(lines), 4)]


def test_no_references_or_consensus_backend_needed(config_file):
    config = replace(
        load_config(config_file),
        consensus=replace(load_config(config_file).consensus, backend="mafft_spoa"),
    )
    with (
        patch("nanopore3.pipeline.read_reference_libraries", side_effect=AssertionError),
        patch("nanopore3.pipeline.require_mafft_spoa", side_effect=AssertionError),
    ):
        assert validate_inputs(config)["references"] == 0
        run = run_pipeline(config, run_id="demux")
    assert {p.name for p in (run / "stages").iterdir()} == {"01_ingest", "02_demux"}
    assert not (run / "consensus_by_plate").exists()
    stage = run / "stages/02_demux"
    validate_stage_directory(stage)
    index = list(csv.DictReader((stage / "reads/index.csv").open()))
    assert {r["scope"]: int(r["reads"]) for r in index} == {
        "plate": 3,
        "well": 2,
        "unresolved_well": 1,
    }
    plate = read_fastq(stage / "reads/by_plate/P1.fastq.gz")
    well = read_fastq(stage / "reads/by_well/P1/A1.fastq.gz")
    assert len(plate) == 3 and len(well) == 2
    assert [r[1] for r in well] == ["AAAACCCCACGTACGTGATTACAGG"] * 2
    assert well[1][2] == "".join(chr(40 + i) for i in range(25))
    assert "read_uid=" in well[1][0] and "reverse" in well[1][0]
    with gzip.open(stage / "demuxed_reads.jsonl.gz", "rt") as handle:
        accepted = [json.loads(line) for line in handle]
    assert len(accepted) == 2 and all(r["well_id"] == "A1" for r in accepted)
    original = (stage / "manifest.json").read_bytes()
    run_pipeline(config, run_id="demux", resume=True)
    assert (stage / "manifest.json").read_bytes() == original
    with (stage / "reads/by_well/P1/A1.fastq.gz").open("ab") as handle:
        handle.write(b"corrupt")
    with pytest.raises(StageValidationError):
        run_pipeline(config, run_id="demux", resume=True)


def test_full_requires_references_but_cli_override_does_not(config_file):
    body = yaml.safe_load(config_file.read_text())
    body.pop("workflow")
    config_file.write_text(yaml.safe_dump(body))
    with pytest.raises(ConfigError):
        load_config(config_file)
    assert load_config(config_file, demux_only=True).reference_sets == {}
    command = [sys.executable, "-m", "nanopore3"]
    validated = subprocess.run(
        command + ["validate", "--config", str(config_file), "--demux-only"],
        capture_output=True,
        text=True,
    )
    assert validated.returncode == 0, validated.stderr
    completed = subprocess.run(
        command + ["run", "--config", str(config_file), "--demux-only", "--run-id", "cli"],
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
    assert "02_demux/reads/index.csv" in completed.stdout
    assert "06_report" not in completed.stdout


def test_missing_reference_file_is_ignored_only_in_demux(config_file):
    body = yaml.safe_load(config_file.read_text())
    body["references"] = {"fasta": "missing.fasta"}
    config_file.write_text(yaml.safe_dump(body))
    config = load_config(config_file)
    assert validate_inputs(config)["references"] == 0
    with pytest.raises(Exception, match="missing input"):
        validate_inputs(replace(config, workflow="full"))


@pytest.mark.parametrize("value", ["bad", True, [], None])
def test_invalid_workflow_rejected(config_file, value):
    body = yaml.safe_load(config_file.read_text())
    body["workflow"] = value
    config_file.write_text(yaml.safe_dump(body))
    with pytest.raises(ConfigError):
        load_config(config_file)


def test_empty_export_is_valid(config_file):
    config = load_config(config_file)
    config = replace(config, library=replace(config.library, minimum_read_length=1000))
    run = run_pipeline(config, run_id="empty")
    root = run / "stages/02_demux/reads"
    assert list(csv.DictReader((root / "index.csv").open())) == []
    assert not list(root.rglob("*.fastq.gz"))
    validate_stage_directory(run / "stages/02_demux")


def test_lru_gzip_members_and_portable_paths(tmp_path):
    with DemuxFastqWriter(tmp_path / "reads", max_open=1) as writer:
        for i in range(5):
            for plate in ["P1", "p1", "../escape", "CON"]:
                writer.write(
                    dict(
                        plate_id=plate,
                        well_id="A1",
                        original_read_id=f"r{i}",
                        read_uid=f"u{i}",
                        sample_id="s",
                        sequence="ACGT",
                        quality="IJKL",
                    ),
                    assigned=True,
                )
                assert len(writer.handles) <= 1
    index = list(csv.DictReader((tmp_path / "reads/index.csv").open()))
    assert len(index) == 8
    assert len({r["file"].casefold() for r in index}) == 8
    for row in index:
        path = tmp_path / "reads" / row["file"]
        assert path.resolve().is_relative_to((tmp_path / "reads").resolve())
        assert len(read_fastq(path)) == int(row["reads"]) == 5


def test_notebook_demux_results_and_archive(config_file):
    import shutil
    from pathlib import Path
    from uuid import uuid4

    pd = pytest.importorskip("pandas")
    notebook = json.loads(
        (Path(__file__).parents[1] / "notebooks/Nanopore3_Colab.ipynb").read_text()
    )
    run = run_pipeline(load_config(config_file), run_id="notebook")
    namespace = dict(
        run=run,
        config=config_file,
        pd=pd,
        display=lambda value: None,
        WORKSPACE=config_file.parent,
        IN_COLAB=False,
        Path=Path,
        shutil=shutil,
        uuid4=uuid4,
    )
    exec("".join(notebook["cells"][23]["source"]), namespace)
    assert len(namespace["reads"]) == 3
    archive = "".join(notebook["cells"][25]["source"]).replace("Path.home()", "WORKSPACE")
    exec(archive, namespace)
    assert (namespace["OUT"] / "demux/reads/by_well/P1/A1.fastq.gz").is_file()
    assert not (namespace["OUT"] / "qc.csv.gz").exists()


def test_plate_only_panel_has_no_placeholder_wells(config_file):
    body = yaml.safe_load(config_file.read_text())
    body["barcodes"].pop("well")
    config_file.write_text(yaml.safe_dump(body))
    run = run_pipeline(load_config(config_file), run_id="plates")
    index = list(csv.DictReader((run / "stages/02_demux/reads/index.csv").open()))
    assert len(index) == 1 and index[0]["scope"] == "plate"
    assert int(index[0]["reads"]) == 3
