"""Insert/vector separation must not promote damaged inserts or hide backbone edits."""

import csv
from pathlib import Path
from tempfile import TemporaryDirectory

import pytest

from nanopore3.config import ConfigError, QcSettings, _parse_qc
from nanopore3.export import write_consensus_tree
from nanopore3.split_qc import insert_spans, split_grades

LEFT = "CCGATGCAGCTT"
RIGHT = "TCGGGCCACTAA"
INSERT = "GACGTTGCCGAG" * 8
REFERENCE = "TTTCCCAA" + LEFT + INSERT + RIGHT + "GGGACCCC"
SPANS = insert_spans(REFERENCE, LEFT, RIGHT)
START, END = SPANS["insert"]


def check(query, **kwargs):
    result = split_grades(
        query,
        REFERENCE,
        SPANS,
        upstream="ATGCAGCTT",
        downstream=RIGHT,
        left=LEFT,
        right=RIGHT,
        **kwargs,
    )
    assert (
        result["vector_edit_distance"]
        == result["flank_5p_edit_distance"] + result["flank_3p_edit_distance"]
    )
    return result


@pytest.mark.parametrize("replacement", ["A", "N", "", "AC"])
def test_backbone_changes_leave_exact_insert_perfect(replacement):
    result = check(replacement + REFERENCE[1:])
    assert result["insert_grade"] == "perfect"
    assert result["vector_status"] == "edited"
    assert result["insert_coding_status"] == "pass"


def test_exact_both():
    result = check(REFERENCE)
    assert result["insert_grade"] == result["vector_status"] == "perfect"


@pytest.mark.parametrize(
    "codon, grade",
    [
        ("AAC", "screenable"),
        ("TAA", "premature_stop"),
        ("GNC", "mixed_variants"),
        ("G", "frameshift"),
    ],
)
def test_insert_changes_do_not_downgrade_vector(codon, grade):
    query = REFERENCE[: START + 12] + codon + REFERENCE[START + 15 :]
    result = check(query)
    assert result["insert_grade"] == grade
    assert result["vector_status"] == "perfect"


def test_final_insert_codon_stop_is_premature():
    assert check(REFERENCE[: END - 3] + "TAA" + REFERENCE[END:])["insert_grade"] == "premature_stop"


def test_large_inframe_deletion_is_truncated():
    assert check(REFERENCE[: START + 12] + REFERENCE[START + 36 :])["insert_grade"] == "truncated"


@pytest.mark.parametrize("position", [START, END])
def test_insertions_at_intact_insert_edges_belong_to_insert(position):
    result = check(REFERENCE[:position] + "GCC" + REFERENCE[position:], minimum_identity=0.95)
    assert result["insert_edit_distance"] == 3
    assert result["insert_grade"] == "screenable"
    assert result["vector_status"] == "perfect"


def test_unknown_coding_is_explicit():
    result = split_grades(REFERENCE, REFERENCE, SPANS)
    assert result["insert_grade"] == "perfect"  # DNA exactness, not protein functionality
    assert result["insert_coding_status"] == "not_evaluable"
    query = REFERENCE[:START] + "AAC" + REFERENCE[START + 3 :]
    assert split_grades(query, REFERENCE, SPANS)["insert_grade"] == "not_evaluable"


def test_boundaries_are_unique_ordered_and_excluded():
    assert REFERENCE[START:END] == INSERT
    for ref, left, right in [
        (REFERENCE, RIGHT, LEFT),
        (REFERENCE, "AAAAAA", RIGHT),
        (REFERENCE + LEFT, LEFT, RIGHT),
    ]:
        with pytest.raises(ValueError):
            insert_spans(ref, left, right)


def test_yaml_pair_validation():
    with pytest.raises(ConfigError):
        QcSettings(insert_left_boundary=LEFT)
    parsed = _parse_qc(dict(insert_left_boundary=LEFT.lower(), insert_right_boundary=RIGHT))
    assert parsed.insert_left_boundary == LEFT


def test_export_scoped_filename_and_summary():
    from test_export import consensus, qc

    with TemporaryDirectory() as temp:
        root = Path(temp) / "tree"
        summary = write_consensus_tree(
            [consensus()],
            {"cons-1": qc(**check("A" + REFERENCE[1:]))},
            {"cons-1": "A" + REFERENCE[1:]},
            root,
        )
        assert summary["insert_grades"] == {"perfect": 1}
        assert summary["vector_statuses"] == {"edited": 1}
        with (root / "index.csv").open() as handle:
            row = next(csv.DictReader(handle))
        assert row["file"].endswith("__insert-perfect__vector-edited.fasta")
        assert "insert_grade=perfect vector_status=edited" in (root / row["file"]).read_text()
