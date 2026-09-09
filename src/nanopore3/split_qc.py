"""Independent grades for an explicitly bounded insert and sequenced vector."""

import re

import edlib

from .qc import translate

SPLIT_FIELDS = (
    "insert_grade",
    "vector_status",
    "vector_edit_distance",
    "insert_ambiguous_bases",
    "vector_ambiguous_bases",
    "insert_query_length",
    "insert_coding_status",
)


def insert_spans(reference, left, right):
    """Boundaries are unique reference motifs, excluded from the insert."""
    reference, left, right = reference.upper(), left.upper(), right.upper()
    if reference.count(left) != 1 or reference.count(right) != 1:
        raise ValueError("Insert boundary motifs must each occur exactly once in every reference")
    start, end = reference.index(left) + len(left), reference.index(right)
    if end <= start:
        raise ValueError("Insert boundary motifs must enclose a nonempty forward insert")
    return {"flank_5p": (0, start), "insert": (start, end), "flank_3p": (end, len(reference))}


def split_grades(
    query,
    reference,
    spans,
    *,
    upstream=None,
    downstream=None,
    left=None,
    right=None,
    minimum_identity=0.98,
    minimum_query_coverage=0.95,
    minimum_reference_coverage=0.95,
    length_tolerance=10,
):
    """Project a global alignment onto the insert; outside edits cannot downgrade it.

    Insertions belong to the region at their next reference position (terminal
    insertions belong to the final region). Ambiguity always prevents exactness.
    Coding checks use the reference frame and reference bases for partial edge
    codons, so these describe the insert itself, not the complete expression ORF.
    The existing whole-ORF checks remain available separately.
    """
    query, reference = query.upper(), reference.upper()
    start, end = spans["insert"]
    edits = [0, 0]  # insert, vector
    ambiguities = [0, 0]
    insert = []
    qpos = rpos = 0
    segments = [(query, reference)]
    if left and right:
        left, right = left.upper(), right.upper()
        if query.count(left) == 1 and query.count(right) == 1:
            qstart, qend = query.index(left) + len(left), query.index(right)
            if qstart <= qend:
                # Pin intact boundary motifs: equally optimal global alignments
                # in repetitive sequence must not move an insert deletion into
                # the vector or vice versa.
                segments = [
                    (query[:qstart], reference[:start]),
                    (query[qstart:qend], reference[start:end]),
                    (query[qend:], reference[end:]),
                ]
    region_edits = {name: 0 for name in spans}
    steps = []
    for index, (q, r) in enumerate(segments):
        cigar = edlib.align(q, r, mode="NW", task="path")["cigar"] or (
            f"{len(r)}D" if r else f"{len(q)}I"
        )
        forced = ("flank_5p", "insert", "flank_3p")[index] if len(segments) == 3 else None
        steps.extend((int(count), op, forced) for count, op in re.findall(r"(\d+)([=XID])", cigar))
    for count, op, forced in steps:
        for _ in range(count):
            name = forced or next(
                name for name, (lo, hi) in spans.items() if lo <= min(rpos, len(reference) - 1) < hi
            )
            region = 0 if name == "insert" else 1
            if op != "=":
                edits[region] += 1
                region_edits[name] += 1
            if op != "D":
                base = query[qpos]
                ambiguities[region] += base not in "ACGT"
                if region == 0:
                    insert.append(base)
                qpos += 1
            if op != "I":
                rpos += 1
    sequence = "".join(insert)
    expected = end - start
    coding = "not_evaluable"
    frame_failure = stop_failure = False
    if upstream and downstream:
        upstream, downstream = upstream.upper(), downstream.upper()
        if reference.count(upstream) == 1 and reference.count(downstream) == 1:
            origin, stop = reference.index(upstream), reference.index(downstream)
            if origin <= start < end <= stop:
                frame_failure = (len(sequence) - expected) % 3 != 0
                phase = (start - origin) % 3
                suffix = (-(end - origin)) % 3
                peptide = translate(
                    reference[start - phase : start] + sequence + reference[end : end + suffix]
                )
                stop_failure = "*" in peptide
                coding = "fail" if frame_failure or stop_failure else "pass"
    if ambiguities[0]:
        grade = "mixed_variants"
        coding = "not_evaluable"
    elif frame_failure:
        grade = "frameshift"
    elif stop_failure:
        grade = "premature_stop"
    elif abs(len(sequence) - expected) > length_tolerance:
        grade = "truncated"
    elif (
        1 - edits[0] / max(expected, len(sequence)) < minimum_identity
        or min(1.0, expected / max(1, len(sequence))) < minimum_query_coverage
        or min(1.0, len(sequence) / expected) < minimum_reference_coverage
    ):
        grade = "mismatched"
    elif edits[0] == 0:
        grade = "perfect"
    else:
        grade = "screenable" if coding == "pass" else "not_evaluable"
    result = dict(
        zip(
            SPLIT_FIELDS,
            (
                grade,
                "perfect" if edits[1] == 0 and not ambiguities[1] else "edited",
                edits[1],
                ambiguities[0],
                ambiguities[1],
                len(sequence),
                coding,
            ),
            strict=True,
        )
    )
    for name, (lo, hi) in spans.items():
        result[f"{name}_length"] = hi - lo
        result[f"{name}_edit_distance"] = region_edits[name]
        result[f"{name}_identity"] = (
            f"{max(0.0, 1 - region_edits[name] / (hi - lo)):.6f}" if hi > lo else ""
        )
    return result
