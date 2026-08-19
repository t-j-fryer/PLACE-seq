"""Dependency-free deterministic HTML reporting."""

from __future__ import annotations

import csv
from collections import Counter
from collections.abc import Iterable, Mapping
from html import escape
from pathlib import Path


def status_counts(rows: Iterable[Mapping[str, object]], field: str) -> Counter[str]:
    return Counter(str(row.get(field, "missing")) for row in rows)


def _format_value(value: object) -> str:
    """Render one cell without assuming the value is a count."""

    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float):
        return f"{value:,.4g}"
    return escape(str(value))


def _rows(values: Mapping[str, object], depth: int = 0) -> str:
    """Render a summary mapping, nesting one level for grouped sections.

    Sections are not all flat counts: a per-plate breakdown is a mapping of
    mappings, and some entries name a plate rather than count one.  Coercing
    every value to an integer made adding such a section a runtime failure in
    the last stage of a twenty-minute run.
    """

    out: list[str] = []
    for key, value in sorted(values.items()):
        label = escape(str(key))
        if isinstance(value, Mapping):
            out.append(
                f'<tr><td colspan="2" style="background:#f6f9fb">'
                f"<strong>{label}</strong></td></tr>"
            )
            out.append(_rows(value, depth + 1))
            continue
        pad = "padding-left:1.6rem" if depth else ""
        out.append(
            f'<tr><td style="{pad}">{label}</td><td>{_format_value(value)}</td></tr>'
        )
    return "".join(out)


def write_html_report(
    path: Path,
    *,
    title: str,
    sections: Mapping[str, Mapping[str, object]],
    provenance: Mapping[str, object],
) -> None:
    """Write a compact report from canonical summaries only."""

    blocks: list[str] = []
    for heading, values in sections.items():
        blocks.append(
            f"<h2>{escape(heading)}</h2>"
            f"<table><tr><th>Item</th><th>Value</th></tr>{_rows(values)}</table>"
        )
    provenance_rows = "".join(
        f"<tr><td>{escape(str(k))}</td><td>{escape(str(v))}</td></tr>"
        for k, v in sorted(provenance.items())
    )
    html = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>{escape(title)}</title>
<style>body{{font:16px system-ui;max-width:960px;margin:2rem auto;padding:0 1rem;color:#18202a}}
table{{border-collapse:collapse;margin-bottom:2rem}}
th,td{{border:1px solid #ccd4dd;padding:.4rem .7rem;text-align:left}}
th{{background:#edf3f7}}code{{background:#edf3f7;padding:.1rem .25rem}}</style></head>
<body><h1>{escape(title)}</h1>{''.join(blocks)}
<h2>Provenance</h2><table>{provenance_rows}</table></body></html>"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(html, encoding="utf-8")


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))

