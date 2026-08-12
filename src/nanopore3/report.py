"""Dependency-free deterministic HTML reporting."""

from __future__ import annotations

from collections import Counter
import csv
from html import escape
from pathlib import Path
from typing import Iterable, Mapping


def status_counts(rows: Iterable[Mapping[str, object]], field: str) -> Counter[str]:
    return Counter(str(row.get(field, "missing")) for row in rows)


def write_html_report(
    path: Path,
    *,
    title: str,
    sections: Mapping[str, Mapping[str, int]],
    provenance: Mapping[str, object],
) -> None:
    """Write a compact report from canonical summaries only."""

    blocks: list[str] = []
    for heading, values in sections.items():
        rows = "".join(
            f"<tr><td>{escape(str(key))}</td><td>{int(value)}</td></tr>"
            for key, value in sorted(values.items())
        )
        blocks.append(f"<h2>{escape(heading)}</h2><table><tr><th>Status</th><th>Count</th></tr>{rows}</table>")
    provenance_rows = "".join(
        f"<tr><td>{escape(str(k))}</td><td>{escape(str(v))}</td></tr>"
        for k, v in sorted(provenance.items())
    )
    html = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>{escape(title)}</title>
<style>body{{font:16px system-ui;max-width:960px;margin:2rem auto;padding:0 1rem;color:#18202a}}
table{{border-collapse:collapse;margin-bottom:2rem}}th,td{{border:1px solid #ccd4dd;padding:.4rem .7rem;text-align:left}}
th{{background:#edf3f7}}code{{background:#edf3f7;padding:.1rem .25rem}}</style></head>
<body><h1>{escape(title)}</h1>{''.join(blocks)}
<h2>Provenance</h2><table>{provenance_rows}</table></body></html>"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(html, encoding="utf-8")


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))

