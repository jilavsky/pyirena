#!/usr/bin/env python3
"""Convert a validation results table (Markdown) into a Word document.

    python3 validationData/md_to_docx.py VALIDATION_RESULTS_paper.md
    python3 validationData/md_to_docx.py VALIDATION_RESULTS.md -o supplement.docx

Produces landscape US Letter with real Word tables (so the nine-column results
tables fit and stay editable), built-in heading styles (so Word's navigation
pane and any table of contents work), and Times New Roman throughout for a
manuscript draft.

Markdown handled: ATX headings, pipe tables, block quotes, bullet lists, fenced
code blocks, horizontal rules, HTML comments (dropped), and inline `code`,
**bold**, *italic* and [links](target).  That is the whole vocabulary these
generated tables use; it is not a general Markdown implementation.

Requires Node with the `docx` package (`npm install docx`).
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

BODY_FONT = "Times New Roman"
BODY_PT = 10
TABLE_PT = 8

# Landscape US Letter, 0.6 in margins, in DXA (1440 = 1 inch)
PAGE_W, PAGE_H, MARGIN = 12240, 15840, 864
CONTENT_W = PAGE_H - 2 * MARGIN          # landscape: the long edge is the width

#: Column widths for the nine-column results table, in DXA, summing to CONTENT_W.
RESULTS_WIDTHS = [4112, 1150, 1500, 1500, 1050, 850, 950, 1750, 1250]


# ---------------------------------------------------------------------------
# Markdown -> document model
# ---------------------------------------------------------------------------

INLINE = re.compile(r"(\*\*.+?\*\*|`[^`]+`|\*[^*]+?\*|\[[^\]]+\]\([^)]+\))")


def runs(text, bold=False, italic=False):
    """Split inline markup into styled runs.

    Recursive, so emphasis nested inside emphasis is handled: `**a `b` c**`
    yields a bold run, a bold monospace run and a bold run rather than one bold
    run with the backticks left in it.
    """
    out = []
    for piece in INLINE.split(text):
        if not piece:
            continue
        if piece.startswith("**") and piece.endswith("**") and len(piece) > 4:
            out.extend(runs(piece[2:-2], True, italic))
        elif piece.startswith("`") and piece.endswith("`"):
            out.append({"t": piece[1:-1], "code": True, "b": bold, "i": italic})
        elif piece.startswith("*") and piece.endswith("*") and len(piece) > 2:
            out.extend(runs(piece[1:-1], bold, True))
        elif piece.startswith("["):
            out.append({"t": piece[1:piece.index("]")], "b": bold, "i": italic})
        else:
            out.append({"t": piece, "b": bold, "i": italic})
    return out or [{"t": ""}]


def plain(rs):
    """Drop monospace styling — used for headings, which keep the body font."""
    return [{k: v for k, v in r.items() if k != "code"} for r in rs]


def split_row(line):
    return [c.strip() for c in line.strip().strip("|").split("|")]


def parse(md: str) -> list:
    lines = md.splitlines()
    blocks, i = [], 0
    while i < len(lines):
        line = lines[i]

        if line.startswith("<!--"):
            while i < len(lines) and "-->" not in lines[i]:
                i += 1
            i += 1
            continue

        if line.startswith("```"):
            i += 1
            code = []
            while i < len(lines) and not lines[i].startswith("```"):
                code.append(lines[i])
                i += 1
            i += 1
            blocks.append({"type": "code", "lines": code})
            continue

        if re.match(r"^(---+|\*\*\*+)\s*$", line):
            blocks.append({"type": "hr"})
            i += 1
            continue

        m = re.match(r"^(#{1,4})\s+(.*)$", line)
        if m:
            blocks.append({"type": f"h{len(m.group(1))}",
                           "runs": plain(runs(m.group(2).strip()))})
            i += 1
            continue

        if line.startswith("|") and i + 1 < len(lines) and re.match(r"^\|[-: |]+\|$", lines[i + 1]):
            header = split_row(line)
            i += 2
            rows = []
            while i < len(lines) and lines[i].startswith("|"):
                rows.append(split_row(lines[i]))
                i += 1
            blocks.append({"type": "table", "header": header, "rows": rows})
            continue

        if line.startswith(">"):
            buf = []
            while i < len(lines) and (lines[i].startswith(">") or
                                      (buf and lines[i].strip() and not lines[i].startswith(("|", "#", ">")))):
                buf.append(lines[i].lstrip("> ").rstrip())
                i += 1
            blocks.append({"type": "quote", "runs": runs(" ".join(buf).strip())})
            continue

        if re.match(r"^[*-]\s+", line):
            buf = [re.sub(r"^[*-]\s+", "", line).rstrip()]
            i += 1
            while i < len(lines) and lines[i].startswith("  ") and lines[i].strip():
                buf.append(lines[i].strip())
                i += 1
            blocks.append({"type": "bullet", "runs": runs(" ".join(buf))})
            continue

        if re.match(r"^\d+\.\s+", line):
            buf = [re.sub(r"^\d+\.\s+", "", line).rstrip()]
            i += 1
            while i < len(lines) and lines[i].startswith("   ") and lines[i].strip():
                buf.append(lines[i].strip())
                i += 1
            blocks.append({"type": "numbered", "runs": runs(" ".join(buf))})
            continue

        if not line.strip():
            i += 1
            continue

        buf = []
        while i < len(lines) and lines[i].strip() and not lines[i].startswith(("|", "#", ">", "```", "---")) \
                and not re.match(r"^([*-]|\d+\.)\s+", lines[i]):
            buf.append(lines[i].strip())
            i += 1
        blocks.append({"type": "p", "runs": runs(" ".join(buf))})
    return blocks


# ---------------------------------------------------------------------------
# Document model -> .docx, via the docx npm package
# ---------------------------------------------------------------------------

BUILDER = r"""
const fs = require('fs');
const d = require('docx');
const {Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
       HeadingLevel, WidthType, ShadingType, AlignmentType, BorderStyle,
       PageOrientation, LevelFormat, Footer, PageNumber} = d;

const cfg    = JSON.parse(fs.readFileSync(process.argv[2], 'utf8'));
const blocks = cfg.blocks;
const F = cfg.font, PT = cfg.body_pt * 2, TPT = cfg.table_pt * 2;

const mkRuns = (rs, size, bold) => rs.map(r => new TextRun({
  text: r.t,
  bold: bold || !!r.b,
  italics: !!r.i,
  font: r.code ? 'Consolas' : F,
  size: size,
}));

const HEAD = {h1: HeadingLevel.TITLE, h2: HeadingLevel.HEADING_1,
              h3: HeadingLevel.HEADING_2, h4: HeadingLevel.HEADING_3};

function widthsFor(nCols) {
  if (nCols === cfg.results_widths.length) return cfg.results_widths;
  const each = Math.floor(cfg.content_w / nCols);
  const w = new Array(nCols).fill(each);
  w[0] += cfg.content_w - each * nCols;
  return w;
}

function cell(text, widths, idx, header) {
  return new TableCell({
    width: {size: widths[idx], type: WidthType.DXA},
    shading: header ? {type: ShadingType.CLEAR, fill: 'E8E8E8'} : undefined,
    margins: {top: 40, bottom: 40, left: 80, right: 80},
    children: [new Paragraph({
      spacing: {before: 0, after: 0},
      alignment: idx === 0 ? AlignmentType.LEFT : AlignmentType.RIGHT,
      children: mkRuns(cfg.parseInline ? [{t: text}] : [{t: text}], TPT, header),
    })],
  });
}

const children = [];
for (const b of blocks) {
  if (b.type === 'hr') {
    children.push(new Paragraph({
      spacing: {before: 160, after: 160},
      border: {bottom: {style: BorderStyle.SINGLE, size: 6, color: 'AAAAAA'}},
      children: [new TextRun({text: '', font: F, size: PT})],
    }));
  } else if (HEAD[b.type]) {
    children.push(new Paragraph({
      heading: HEAD[b.type],
      spacing: {before: b.type === 'h1' ? 0 : 240, after: 120},
      children: mkRuns(b.runs, b.type === 'h1' ? PT + 10 : PT + 4, true),
    }));
  } else if (b.type === 'table') {
    const widths = widthsFor(b.header.length);
    const rows = [new TableRow({
      tableHeader: true,
      children: b.header.map((t, i) => cell(t.replace(/\*\*/g, ''), widths, i, true)),
    })];
    for (const r of b.rows) {
      rows.push(new TableRow({
        children: widths.map((_, i) => cell((r[i] || '').replace(/`/g, ''), widths, i, false)),
      }));
    }
    children.push(new Table({
      columnWidths: widths,
      width: {size: cfg.content_w, type: WidthType.DXA},
      rows,
    }));
    children.push(new Paragraph({spacing: {after: 120},
      children: [new TextRun({text: '', font: F, size: PT})]}));
  } else if (b.type === 'quote') {
    children.push(new Paragraph({
      spacing: {before: 60, after: 60}, indent: {left: 360},
      children: mkRuns(b.runs.map(r => Object.assign({}, r, {i: true})), PT - 2),
    }));
  } else if (b.type === 'bullet' || b.type === 'numbered') {
    children.push(new Paragraph({
      numbering: {reference: b.type === 'bullet' ? 'bullets' : 'numbers', level: 0},
      spacing: {before: 40, after: 40},
      children: mkRuns(b.runs, PT),
    }));
  } else if (b.type === 'code') {
    for (const ln of b.lines) {
      children.push(new Paragraph({
        spacing: {before: 0, after: 0}, indent: {left: 360},
        children: [new TextRun({text: ln, font: 'Consolas', size: PT - 2})],
      }));
    }
    children.push(new Paragraph({spacing: {after: 120},
      children: [new TextRun({text: '', font: F, size: PT})]}));
  } else {
    children.push(new Paragraph({
      spacing: {before: 60, after: 60},
      children: mkRuns(b.runs, PT),
    }));
  }
}

const doc = new Document({
  styles: {default: {document: {run: {font: F, size: PT}}}},
  numbering: {config: [
    {reference: 'bullets', levels: [{level: 0, format: LevelFormat.BULLET, text: '•',
      style: {paragraph: {indent: {left: 460, hanging: 260}}}}]},
    {reference: 'numbers', levels: [{level: 0, format: LevelFormat.DECIMAL, text: '%1.',
      style: {paragraph: {indent: {left: 460, hanging: 260}}}}]},
  ]},
  sections: [{
    footers: {default: new Footer({children: [new Paragraph({
      alignment: AlignmentType.CENTER,
      children: [new TextRun({children: [PageNumber.CURRENT], font: F, size: PT - 2})],
    })]})},
    properties: {page: {
      size: {width: cfg.page_w, height: cfg.page_h, orientation: PageOrientation.LANDSCAPE},
      margin: {top: cfg.margin, bottom: cfg.margin, left: cfg.margin, right: cfg.margin},
    }},
    children,
  }],
});

Packer.toBuffer(doc).then(buf => fs.writeFileSync(process.argv[3], buf));
"""


def build(md_path: Path, out_path: Path) -> None:
    blocks = parse(md_path.read_text())
    cfg = dict(blocks=blocks, font=BODY_FONT, body_pt=BODY_PT, table_pt=TABLE_PT,
               page_w=PAGE_W, page_h=PAGE_H, margin=MARGIN,
               content_w=CONTENT_W, results_widths=RESULTS_WIDTHS)

    with tempfile.TemporaryDirectory() as tmp:
        cfg_file = Path(tmp) / "doc.json"
        js_file = Path(tmp) / "build.js"
        cfg_file.write_text(json.dumps(cfg))
        js_file.write_text(BUILDER)
        env = dict(os.environ)
        # the docx package may be installed in the home directory rather than here
        env["NODE_PATH"] = os.pathsep.join(
            filter(None, [env.get("NODE_PATH"), str(Path.home() / "node_modules")]))
        subprocess.run(["node", str(js_file), str(cfg_file), str(out_path)],
                       check=True, env=env)

    n_tables = sum(1 for b in blocks if b["type"] == "table")
    n_rows = sum(len(b["rows"]) for b in blocks if b["type"] == "table")
    print(f"{out_path.name}: {len(blocks)} blocks, {n_tables} tables, {n_rows} data rows")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    out = None
    if "-o" in sys.argv:
        out = Path(sys.argv[sys.argv.index("-o") + 1])
        args = [a for a in args if a != str(out)]
    if len(args) != 1:
        sys.exit(__doc__)
    src = Path(args[0])
    if not src.is_absolute():
        src = Path(__file__).resolve().parent / src
    build(src, out or src.with_suffix(".docx"))
