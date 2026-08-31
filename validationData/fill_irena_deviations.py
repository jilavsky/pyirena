#!/usr/bin/env python3
"""Fill the ``Irena dev %`` column of a validation results table.

    python3 validationData/fill_irena_deviations.py <table>.md
    python3 validationData/fill_irena_deviations.py <table>.md --export-csv

The Irena values are entered by hand after analysing the same files in the Igor
Pro Irena package.  The first form computes the last column in place, so it can
be re-run whenever more Irena numbers are filled in; the second additionally
writes the table's Irena column back into ``irena_values.csv``, which is what
``run_validation_report.py`` reads, so the numbers survive in the repository as
data rather than only inside a document.

    Irena dev % = 100 (Irena - true) / true

Rows are left blank when the Irena cell is empty, when it is a marker such as
``*`` (a quantity Irena defines differently and that is therefore not directly
comparable), or when the true value is zero (the model-versus-exact-curve rows,
whose target is zero by construction, so a relative deviation is undefined).

The script also prints agreement statistics between pyIrena and Irena over
every row where both packages reported a number.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

NUM = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$")

#: "2.3206-7" — a number whose exponent lost its "e" on the way out of Igor or
#: into the table.  Unambiguous (a bare mantissa cannot be followed by a signed
#: integer), so it is repaired and the cell is rewritten in explicit scientific
#: notation rather than being silently skipped.
LOST_E = re.compile(r"^([+-]?(?:\d+\.?\d*|\.\d+))([-+]\d+)$")

#: Irena values entered without their exponent.  Key is (dataset, quantity);
#: value is (value exactly as entered, value to use instead, why).  The fix is
#: applied only when the cell still holds the raw entered value, so re-running
#: this script cannot apply it twice.
SCALE_FIXES = {
    ("simple_porod", "Kp"): (
        2.5137, 2.5137e-06,
        "Irena value entered as 2.5137; read as 2.5137e-06, since Irena reports "
        "the Porod constant in the same units as the true value (cm^-1 A^-4)",
    ),
}


def _num(cell):
    """Parse a table cell, repairing an exponent that lost its "e".

    Returns (value, repaired_text or None).
    """
    cell = cell.strip().replace("**", "")
    if NUM.match(cell):
        return float(cell), None
    m = LOST_E.match(cell)
    if m:
        fixed = f"{m.group(1)}e{m.group(2)}"
        return float(fixed), fixed
    return None, None


def process(path: Path) -> None:
    lines = path.read_text().splitlines()
    out, dataset = [], ""
    filled = skipped = 0
    pairs = []          # (dataset, quantity, true, pyirena, irena)
    applied_fixes = []

    for line in lines:
        m = re.match(r"^### `([^`]+)`", line)
        if m:
            dataset = m.group(1)

        if not (line.startswith("|") and line.count("|") == 10
                and not line.startswith("|---")
                and "quantity | unit" not in line):
            out.append(line)
            continue

        cells = [c.strip() for c in line.split("|")[1:-1]]
        quantity, _unit, true_s, py_s, _dev, _tol, _ok, irena_s, _ = cells

        true_v, _ = _num(true_s)
        py_v, _ = _num(py_s)
        irena_v, repaired = _num(irena_s)
        if repaired is not None:
            cells[7] = repaired
            applied_fixes.append((dataset, quantity,
                                  f"exponent had lost its 'e': {irena_s} read as {repaired}"))
        fix = SCALE_FIXES.get((dataset, quantity))
        if fix and irena_v == fix[0]:
            irena_v = fix[1]
            cells[7] = f"{irena_v:.6g}"
            applied_fixes.append((dataset, quantity, fix[2]))

        if irena_v is None or true_v is None or true_v == 0:
            cells[8] = ""
            if irena_s and irena_s != "*":
                skipped += 1
        else:
            cells[8] = f"{100.0 * (irena_v - true_v) / true_v:+.3f}"
            filled += 1
            if py_v is not None:
                pairs.append((dataset, quantity, true_v, py_v, irena_v))

        out.append("| " + " | ".join(cells) + " |")

    path.write_text("\n".join(out) + "\n")

    print(f"{path.name}: filled {filled} deviations"
          + (f", {skipped} left blank" if skipped else ""))
    for ds, q, note in applied_fixes:
        print(f"  scale fix applied: {ds} / {q} — {note}")

    if not pairs:
        return
    both = [(abs(100 * (p - t) / t), abs(100 * (i - t) / t),
             abs(100 * (p - i) / i) if i else float("nan"), ds, q)
            for ds, q, t, p, i in pairs]
    n = len(both)
    med = lambda xs: sorted(xs)[n // 2]                          # noqa: E731
    print(f"\n{n} quantities reported by both packages:")
    print(f"  median |pyIrena - true|  = {med([b[0] for b in both]):.3f} %")
    print(f"  median |Irena  - true|   = {med([b[1] for b in both]):.3f} %")
    print(f"  median |pyIrena - Irena| = {med([b[2] for b in both]):.3f} %")
    worst = sorted(both, key=lambda b: -b[2])[:6]
    print("  largest pyIrena-vs-Irena differences:")
    for _, _, d, ds, q in worst:
        print(f"    {d:8.2f} %  {ds} / {q}")


def export_csv(path: Path, csv_path: Path) -> None:
    """Write the table's Irena column into irena_values.csv.

    Existing entries for other datasets are preserved, so exporting from a
    table that covers part of the set never discards the rest.
    """
    import csv as _csv

    existing = {}
    if csv_path.exists():
        with csv_path.open(newline="") as fh:
            for r in _csv.DictReader(fh):
                existing[(r["dataset"], r["quantity"])] = r

    added = updated = 0
    dataset = ""
    for line in path.read_text().splitlines():
        m = re.match(r"^### `([^`]+)`", line)
        if m:
            dataset = m.group(1)
            continue
        if not (line.startswith("|") and line.count("|") == 10
                and not line.startswith("|---") and "quantity | unit" not in line):
            continue
        c = [x.strip() for x in line.split("|")[1:-1]]
        quantity, unit, irena_s = c[0], c[1], c[7]
        if not irena_s:
            continue
        if irena_s == "*":
            row = dict(dataset=dataset, quantity=quantity, unit=unit,
                       irena_value="", status="not comparable",
                       note="Irena parameterises this differently; see irena_notes.md")
        else:
            value, _ = _num(irena_s)
            if value is None:
                continue
            row = dict(dataset=dataset, quantity=quantity, unit=unit,
                       irena_value=irena_s, status="", note="")
        key = (dataset, quantity)
        if key in existing and existing[key] == row:
            continue
        if key in existing:
            updated += 1
        else:
            added += 1
        existing[key] = row

    with csv_path.open("w", newline="") as fh:
        w = _csv.DictWriter(fh, ["dataset", "quantity", "unit", "irena_value",
                                 "status", "note"])
        w.writeheader()
        for row in existing.values():
            w.writerow(row)
    print(f"{csv_path.name}: {added} added, {updated} updated, "
          f"{len(existing)} entries total")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    flags = {a for a in sys.argv[1:] if a.startswith("--")}
    if len(args) != 1 or flags - {"--export-csv"}:
        sys.exit(__doc__)
    target = Path(args[0])
    if not target.is_absolute():
        target = Path(__file__).resolve().parent / target
    process(target)
    if "--export-csv" in flags:
        export_csv(target, Path(__file__).resolve().parent / "irena_values.csv")
