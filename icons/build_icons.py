"""Build the shared app-icon set (direction A, pictogram family) for Jan's SAS tools.

Kept in pyirena/icons/ as the single source for all five apps (pyIrena, Matilda, Anisa,
Bernardyn, AIDA). Change a glyph or color here, re-run, and copy the outputs into each repo:
    <repo>/icons/                      <app>_icon.svg, <app>.icns, <app>.ico
    <package GUI folder>/.../app_icon.png   (loaded by the package's gui/icon.py)

Writes, per app, into <out>/<app>/:
    <app>_icon.svg        master, 1024 macOS grid (824 px tile, 100 px margin)
    app_icon.png          1024 px render of the master (used by Qt setWindowIcon)
    <app>.icns            macOS bundle icon; 16/32 px slots use the small variant
    <app>.ico             Windows icon, full-bleed tile; 16/24/32 px use the small variant
The small variant drops the faint axes and thickens strokes so it survives 16 px.

Needs rsvg-convert (librsvg) on PATH and Pillow.

Usage:
    python build_icons.py [--out DIR]
"""

from __future__ import annotations

import argparse
import subprocess
from io import BytesIO
from pathlib import Path

from PIL import Image

APPS = {
    "pyIrena": "#2D5BD0",
    "Matilda": "#0E7F73",
    "Anisa": "#6A45C2",
    "Bernardyn": "#B4570F",
    "Aida": "#B8325F",
}


def spark(cx: float, cy: float, r: float, color: str) -> str:
    """Four-point spark with deeply concave sides (so it never reads as a 5-point star)."""
    a, b = 0.10 * r, 0.24 * r
    d = (
        f"M{cx} {cy - r} "
        f"C{cx + a} {cy - b} {cx + b} {cy - a} {cx + r} {cy} "
        f"C{cx + b} {cy + a} {cx + a} {cy + b} {cx} {cy + r} "
        f"C{cx - a} {cy + b} {cx - b} {cy + a} {cx - r} {cy} "
        f"C{cx - b} {cy - a} {cx - a} {cy - b} {cx} {cy - r} Z"
    )
    return f'<path d="{d}" fill="{color}"/>'


def glyph(app: str, small: bool) -> str:
    """White glyph in a 100x100 box."""
    k = 1.3 if small else 1.0
    c = "#FFFFFF"
    s = f'fill="none" stroke="{c}" stroke-linecap="round" stroke-linejoin="round"'
    axes = "" if small else f'<path d="M24 22 V78 H80" {s} stroke-width="4" opacity="0.55"/>'
    if app == "pyIrena":
        return (
            axes
            + f'<path d="M28 30 C42 31 47 36 52 46 C58 60 66 70 78 72" {s} stroke-width="{7 * k}"/>'
        )
    if app == "Matilda":
        return (
            f'<circle cx="50" cy="50" r="14" {s} stroke-width="{6 * k}"/>'
            f'<circle cx="50" cy="50" r="29" {s} stroke-width="{6 * k}"/>'
            f'<circle cx="50" cy="50" r="{4 * k}" fill="{c}"/>'
        )
    if app == "Anisa":
        return (
            f'<g transform="rotate(-35 50 50)">'
            f'<ellipse cx="50" cy="50" rx="35" ry="13" {s} stroke-width="{6 * k}"/>'
            f'<ellipse cx="50" cy="50" rx="21" ry="9" {s} stroke-width="{5 * k}"/></g>'
            f'<circle cx="50" cy="50" r="3" fill="{c}"/>'
        )
    if app == "Bernardyn":
        return (
            axes
            + f'<path d="M30 28 C44 34 56 50 78 56" {s} stroke-width="{6 * k}"/>'
            + f'<path d="M30 46 C44 52 56 64 78 70" {s} stroke-width="{6 * k}" opacity="0.75"/>'
        )
    if app == "Aida":
        return (
            f'<path d="M18 78 C34 78 46 72 54 64 C62 56 70 50 84 50" {s} stroke-width="{6 * k}"/>'
            + spark(44, 40, 27, c)
        )
    raise ValueError(app)


def svg(app: str, small: bool = False, bleed: bool = False) -> str:
    """Icon SVG. bleed=False: macOS grid (824 tile, 100 margin). bleed=True: tile fills the
    canvas, for Windows .ico where the macOS margin would only make the icon look small."""
    color = APPS[app]
    if bleed:
        tile = f'<rect width="1024" height="1024" rx="225" fill="{color}"/>'
        g = f'<g transform="scale(10.24)">{glyph(app, small)}</g>'
    else:
        tile = f'<rect x="100" y="100" width="824" height="824" rx="185" fill="{color}"/>'
        g = f'<g transform="translate(100 100) scale(8.24)">{glyph(app, small)}</g>'
    return (
        '<svg xmlns="http://www.w3.org/2000/svg" width="1024" height="1024" viewBox="0 0 1024 1024">'
        f"<title>{app}</title>{tile}{g}</svg>"
    )


def render(svg_text: str, size: int) -> Image.Image:
    png = subprocess.run(
        ["rsvg-convert", "-w", str(size), "-h", str(size)],
        input=svg_text.encode(),
        capture_output=True,
        check=True,
    ).stdout
    return Image.open(BytesIO(png)).convert("RGBA")


def build(app: str, out: Path) -> Path:
    """Write the full icon set for *app* into out/<app>/ and return that folder."""
    d = out / app
    d.mkdir(parents=True, exist_ok=True)
    big, small = svg(app), svg(app, small=True)
    (d / f"{app}_icon.svg").write_text(big, encoding="utf-8")
    render(big, 1024).save(d / "app_icon.png")
    # ICNS: Pillow takes the largest image plus append_images for the other slots.
    sizes = [16, 32, 64, 128, 256, 512, 1024]
    imgs = [render(small if s <= 32 else big, s) for s in sizes]
    imgs[-1].save(d / f"{app}.icns", append_images=imgs[:-1])
    ico_sizes = [16, 24, 32, 48, 64, 128, 256]
    ico_big, ico_small = svg(app, bleed=True), svg(app, small=True, bleed=True)
    ico_imgs = [render(ico_small if s <= 32 else ico_big, s) for s in ico_sizes]
    ico_imgs[-1].save(
        d / f"{app}.ico", sizes=[(s, s) for s in ico_sizes], append_images=ico_imgs[:-1]
    )
    return d


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=Path(__file__).resolve().parent / "out")
    args = parser.parse_args(argv)
    for name in APPS:
        print(f"Wrote {build(name, args.out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
