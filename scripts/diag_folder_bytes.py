"""Scan a .pxp for FolderStart records and dump raw bytes for any whose
payload contains a given substring. Used to confirm whether 'missing'
long-named folders are even present in the file and what their on-disk
encoding looks like.

Usage:
    python scripts/diag_folder_bytes.py /path/to/Big.pxp "crystalline"
"""
from __future__ import annotations

import argparse
from pathlib import Path

from igor2.packed import (
    setup_packed_file_record_header,
    PACKEDRECTYPE_MASK,
)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("pxp", help="Path to the .pxp file")
    ap.add_argument("needle", help="Substring to search for in folder-record payloads")
    ap.add_argument("--show-context", type=int, default=2,
                    help="How many records before/after to show (default 2)")
    args = ap.parse_args()

    pxp = Path(args.pxp)
    print(f"Loading {pxp} ({pxp.stat().st_size:,} bytes)\n")

    with open(pxp, "rb") as f:
        peek = f.read(2)
    initial = "<" if peek[0] != 0 else ">"
    hs = setup_packed_file_record_header(byte_order=initial)
    print(f"byte order: {initial}, header size: {hs.size}")

    needle = args.needle.encode("latin-1")
    needle_low = args.needle.lower().encode("latin-1")

    # First pass: collect all records as (index, offset, recordType, version, num_bytes, data)
    records: list[tuple[int, int, int, int, int, bytes]] = []
    with open(pxp, "rb") as f:
        idx = 0
        while True:
            pos = f.tell()
            b = f.read(hs.size)
            if len(b) < hs.size:
                break
            h = hs.unpack_from(b)
            n = h["numDataBytes"]
            if n < 0 or n > 256 * 1024 * 1024:
                break
            data = f.read(n)
            if len(data) < n:
                break
            records.append((idx, pos, h["recordType"] & PACKEDRECTYPE_MASK,
                             h.get("version", 0), n, data))
            idx += 1

    print(f"Read {len(records)} records total\n")

    # Second pass: find matches and show context
    matches = []
    for i, (idx, pos, rtype, ver, n, data) in enumerate(records):
        if needle in data or needle_low in data.lower():
            matches.append(i)

    if not matches:
        print(f"No record payload contains {args.needle!r}")
        # As a sanity check, also search for the FIRST 8 chars only — long
        # names might be wrapped or stored with prefixes.
        short_needle = args.needle.encode("latin-1")[:8]
        for i, (idx, pos, rtype, ver, n, data) in enumerate(records):
            if short_needle in data:
                print(f"  (but record {idx} contains the 8-char prefix {short_needle!r})")
                matches.append(i)
                break
        if not matches:
            return 0

    print(f"Found {len(matches)} matching record(s); showing first 5:\n")
    for hit in matches[:5]:
        lo = max(0, hit - args.show_context)
        hi = min(len(records), hit + args.show_context + 1)
        print(f"=== match at record index {hit} ===")
        for j in range(lo, hi):
            idx, pos, rtype, ver, n, data = records[j]
            tag = ">>>" if j == hit else "   "
            rtname = _rtype_name(rtype)
            preview_bytes = data[:80]
            preview_repr = preview_bytes.decode("latin-1", errors="replace")
            print(f"{tag} #{idx:5d} @0x{pos:08x} type={rtype:2d}({rtname:18s}) "
                  f"ver={ver:3d} n={n:6d}  data[:80]={preview_repr!r}")
            if j == hit and rtype != 9:
                # If the match isn't a folder, also show raw hex of first 64 bytes
                print(f"            raw hex: {data[:64].hex(' ')}")
        print()
    return 0


def _rtype_name(rt: int) -> str:
    return {
        0: "Unused", 1: "Variables", 2: "History", 3: "Wave",
        4: "Recreation", 5: "Procedure", 6: "Unused",
        7: "GetHistory", 8: "PackedFile",
        9: "FolderStart", 10: "FolderEnd",
    }.get(rt, f"Unknown{rt}")


if __name__ == "__main__":
    raise SystemExit(main())
