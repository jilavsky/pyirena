"""Diagnose how many records in an Igor .pxp file pyirena can decode.

Mirrors the same retry logic that ``pyirena.io.pxp_to_nexus`` uses, so the
SKIPPED count printed here matches what the importer will report.

Usage::

    python scripts/diag_pxp.py /path/to/Big.pxp
    python scripts/diag_pxp.py /path/to/Big.pxp --show-folders
"""
from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

from igor2.packed import (
    setup_packed_file_record_header, _RECORD_TYPE,
    _UnknownRecord, PACKEDRECTYPE_MASK,
    _byte_order, _need_to_reorder_bytes,
)
from igor2.record.wave import WaveRecord
from igor2.record.folder import FolderStartRecord, FolderEndRecord

from pyirena.io.pxp_to_nexus import _patch_v7_wave_to_v5


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("pxp", help="Path to the .pxp file")
    ap.add_argument("--show-folders", action="store_true",
                    help="Print the top-level folder structure as discovered")
    args = ap.parse_args()

    pxp = Path(args.pxp)
    print(f"Loading {pxp} ({pxp.stat().st_size:,} bytes)\n")

    # Match the loader's byte-order detection (peek first 2 bytes).
    with open(pxp, "rb") as f:
        peek = f.read(2)
    initial_byte_order = "<" if peek[0] != 0 else ">"
    print(f"Detected byte order: {initial_byte_order}")
    header_struct = setup_packed_file_record_header(byte_order=initial_byte_order)
    byte_order = initial_byte_order
    print(f"Header size: {header_struct.size}")

    n_ok = 0
    n_skipped = 0
    n_recovered_v7 = 0
    by_type_ok: Counter[str] = Counter()
    by_type_skipped: Counter[str] = Counter()
    size_buckets_skipped: Counter[str] = Counter()
    errors_by_type: dict[str, list[str]] = {}
    total_size_skipped = 0
    total_size_ok = 0
    folder_stack: list[str] = []
    top_folders: Counter[str] = Counter()

    with open(pxp, "rb") as f:
        while True:
            pos = f.tell()
            b = f.read(header_struct.size)
            if len(b) < header_struct.size:
                break
            header = header_struct.unpack_from(b)
            if header["version"]:
                need = _need_to_reorder_bytes(header["version"])
                refined = _byte_order(need)
                if refined != byte_order:
                    byte_order = refined
                    header_struct = setup_packed_file_record_header(byte_order=byte_order)
                    header = header_struct.unpack_from(b)
            num_bytes = header["numDataBytes"]
            if num_bytes < 0 or num_bytes > 256 * 1024 * 1024:
                print(f"  bad size {num_bytes} at offset 0x{pos:x}; stopping")
                break
            data = bytes(f.read(num_bytes))
            if len(data) < num_bytes:
                print(f"  truncated at offset 0x{pos:x}; stopping")
                break

            rtype_id = header["recordType"] & PACKEDRECTYPE_MASK
            rtype = _RECORD_TYPE.get(rtype_id, _UnknownRecord)
            tname = rtype.__name__

            # Track folder context so failures can be attributed to a tree path
            if rtype is FolderStartRecord:
                try:
                    rec = rtype(header, data, byte_order=byte_order)
                    name = rec.null_terminated_text
                    if isinstance(name, bytes):
                        name = name.decode("utf-8", errors="replace")
                    folder_stack.append(name)
                    if len(folder_stack) == 1:
                        top_folders[name] += 1
                except Exception:
                    folder_stack.append("?")
            elif rtype is FolderEndRecord:
                if folder_stack:
                    folder_stack.pop()

            try:
                rec = rtype(header, data, byte_order=byte_order)
                n_ok += 1
                by_type_ok[tname] += 1
                total_size_ok += num_bytes
                continue
            except Exception as exc:
                # Match the loader's v7 retry path.
                if rtype is WaveRecord:
                    patched = _patch_v7_wave_to_v5(data, byte_order)
                    if patched is not data:
                        try:
                            WaveRecord(header, patched, byte_order=byte_order)
                            n_ok += 1
                            n_recovered_v7 += 1
                            by_type_ok[tname] += 1
                            total_size_ok += num_bytes
                            continue
                        except Exception as exc2:
                            exc = exc2

                n_skipped += 1
                by_type_skipped[tname] += 1
                total_size_skipped += num_bytes
                if num_bytes < 1024:
                    bucket = "<1 KB"
                elif num_bytes < 64 * 1024:
                    bucket = "1-64 KB"
                elif num_bytes < 1024 * 1024:
                    bucket = "64 KB - 1 MB"
                else:
                    bucket = ">1 MB"
                size_buckets_skipped[bucket] += 1
                key = f"{tname}:{type(exc).__name__}"
                msgs = errors_by_type.setdefault(key, [])
                short = f"{exc}".splitlines()[0][:120]
                if short not in msgs and len(msgs) < 3:
                    msgs.append(short)

    print(f"\nRecords OK:        {n_ok:6d}  ({total_size_ok/1e6:.1f} MB)")
    if n_recovered_v7:
        print(f"  of which v7→v5 recovered: {n_recovered_v7}")
    print(f"Records skipped:   {n_skipped:6d}  ({total_size_skipped/1e6:.1f} MB)")

    print("\nBy record type (OK):")
    for t, n in by_type_ok.most_common():
        print(f"  {n:6d}  {t}")
    print("\nBy record type (SKIPPED):")
    for t, n in by_type_skipped.most_common():
        print(f"  {n:6d}  {t}")
    print("\nSkipped record size distribution:")
    for b, n in size_buckets_skipped.most_common():
        print(f"  {n:6d}  {b}")
    print("\nSample error messages (first 3 distinct per type):")
    for key, msgs in errors_by_type.items():
        print(f"  [{key}]")
        for m in msgs:
            print(f"      {m}")
    if args.show_folders:
        print("\nTop-level folders seen:")
        for name, n in top_folders.most_common():
            print(f"  {n:4d}  {name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
