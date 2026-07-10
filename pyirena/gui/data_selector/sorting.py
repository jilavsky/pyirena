"""
pyirena.gui.data_selector.sorting — file-list sort keys.

Split from the original monolithic data_selector.py (no behavior change).
"""

import re





# ── Filename sort-key extractors ───────────────────────────────────────────
def _sort_key_name(name: str) -> str:
    return name.lower()


def _sort_key_temperature(name: str) -> float:
    m = re.search(r'_(-?\d+(?:\.\d+)?)C(?=_|\.|$)', name, re.IGNORECASE)
    return float(m.group(1)) if m else float('inf')


def _sort_key_time(name: str) -> float:
    m = re.search(r'_(\d+(?:\.\d+)?)min(?=_|\.|$)', name, re.IGNORECASE)
    return float(m.group(1)) if m else float('inf')


def _sort_key_order(name: str) -> float:
    # Strip extension then scan _-segments right-to-left for a bare integer
    # (digits only).  Skips any suffix that contains letters, including
    # _merged, _mrg, _scaled, and unit-bearing tokens like _10min or _5C.
    stem = re.sub(r'\.[^.]+$', '', name)
    for part in reversed(stem.split('_')):
        if re.fullmatch(r'\d+', part):
            return float(part)
    return float('inf')


def _sort_key_pressure(name: str) -> float:
    m = re.search(r'_(\d+(?:\.\d+)?)PSI(?=_|\.|$)', name, re.IGNORECASE)
    return float(m.group(1)) if m else float('inf')


_SORT_KEYS = [
    _sort_key_name,        # 0 Filename A→Z
    _sort_key_name,        # 1 Filename Z→A
    _sort_key_temperature, # 2 Temperature ↑
    _sort_key_temperature, # 3 Temperature ↓
    _sort_key_time,        # 4 Time ↑
    _sort_key_time,        # 5 Time ↓
    _sort_key_order,       # 6 Order number ↑
    _sort_key_order,       # 7 Order number ↓
    _sort_key_pressure,    # 8 Pressure ↑
    _sort_key_pressure,    # 9 Pressure ↓
]
