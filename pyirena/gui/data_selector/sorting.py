"""
pyirena.gui.data_selector.sorting — file-list sort keys.

The keys moved to :mod:`pyirena.core.file_sorting`, which every file browser
now shares (they had drifted into four near-copies).  This module re-exports
them under their historical private names so existing imports keep working;
new code should import from ``pyirena.core.file_sorting``.
"""

from pyirena.core.file_sorting import (
    SORT_KEYS as _SORT_KEYS,
)
from pyirena.core.file_sorting import (
    sort_key_name as _sort_key_name,
)
from pyirena.core.file_sorting import (
    sort_key_order as _sort_key_order,
)
from pyirena.core.file_sorting import (
    sort_key_pressure as _sort_key_pressure,
)
from pyirena.core.file_sorting import (
    sort_key_temperature as _sort_key_temperature,
)
from pyirena.core.file_sorting import (
    sort_key_time as _sort_key_time,
)

__all__ = [
    "_SORT_KEYS",
    "_sort_key_name",
    "_sort_key_temperature",
    "_sort_key_time",
    "_sort_key_order",
    "_sort_key_pressure",
]
