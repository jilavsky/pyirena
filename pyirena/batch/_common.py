"""
pyirena.batch._common — Shared config/data-loading helpers for the batch API.

Split from the original monolithic pyirena/batch.py (no behavior change).
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Optional, Union



if TYPE_CHECKING:
    pass

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _load_config(config_file: Union[str, Path]) -> Optional[Dict]:
    """Load and validate a pyIrena JSON config file.  Returns None on failure."""
    config_file = Path(config_file)
    try:
        with open(config_file, 'r') as f:
            config = json.load(f)
    except Exception as e:
        log.error(f"[pyirena.batch] Cannot read config file '{config_file}': {e}")
        return None

    if '_pyirena_config' not in config:
        log.error(f"[pyirena.batch] '{config_file}' is not a pyIrena configuration file "
              f"(missing '_pyirena_config' header).")
        return None

    return config


def _load_data(data_file: Union[str, Path]) -> Optional[Dict]:
    """Load SAS data from a text (.dat/.txt) or HDF5 (.h5/.hdf5) file.

    Text files are converted to a cleaned NXcanSAS HDF5 sibling via
    ``ensure_nxcansas_sibling`` before loading, so all callers always
    receive ``is_nxcansas=True`` and a valid HDF5 filepath.  The sibling
    is cached by mtime and reused on subsequent calls.

    Returns a dict with keys: Q, Intensity, Error (may be None).
    """
    from pyirena.io.hdf5 import readGenericNXcanSAS

    data_file = Path(data_file)
    if not data_file.exists():
        log.error(f"[pyirena.batch] Data file not found: '{data_file}'")
        return None

    ext = data_file.suffix.lower()

    try:
        if ext in ('.txt', '.dat'):
            from pyirena.io.text_import import ensure_nxcansas_sibling
            h5_file = ensure_nxcansas_sibling(data_file)
            data = readGenericNXcanSAS(str(h5_file.parent), h5_file.name)
            actual_file = h5_file
        else:
            data = readGenericNXcanSAS(str(data_file.parent), data_file.name)
            actual_file = data_file
    except Exception as e:
        log.error(f"[pyirena.batch] Error reading '{data_file}': {e}")
        return None

    if data is None:
        log.error(f"[pyirena.batch] Could not read data from '{data_file}'")
        return None

    data['filepath'] = str(actual_file)
    data['is_nxcansas'] = True
    data.setdefault('label', data_file.stem)
    return data
