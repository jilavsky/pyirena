"""Number formatting utilities for the pyirena GUI.

Two functions cover all display contexts:

eng_fmt(v, sig=4)
    Engineering notation for display-only widgets (labels, graph annotations,
    report tables).  Values in [0.001, 999] are shown as plain decimals; outside
    that range the exponent is chosen as the nearest multiple of 3 and the value
    is rendered as  "m.mmm×10^N"  (where N ∈ {…, -6, -3, 3, 6, …}).

eng_fmt_edit(v, sig=4)
    Same engineering-exponent logic but emits parseable text ("m.mmme+N") for
    editable QLineEdit / wheel-scroll fields, so Python's float() can round-trip
    the displayed string without modification.
"""

import math
import numpy as np


def eng_fmt(v: float, sig: int = 4) -> str:
    """Format *v* for display using engineering notation (×10^N, N multiple of 3).

    Values whose magnitude is in [0.001, 999] are shown as plain decimals
    (no exponent).  Everything else uses engineering notation with the
    multiplication-sign form, e.g. '1.234×10^-6'.
    """
    if v == 0.0:
        return '0'
    try:
        if not np.isfinite(v):
            return str(v)
    except (TypeError, ValueError):
        return str(v)

    abs_v = abs(v)
    if 0.001 <= abs_v < 1000:
        return f'{v:.{sig}g}'

    exp = int(math.floor(math.log10(abs_v) / 3)) * 3
    mantissa = v / (10.0 ** exp)
    return f'{mantissa:.{sig}g}×10^{exp}'


def eng_fmt_edit(v: float, sig: int = 4) -> str:
    """Format *v* for editable fields using parseable engineering notation.

    Like eng_fmt() but emits 'me+N' / 'me-N' style text that Python's float()
    can parse directly.  Values in [0.001, 999] are still shown as plain
    decimals.
    """
    if v == 0.0:
        return '0'
    try:
        if not np.isfinite(v):
            return str(v)
    except (TypeError, ValueError):
        return str(v)

    abs_v = abs(v)
    if 0.001 <= abs_v < 1000:
        return f'{v:.{sig}g}'

    exp = int(math.floor(math.log10(abs_v) / 3)) * 3
    mantissa = v / (10.0 ** exp)
    sign = '+' if exp >= 0 else '-'
    return f'{mantissa:.{sig}g}e{sign}{abs(exp):02d}'
