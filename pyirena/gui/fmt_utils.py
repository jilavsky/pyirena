"""Number formatting utilities — moved to :mod:`pyirena.core.fmt_utils`.

The formatters are pure numeric code (math + numpy, no Qt) and are needed by
``pyirena.core.reporting``, which ``core`` may not import from ``gui``.  They
now live in ``core``; this module re-exports them so the ~20 GUI call sites
keep working.  New code should import from ``pyirena.core.fmt_utils``.
"""

from pyirena.core.fmt_utils import eng_fmt, eng_fmt_edit

__all__ = ["eng_fmt", "eng_fmt_edit"]
