"""Backwards-compatible alias for :mod:`pyirena.api.dispatch`.

The dispatcher moved to ``pyirena/api/`` when a second transport (the ZMQ
service) started using it — it never depended on MCP, and living under
``mcp/`` implied otherwise. Importing it from here keeps working; new code
should import ``pyirena.api.dispatch``.
"""
from __future__ import annotations

from pyirena.api.calculator_schemas import CALCULATOR_SCHEMA_BY_NAME  # noqa: F401
from pyirena.api.control.schemas import TOOL_SCHEMA_BY_NAME  # noqa: F401
from pyirena.api.data_op_schemas import DATA_OP_SCHEMA_BY_NAME  # noqa: F401
from pyirena.api.dispatch import (  # noqa: F401
    _CATEGORY_BLURBS,
    _CATEGORY_BY_MODULE,
    _REGISTRY,
    _SOURCES,
    PROFILE_ALL,
    PROFILE_JSON_IMAGES,
    PROFILE_JSON_ONLY,
    SESSION_LIFECYCLE_NAMES,
    Profile,
    call_tool,
    describe_tool,
    list_categories,
    list_tools,
    resolve_profile,
    tool_flags,
)

__all__ = [
    "SESSION_LIFECYCLE_NAMES",
    "Profile",
    "PROFILE_ALL",
    "PROFILE_JSON_ONLY",
    "PROFILE_JSON_IMAGES",
    "resolve_profile",
    "tool_flags",
    "list_categories",
    "list_tools",
    "describe_tool",
    "call_tool",
]
