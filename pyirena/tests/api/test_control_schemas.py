"""Contract tests for the control-surface JSON schemas (api.control.schemas).

The control tools are exposed to LLM/automation clients two ways: directly via
``pyirena.api.control`` (each schema in ``TOOL_SCHEMAS`` describes one callable)
and over MCP as ``pyirena_ctrl_*`` tools. A schema that drifts from its callable
silently breaks schema-driven clients — this is exactly how the ``open_dataset``
``use_slit_smeared`` omission slipped through CI. These tests lock the contract:

* every schema is structurally well-formed;
* every schema name resolves to a real ``pyirena.api.control`` callable;
* schema properties/required match the callable's signature (parity);
* every schema is registered as an MCP tool (and the count is locked).
"""
from __future__ import annotations

import inspect

import pytest

import pyirena.api.control as ctrl
from pyirena.api.control.schemas import TOOL_SCHEMA_BY_NAME, TOOL_SCHEMAS

_JSON_TYPES = {"string", "number", "integer", "boolean", "array", "object", "null"}


def _sig_params(fn):
    """Named (non-var) parameters of *fn* as {name: has_default}."""
    sig = inspect.signature(fn)
    return {
        p.name: (p.default is not p.empty)
        for p in sig.parameters.values()
        if p.kind in (p.POSITIONAL_OR_KEYWORD, p.KEYWORD_ONLY)
    }


def test_schema_names_unique_and_indexed():
    names = [s["name"] for s in TOOL_SCHEMAS]
    assert len(names) == len(set(names)), "duplicate schema names in TOOL_SCHEMAS"
    assert set(TOOL_SCHEMA_BY_NAME) == set(names)
    # Lock the count so an accidentally-dropped schema is caught.
    assert len(names) == 51, f"expected 51 control schemas, found {len(names)}"


@pytest.mark.parametrize("schema", TOOL_SCHEMAS, ids=lambda s: s["name"])
def test_schema_is_well_formed(schema):
    assert isinstance(schema.get("name"), str) and schema["name"]
    assert isinstance(schema.get("description"), str) and schema["description"]
    ins = schema.get("input_schema")
    assert isinstance(ins, dict) and ins.get("type") == "object"
    props = ins.get("properties", {})
    assert isinstance(props, dict)
    for pname, pspec in props.items():
        assert isinstance(pspec, dict), f"{schema['name']}.{pname} spec not a dict"
        # 'type' may be a single JSON type or a list (nullable form, e.g.
        # ["integer", "null"]). Both must contain only valid JSON types.
        ptype = pspec.get("type")
        types = ptype if isinstance(ptype, list) else [ptype]
        assert types and all(t in _JSON_TYPES for t in types), \
            f"{schema['name']}.{pname} has invalid JSON type {ptype!r}"
    required = ins.get("required", [])
    assert isinstance(required, list)
    assert set(required) <= set(props), \
        f"{schema['name']}: required lists non-property {set(required) - set(props)}"


@pytest.mark.parametrize("schema", TOOL_SCHEMAS, ids=lambda s: s["name"])
def test_schema_matches_callable_signature(schema):
    fn = getattr(ctrl, schema["name"], None)
    assert callable(fn), f"schema '{schema['name']}' has no api.control callable"

    params = _sig_params(fn)
    mandatory = {n for n, has_default in params.items() if not has_default}

    ins = schema["input_schema"]
    props = set(ins.get("properties", {}))
    required = set(ins.get("required", []))

    # No schema property that isn't a real parameter.
    assert props <= set(params), \
        f"{schema['name']}: schema exposes non-parameters {props - set(params)}"
    # Required set must be exactly the parameters with no default.
    assert required == mandatory, (
        f"{schema['name']}: required {required} != mandatory params {mandatory}"
    )
    # Every real parameter must be exposed (this is the check that would have
    # caught the missing open_dataset 'use_slit_smeared' property).
    assert set(params) <= props, \
        f"{schema['name']}: parameters missing from schema {set(params) - props}"


def test_open_dataset_exposes_use_slit_smeared():
    """Regression guard for the specific gap the review found."""
    props = TOOL_SCHEMA_BY_NAME["open_dataset"]["input_schema"]["properties"]
    assert "use_slit_smeared" in props
    assert props["use_slit_smeared"]["type"] == "boolean"


def test_mcp_registers_all_tools_structurally():
    """Validate the *whole* registered MCP surface, not just the read tools.

    The old smoke test only name-checked 18 read/discovery tools while the
    server registers 68 (18 read + 50 control). Here we lock the counts and
    assert every registered tool is structurally complete (has a description
    and a well-formed object input schema).
    """
    pytest.importorskip("mcp")
    import asyncio

    from pyirena.mcp.server import mcp

    tools = asyncio.run(mcp.list_tools())
    names = [t.name for t in tools]
    ctrl_tools = [n for n in names if n.startswith("pyirena_ctrl_")]
    read_tools = [n for n in names if not n.startswith("pyirena_ctrl_")]

    # Locked counts — adding/removing a tool is an intentional change that must
    # update this test (mirrors the 51-schema lock above).
    assert len(names) == 68, f"expected 68 registered MCP tools, found {len(names)}"
    assert len(ctrl_tools) == 50, f"expected 50 control tools, found {len(ctrl_tools)}"
    assert len(read_tools) == 18, f"expected 18 read tools, found {len(read_tools)}"

    for t in tools:
        assert t.description, f"MCP tool {t.name} has no description"
        ins = getattr(t, "inputSchema", None)
        assert isinstance(ins, dict) and ins.get("type") == "object", \
            f"MCP tool {t.name} has a malformed inputSchema"
