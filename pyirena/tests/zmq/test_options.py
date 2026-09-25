"""Tests for the ZMQ service's options surface.

Three layers — config file, CLI, per-call — with a strict rule: a per-call
override may narrow what the operator allowed and never widen it. The value of
the whole design is that ``server_info`` tells an agent what a deployment
permits, so the effective options and what is reported must agree.
"""
from __future__ import annotations

import json

import pytest

from pyirena.zmq.options import (
    DEFAULT_PORT,
    PER_CALL_OPTIONS,
    ServerOptions,
    load_config_file,
    options_from_args,
)

# ---------------------------------------------------------------------------
# Defaults and derived values
# ---------------------------------------------------------------------------

def test_defaults_are_the_agreed_contract():
    o = ServerOptions()
    assert o.bind == f"tcp://0.0.0.0:{DEFAULT_PORT}"
    assert DEFAULT_PORT == 9865
    # Inside the orchestrator's 60 s client timeout, so the server answers
    # rather than letting the caller's socket die.
    assert o.request_budget_s < 60
    assert o.profile == "json_only"
    assert o.allow_files is False and o.allow_images is False
    assert o.validate() is None


@pytest.mark.parametrize(
    ("allow_files", "allow_images", "profile"),
    [(False, False, "json_only"), (False, True, "json_images"),
     (True, False, "all"), (True, True, "all")],
)
def test_profile_is_derived_from_the_flags(allow_files, allow_images, profile):
    """Derived, not stored: a profile field could contradict the flags."""
    o = ServerOptions(allow_files=allow_files, allow_images=allow_images,
                      data_root="/data" if allow_files else None)
    assert o.profile == profile
    assert o.to_dict()["profile"] == profile


def test_allow_files_without_a_data_root_is_refused():
    assert "data_root" in ServerOptions(allow_files=True).validate()
    assert ServerOptions(allow_files=True, data_root="/data").validate() is None


@pytest.mark.parametrize("kwargs", [
    {"request_budget_s": 0},
    {"max_message_mb": 0},
    {"max_sessions": 0},
    {"max_points": 0},
    {"bind": "workstation:9865"},        # not a ZMQ endpoint
])
def test_nonsense_options_are_rejected_with_a_reason(kwargs):
    assert ServerOptions(**kwargs).validate()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def test_port_is_shorthand_for_bind():
    o, err = options_from_args(["--port", "9999"])
    assert err is None and o.bind == "tcp://0.0.0.0:9999"


def test_explicit_bind_beats_port():
    o, _ = options_from_args(["--port", "9999", "--bind", "tcp://127.0.0.1:1234"])
    assert o.bind == "tcp://127.0.0.1:1234"


def test_flags_map_onto_the_dataclass():
    o, err = options_from_args([
        "--allow-images", "--request-budget-s", "10", "--max-sessions", "4",
        "--include-arrays", "--log-level", "DEBUG",
    ])
    assert err is None
    assert o.allow_images and o.include_arrays
    assert o.request_budget_s == 10.0 and o.max_sessions == 4
    assert o.profile == "json_images"


def test_a_bad_combination_is_reported_not_raised():
    o, err = options_from_args(["--allow-files"])
    assert err and "data_root" in err


# ---------------------------------------------------------------------------
# Config file
# ---------------------------------------------------------------------------

def test_json_config_file(tmp_path):
    path = tmp_path / "zmq.json"
    path.write_text(json.dumps({"max_sessions": 2, "allow_images": True,
                                "request_budget_s": 30}))
    values, err = load_config_file(path)
    assert err is None and values["max_sessions"] == 2

    o, err = options_from_args(["--config", str(path)])
    assert err is None
    assert o.max_sessions == 2 and o.allow_images is True


def test_cli_overrides_the_config_file(tmp_path):
    path = tmp_path / "zmq.json"
    path.write_text(json.dumps({"max_sessions": 2, "request_budget_s": 30}))
    o, _ = options_from_args(["--config", str(path), "--max-sessions", "9"])
    assert o.max_sessions == 9
    assert o.request_budget_s == 30          # untouched by the CLI


def test_config_from_the_environment(tmp_path, monkeypatch):
    path = tmp_path / "zmq.json"
    path.write_text(json.dumps({"max_message_mb": 3}))
    monkeypatch.setenv("PYIRENA_ZMQ_CONFIG", str(path))
    o, err = options_from_args([])
    assert err is None and o.max_message_mb == 3


def test_an_unknown_key_in_the_config_is_an_error(tmp_path):
    path = tmp_path / "zmq.json"
    path.write_text(json.dumps({"max_sesions": 2}))     # typo
    _, err = load_config_file(path)
    assert err and "max_sesions" in err and "max_sessions" in err


def test_a_missing_or_broken_config_file_says_which(tmp_path):
    _, err = load_config_file(tmp_path / "absent.json")
    assert "not found" in err

    broken = tmp_path / "broken.json"
    broken.write_text("{not json")
    _, err = load_config_file(broken)
    assert "parse" in err


# ---------------------------------------------------------------------------
# Per-call overrides — narrowing only
# ---------------------------------------------------------------------------

def test_per_call_options_are_a_short_explicit_list():
    assert set(PER_CALL_OPTIONS) == {
        "include_arrays", "max_points", "request_budget_s", "allow_images",
    }


def test_narrowing_is_allowed():
    base = ServerOptions(allow_images=True, request_budget_s=55, include_arrays=True)
    narrowed, err = base.with_overrides({
        "allow_images": False, "request_budget_s": 5, "include_arrays": False,
    })
    assert err is None
    assert narrowed.allow_images is False
    assert narrowed.request_budget_s == 5
    assert narrowed.profile == "json_only"          # the derived value follows
    # The server's own options are untouched — overrides are per request.
    assert base.allow_images is True


def test_widening_is_refused():
    base = ServerOptions()
    for override in ({"allow_images": True}, {"request_budget_s": 600}):
        _, err = base.with_overrides(override)
        assert err and err["code"] == "BAD_OPTION"


def test_an_operator_only_option_cannot_be_set_per_call():
    _, err = ServerOptions().with_overrides({"max_sessions": 1000})
    assert err["code"] == "BAD_OPTION"
    assert "server_info" in err["suggestion"]


def test_wrong_types_are_refused():
    base = ServerOptions()
    for override in ({"include_arrays": "yes"}, {"max_points": "lots"},
                     {"request_budget_s": "soon"}, {"max_points": -1}):
        _, err = base.with_overrides(override)
        assert err and err["code"] == "BAD_OPTION"


def test_max_points_null_means_no_cap():
    narrowed, err = ServerOptions().with_overrides({"max_points": None})
    assert err is None and narrowed.max_points is None


def test_options_must_be_an_object():
    _, err = ServerOptions().with_overrides(["include_arrays"])
    assert err["code"] == "BAD_OPTION"


def test_no_overrides_returns_the_same_object():
    base = ServerOptions()
    assert base.with_overrides(None) == (base, None)
    assert base.with_overrides({}) == (base, None)
