"""Server options for the ZMQ service — one definition, three ways to set them.

The project owner asked for "reasonable variation of options for the server",
so the options are a designed surface rather than flags accreted onto
``main()``. There is exactly one list of them (``ServerOptions``), and three
layers, each overriding the one above:

1. a **config file** (``--config``, or ``$PYIRENA_ZMQ_CONFIG``),
2. **CLI flags**, built from the same dataclass so they cannot drift,
3. **per-call ``options``** in the request envelope, for the few that make
   sense per request.

A per-call override may only narrow what the operator allowed: a request can
turn images off, never on, and can shorten its own deadline, never extend it.
Anything else would let a caller talk the server out of its own policy.

``server_info`` returns the effective values, so an agent discovers what a
given deployment permits instead of guessing — and the same orchestrator code
works against a locked-down and a permissive server.

Unknown keys are an error, never ignored: silently dropping an option the
caller believes it set is the kind of bug that is only found in the data,
weeks later.
"""
from __future__ import annotations

import argparse
import json
import os
from dataclasses import asdict, dataclass, fields, replace
from pathlib import Path
from typing import Any, Optional, Tuple

DEFAULT_PORT = 9865
DEFAULT_BIND = f"tcp://0.0.0.0:{DEFAULT_PORT}"

# The owner's client waits 60 s and then raises, and its REQ socket is dead
# after that. The server's own budget sits just inside that so the caller
# always gets an answer it can act on rather than a dead socket.
DEFAULT_REQUEST_BUDGET_S = 55.0


@dataclass(frozen=True)
class ServerOptions:
    """Everything the ZMQ service can be told to do differently."""

    # --- transport ---------------------------------------------------------
    bind: str = DEFAULT_BIND
    max_message_mb: int = 16

    # --- limits ------------------------------------------------------------
    request_budget_s: float = DEFAULT_REQUEST_BUDGET_S
    max_input_points: int = 100_000
    session_ttl_min: float = 60.0
    max_sessions: int = 32

    # --- what the caller gets back -----------------------------------------
    include_arrays: bool = False
    max_points: Optional[int] = 2000

    # --- what the service is allowed to expose -----------------------------
    # Both off by default: the service exists because the caller has no access
    # to this machine's filesystem, so a path or a rendered-PNG path is at
    # best useless and at worst misleading. allow_images returns base64 only.
    allow_images: bool = False
    allow_files: bool = False
    data_root: Optional[str] = None

    # --- logging -----------------------------------------------------------
    log_file: Optional[str] = None
    log_level: str = "INFO"
    log_stderr: bool = False

    # -----------------------------------------------------------------------

    @property
    def profile(self) -> str:
        """Dispatcher profile implied by allow_files / allow_images.

        Derived rather than stored: a ``profile`` field that could disagree
        with the two flags is a bug waiting to happen.
        """
        if self.allow_files:
            return "all"
        return "json_images" if self.allow_images else "json_only"

    @property
    def max_message_bytes(self) -> int:
        return int(self.max_message_mb * 1024 * 1024)

    def validate(self) -> Optional[str]:
        """Return a human-readable reason these options are unusable, or None."""
        if self.allow_files and not self.data_root:
            return "allow_files needs data_root: refusing to expose the whole filesystem."
        if self.request_budget_s <= 0:
            return "request_budget_s must be positive."
        if self.max_message_mb <= 0:
            return "max_message_mb must be positive."
        if self.max_sessions <= 0:
            return "max_sessions must be positive."
        if self.max_points is not None and self.max_points <= 0:
            return "max_points must be positive, or null for no cap."
        if not self.bind.startswith(("tcp://", "ipc://", "inproc://")):
            return f"bind must be a ZMQ endpoint such as {DEFAULT_BIND!r}, got {self.bind!r}."
        return None

    def to_dict(self) -> dict:
        d = asdict(self)
        d["profile"] = self.profile
        return d

    # --- per-call overrides -------------------------------------------------

    def with_overrides(self, overrides: Optional[dict]) -> Tuple["ServerOptions", Optional[dict]]:
        """Apply a request's ``options`` block. Returns (options, error or None)."""
        if not overrides:
            return self, None
        if not isinstance(overrides, dict):
            return self, {
                "code": "BAD_OPTION",
                "message": "'options' must be an object.",
                "suggestion": f"Valid per-call options: {', '.join(sorted(PER_CALL_OPTIONS))}.",
            }

        changes: dict[str, Any] = {}
        for key, value in overrides.items():
            rule = PER_CALL_OPTIONS.get(key)
            if rule is None:
                return self, {
                    "code": "BAD_OPTION",
                    "message": f"'{key}' cannot be set per call.",
                    "suggestion": (
                        f"Per-call options: {', '.join(sorted(PER_CALL_OPTIONS))}. "
                        "Everything else is set by the service operator; call "
                        "server_info to see the effective values."
                    ),
                }
            narrowed, problem = rule(getattr(self, key), value)
            if problem:
                return self, {
                    "code": "BAD_OPTION",
                    "message": f"option '{key}': {problem}",
                    "suggestion": "Call server_info to see what this deployment allows.",
                }
            changes[key] = narrowed
        return replace(self, **changes), None


# --- per-call override rules -----------------------------------------------
#
# Each rule takes (server_value, requested_value) and returns
# (value_to_use, problem_or_None). A rule may narrow the server's policy and
# never widen it.

def _rule_bool_off_only(server_value, requested):
    if not isinstance(requested, bool):
        return None, "expected true or false."
    if requested and not server_value:
        return None, "this server does not allow it; it can only be turned off per call."
    return requested, None


def _rule_bool(server_value, requested):
    if not isinstance(requested, bool):
        return None, "expected true or false."
    return requested, None


def _rule_budget(server_value, requested):
    if not isinstance(requested, (int, float)) or isinstance(requested, bool):
        return None, "expected a number of seconds."
    if requested <= 0:
        return None, "must be positive."
    if requested > server_value:
        return None, (
            f"the server's budget is {server_value:g} s; a call may shorten its "
            "own deadline but not extend it."
        )
    return float(requested), None


def _rule_max_points(server_value, requested):
    if requested is None:
        return None, None                       # no cap: allowed, arrays are capped by size anyway
    if not isinstance(requested, int) or isinstance(requested, bool):
        return None, "expected an integer, or null for no cap."
    if requested <= 0:
        return None, "must be positive, or null for no cap."
    return requested, None


PER_CALL_OPTIONS = {
    "include_arrays": _rule_bool,
    "max_points": _rule_max_points,
    "request_budget_s": _rule_budget,
    "allow_images": _rule_bool_off_only,
}

#: Options that are also arguments of some tools; when the caller does not
#: name them explicitly, the effective option supplies the default.
TOOL_ARGUMENT_OPTIONS = ("include_arrays", "max_points")


# --- construction -----------------------------------------------------------

def _coerce(name: str, raw: Any) -> Any:
    """Coerce a config-file value to the type its dataclass field declares."""
    field = {f.name: f for f in fields(ServerOptions)}[name]
    if raw is None:
        return None
    annotation = str(field.type)
    if "bool" in annotation and not isinstance(raw, bool):
        if isinstance(raw, str):
            return raw.strip().lower() in ("1", "true", "yes", "on")
        return bool(raw)
    if "int" in annotation and "Optional" not in annotation:
        return int(raw)
    if "float" in annotation:
        return float(raw)
    if name == "max_points":
        return int(raw)
    return raw


def load_config_file(path: str | Path) -> Tuple[dict, Optional[str]]:
    """Read a JSON or TOML options file. Returns (values, error or None)."""
    path = Path(path).expanduser()
    if not path.exists():
        return {}, f"Config file not found: {path}"
    text = path.read_text(encoding="utf-8")

    try:
        if path.suffix.lower() == ".toml":
            try:
                import tomllib  # Python 3.11+
            except ModuleNotFoundError:
                return {}, (
                    "TOML config needs Python 3.11 or newer; "
                    "use a .json config file instead."
                )
            values = tomllib.loads(text)
        else:
            values = json.loads(text)
    except Exception as exc:
        return {}, f"Could not parse {path}: {exc}"

    if not isinstance(values, dict):
        return {}, f"{path} must contain an object of option names to values."

    known = {f.name for f in fields(ServerOptions)}
    unknown = sorted(set(values) - known)
    if unknown:
        return {}, (
            f"Unknown option(s) in {path}: {', '.join(unknown)}. "
            f"Valid options: {', '.join(sorted(known))}."
        )
    return {k: _coerce(k, v) for k, v in values.items()}, None


def build_arg_parser() -> argparse.ArgumentParser:
    """CLI flags, generated from ServerOptions so the two cannot drift."""
    defaults = ServerOptions()
    p = argparse.ArgumentParser(
        prog="pyirena-zmq",
        description=(
            "Serve pyirena.api to a remote orchestrator over ZMQ REQ/REP, "
            "JSON in and JSON out. Options may also come from a config file "
            "(--config); CLI flags override it."
        ),
    )
    p.add_argument("--config", metavar="FILE",
                   help="JSON (or TOML, on Python 3.11+) file of options.")
    p.add_argument("--bind", help=f"ZMQ endpoint to bind (default {defaults.bind}).")
    p.add_argument("--port", type=int,
                   help=f"Shorthand for --bind tcp://0.0.0.0:PORT (default {DEFAULT_PORT}).")
    p.add_argument("--max-message-mb", type=int,
                   help=f"Reject messages larger than this (default {defaults.max_message_mb}).")
    p.add_argument("--request-budget-s", type=float,
                   help=("Reply TIMEOUT rather than exceed this many seconds "
                         f"(default {defaults.request_budget_s:g}); keep it under the "
                         "client's own timeout."))
    p.add_argument("--max-input-points", type=int,
                   help=f"Largest curve accepted (default {defaults.max_input_points}).")
    p.add_argument("--session-ttl-min", type=float,
                   help=f"Evict sessions idle this long (default {defaults.session_ttl_min:g}).")
    p.add_argument("--max-sessions", type=int,
                   help=f"Cap on open sessions (default {defaults.max_sessions}).")
    p.add_argument("--include-arrays", action="store_true", default=None,
                   help="Return model curves by default (callers can still ask per call).")
    p.add_argument("--max-points", type=int,
                   help=f"Cap on each returned array (default {defaults.max_points}).")
    p.add_argument("--allow-images", action="store_true", default=None,
                   help="Expose image tools, returning base64 PNGs (no server paths).")
    p.add_argument("--allow-files", action="store_true", default=None,
                   help="Expose file-reading and -writing tools; requires --data-root.")
    p.add_argument("--data-root", help="Confine all file access to this directory.")
    p.add_argument("--log-file", help="Rotating log file (default: the user log dir).")
    p.add_argument("--log-level", help=f"Logging level (default {defaults.log_level}).")
    p.add_argument("--log-stderr", action="store_true", default=None,
                   help="Also log to stderr. Never to stdout.")
    return p


def options_from_args(argv: Optional[list[str]] = None) -> Tuple[ServerOptions, Optional[str]]:
    """Build options from a config file and CLI flags. Returns (options, error)."""
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    values: dict[str, Any] = {}
    config_path = args.config or os.environ.get("PYIRENA_ZMQ_CONFIG")
    if config_path:
        values, error = load_config_file(config_path)
        if error:
            return ServerOptions(), error

    known = {f.name for f in fields(ServerOptions)}
    for name, value in vars(args).items():
        if name in ("config", "port") or value is None:
            continue
        if name in known:
            values[name] = value

    # --port is a convenience over --bind; an explicit --bind wins.
    if args.port is not None and not args.bind:
        values["bind"] = f"tcp://0.0.0.0:{args.port}"

    options = ServerOptions(**values)
    return options, options.validate()
