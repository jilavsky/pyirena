"""Safe path resolution for the pyirena.api surface.

The api functions accept user-supplied paths. When the api is exposed over
MCP (or any RPC), an upstream LLM may pass arbitrary strings. We restrict
file access to an optional PYIRENA_DATA_ROOT subtree so the agent cannot
read /etc/passwd or wander outside an authorised directory.

Without PYIRENA_DATA_ROOT set: paths are accepted as long as they resolve
to an absolute path and exist. With it set: paths are resolved against the
root and rejected if they escape via .. traversal or symlinks.

Read-only: there are no write operations in pyirena.api v1, so the threat
model is information disclosure only.
"""
from __future__ import annotations

import os
from pathlib import Path


class PathSecurityError(ValueError):
    """Raised when a supplied path escapes the configured PYIRENA_DATA_ROOT."""


def _data_root() -> Path | None:
    root = os.environ.get("PYIRENA_DATA_ROOT")
    if not root:
        return None
    return Path(root).expanduser().resolve()


def resolve_safe(path: str | os.PathLike, must_exist: bool = True) -> Path:
    """Resolve *path* to an absolute Path, enforcing PYIRENA_DATA_ROOT if set.

    Parameters
    ----------
    path : str or PathLike
        User-supplied path. Relative paths are resolved against the data
        root (when set) or against the current working directory.
    must_exist : bool
        If True, raises FileNotFoundError when the resolved path does not
        exist. Set False for paths that are about to be created (plot
        output destinations).

    Raises
    ------
    PathSecurityError
        If PYIRENA_DATA_ROOT is set and the resolved path is not inside it.
    FileNotFoundError
        If must_exist is True and the path does not exist.
    """
    root = _data_root()
    p = Path(path).expanduser()
    if not p.is_absolute() and root is not None:
        p = root / p
    p = p.resolve()

    if root is not None:
        try:
            p.relative_to(root)
        except ValueError as exc:
            raise PathSecurityError(
                f"Path {p} is outside PYIRENA_DATA_ROOT={root}"
            ) from exc

    if must_exist and not p.exists():
        raise FileNotFoundError(f"No such path: {p}")

    return p


def resolve_safe_folder(path: str | os.PathLike) -> Path:
    """Resolve *path* and verify it is an existing directory."""
    p = resolve_safe(path, must_exist=True)
    if not p.is_dir():
        raise NotADirectoryError(f"Not a directory: {p}")
    return p


def resolve_safe_file(path: str | os.PathLike) -> Path:
    """Resolve *path* and verify it is an existing file."""
    p = resolve_safe(path, must_exist=True)
    if not p.is_file():
        raise IsADirectoryError(f"Not a file: {p}")
    return p
