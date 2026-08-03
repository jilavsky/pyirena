"""
pyirena.version_check — lightweight GitHub-based update notification.

Replaces the old Igor Pro APS-server version check: reads the latest *stable*
release straight from GitHub's public Releases API (pre-releases like betas
are excluded automatically by GitHub's "latest" designation), throttled to
roughly once a week and fully silent on any failure (offline, timeout,
rate-limited, malformed response, etc.) so it can never disrupt startup.

Stdlib-only by design (no `requests`/`httpx`, no `packaging`) so it works even
in a minimal, non-GUI install.
"""

import json
import re
import urllib.request
from datetime import datetime, timedelta
from typing import Optional

GITHUB_LATEST_RELEASE_URL = "https://api.github.com/repos/jilavsky/pyirena/releases/latest"

_VERSION_RE = re.compile(r"^(\d+)\.(\d+)\.(\d+)(?:(a|b|rc)(\d+))?")
_PRERELEASE_RANK = {"a": 0, "b": 1, "rc": 2}


def _parse_version(v: str) -> tuple:
    """
    Parse a 'X.Y.Z' or 'X.Y.Z{a|b|rc}N' string into a sortable tuple.

    A final release always sorts after any prerelease of the same X.Y.Z
    (so '1.1.0' > '1.1.0b2'). Unparsable input sorts as the oldest possible
    version so garbage never triggers a false update notice.
    """
    m = _VERSION_RE.match(v.strip().lstrip("vV"))
    if not m:
        return ((0, 0, 0), 0, 0)
    major, minor, micro, pre, pre_n = m.groups()
    release = (int(major), int(minor), int(micro))
    if pre is None:
        return (release, 1, 0)
    return (release, 0, _PRERELEASE_RANK.get(pre, 0) * 1000 + int(pre_n))


def is_newer(remote: str, local: str) -> bool:
    """Return True if `remote` version string is newer than `local`."""
    return _parse_version(remote) > _parse_version(local)


def fetch_latest_release_tag(timeout: float = 3.0) -> Optional[str]:
    """
    Return the latest stable release tag from GitHub (leading 'v' stripped),
    or None on any failure. Never raises.
    """
    try:
        request = urllib.request.Request(
            GITHUB_LATEST_RELEASE_URL,
            headers={
                "User-Agent": "pyIrena-update-check",
                "Accept": "application/vnd.github+json",
            },
        )
        with urllib.request.urlopen(request, timeout=timeout) as response:
            data = json.load(response)
        tag = data.get("tag_name", "")
        return tag.lstrip("vV") or None
    except Exception:
        return None


def should_check_now(state_manager, min_interval_days: int = 7) -> bool:
    """Return True if it's been at least `min_interval_days` since the last check."""
    last_check_str = state_manager.get("data_selector", "last_update_check", "") or ""
    if not last_check_str:
        return True
    try:
        last_check = datetime.fromisoformat(last_check_str)
    except ValueError:
        return True
    return datetime.now() - last_check >= timedelta(days=min_interval_days)
