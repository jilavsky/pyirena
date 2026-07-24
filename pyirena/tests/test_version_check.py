"""
Tests for pyirena.version_check — GitHub-based update notification.

Covers version comparison (including beta-vs-final ordering, which the
project's own current version, e.g. 1.1.0b2, depends on) and the network
fetch failing silently under every error condition.
"""

import socket
from unittest.mock import patch

from pyirena.version_check import (
    _parse_version,
    is_newer,
    fetch_latest_release_tag,
    should_check_now,
)


def test_final_beats_its_own_beta():
    assert is_newer("1.1.0", "1.1.0b2")
    assert not is_newer("1.1.0b2", "1.1.0")


def test_simple_patch_ordering():
    assert is_newer("1.0.2", "1.0.1")
    assert not is_newer("1.0.1", "1.0.2")


def test_equal_versions_not_newer():
    assert not is_newer("1.0.1", "1.0.1")
    assert not is_newer("1.1.0b2", "1.1.0b2")


def test_prerelease_ordering():
    # a < b < rc < final, for the same release tuple
    assert is_newer("1.1.0b1", "1.1.0a3")
    assert is_newer("1.1.0rc1", "1.1.0b9")


def test_leading_v_is_stripped():
    assert is_newer("v1.2.0", "1.1.0")


def test_unparsable_input_never_raises_or_looks_newer():
    assert _parse_version("not-a-version") == ((0, 0, 0), 0, 0)
    assert not is_newer("not-a-version", "1.0.0")
    assert is_newer("1.0.0", "not-a-version")


def test_fetch_returns_tag_on_success():
    class _FakeResponse:
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def read(self):
            return b'{"tag_name": "v1.2.3"}'

    with patch("urllib.request.urlopen", return_value=_FakeResponse()):
        assert fetch_latest_release_tag() == "1.2.3"


def test_fetch_returns_none_on_timeout():
    with patch("urllib.request.urlopen", side_effect=socket.timeout("timed out")):
        assert fetch_latest_release_tag() is None


def test_fetch_returns_none_on_malformed_response():
    class _FakeResponse:
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def read(self):
            return b"not json"

    with patch("urllib.request.urlopen", return_value=_FakeResponse()):
        assert fetch_latest_release_tag() is None


def test_fetch_returns_none_when_offline():
    with patch("urllib.request.urlopen", side_effect=OSError("network unreachable")):
        assert fetch_latest_release_tag() is None


class _FakeStateManager:
    def __init__(self, last_check=""):
        self._last_check = last_check

    def get(self, tool, key, default=None):
        return self._last_check if key == "last_update_check" else default


def test_should_check_now_when_never_checked():
    assert should_check_now(_FakeStateManager(""))


def test_should_check_now_respects_weekly_interval():
    from datetime import datetime, timedelta

    recent = (datetime.now() - timedelta(days=1)).isoformat()
    stale = (datetime.now() - timedelta(days=8)).isoformat()
    assert not should_check_now(_FakeStateManager(recent))
    assert should_check_now(_FakeStateManager(stale))
