"""Tests for the cross-platform GUI launcher configuration."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.__main__ import _bytecode_cache_root


class BytecodeCacheRootTest(unittest.TestCase):
    def test_linux_uses_xdg_cache_home_when_absolute(self):
        with tempfile.TemporaryDirectory() as directory:
            root = _bytecode_cache_root(
                platform="linux",
                environ={"XDG_CACHE_HOME": directory},
                home=Path("/home/test-user"),
            )
            self.assertEqual(root, Path(directory) / "oneflow" / "pycache")

    def test_linux_defaults_to_user_cache_directory(self):
        root = _bytecode_cache_root(
            platform="linux",
            environ={},
            home=Path("/home/test-user"),
        )
        self.assertEqual(
            root, Path("/home/test-user/.cache/oneflow/pycache")
        )

    def test_linux_ignores_relative_xdg_cache_home(self):
        root = _bytecode_cache_root(
            platform="linux",
            environ={"XDG_CACHE_HOME": "relative-cache"},
            home=Path("/home/test-user"),
        )
        self.assertEqual(
            root, Path("/home/test-user/.cache/oneflow/pycache")
        )

    def test_windows_uses_local_app_data(self):
        root = _bytecode_cache_root(
            platform="win32",
            environ={"LOCALAPPDATA": "C:/Users/test/AppData/Local"},
            home=Path("C:/Users/test"),
        )
        self.assertEqual(
            root.as_posix(),
            "C:/Users/test/AppData/Local/OneFLOW/pycache",
        )

    def test_macos_uses_library_caches(self):
        root = _bytecode_cache_root(
            platform="darwin",
            environ={},
            home=Path("/Users/test-user"),
        )
        self.assertEqual(
            root, Path("/Users/test-user/Library/Caches/OneFLOW/pycache")
        )


if __name__ == "__main__":
    unittest.main()
