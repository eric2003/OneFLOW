"""Tests for cross-platform filesystem path helpers."""

from __future__ import annotations

import ctypes
import os
import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.path_utils import (
    is_path_within,
    normalize_path,
    paths_equal,
    relative_path_within,
)


class PathUtilsTest(unittest.TestCase):
    def test_normalize_path_expands_and_resolves(self) -> None:
        # Keep the temporary directory on the checkout drive: Windows cannot
        # compute a relative path between different drive letters.
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as temporary:
            root = Path(temporary)
            relative = Path(os.path.relpath(root, Path.cwd()))
            self.assertEqual(normalize_path(relative), root.resolve())

    def test_paths_equal_for_same_existing_path(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)
            self.assertTrue(paths_equal(path, path / "."))

    def test_paths_equal_for_nonexistent_descendants(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            left = root / "not-created" / "mesh.msh"
            right = root / "." / "not-created" / "mesh.msh"
            self.assertFalse(left.exists())
            self.assertTrue(paths_equal(left, right))

    @unittest.skipUnless(os.name == "nt", "Windows path case-folding only")
    def test_paths_equal_ignores_windows_case(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)
            self.assertTrue(paths_equal(path, Path(str(path).swapcase())))

    @unittest.skipUnless(os.name == "nt", "Windows 8.3 aliases only")
    def test_paths_equal_handles_windows_short_path_alias(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            long_directory = Path(temporary)
            buffer = ctypes.create_unicode_buffer(32768)
            get_short_path = ctypes.windll.kernel32.GetShortPathNameW
            get_short_path.argtypes = [ctypes.c_wchar_p, ctypes.c_wchar_p, ctypes.c_uint]
            get_short_path.restype = ctypes.c_uint
            length = get_short_path(str(long_directory), buffer, len(buffer))
            if not length or length >= len(buffer):
                self.skipTest("Windows did not provide an 8.3 short path")
            short_directory = Path(buffer.value)
            if os.path.normcase(str(short_directory)) == os.path.normcase(str(long_directory)):
                self.skipTest("8.3 short names are disabled for this temporary directory")

            child = long_directory / "mesh.msh"
            child.write_text("mesh", encoding="utf-8")
            self.assertTrue(paths_equal(child, short_directory / child.name))
            missing_mesh = short_directory / "not-created" / "mesh.msh"
            self.assertTrue(
                paths_equal(
                    long_directory / "not-created" / "mesh.msh",
                    missing_mesh,
                )
            )
            self.assertTrue(is_path_within(missing_mesh, long_directory))
            self.assertEqual(
                relative_path_within(missing_mesh, long_directory),
                Path("not-created/mesh.msh"),
            )

    def test_paths_equal_handles_symlink_alias(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / "target"
            target.mkdir()
            alias = root / "alias"
            try:
                alias.symlink_to(target, target_is_directory=True)
            except (OSError, NotImplementedError) as error:
                self.skipTest(f"Directory symlinks are unavailable: {error}")
            self.assertTrue(paths_equal(target, alias))

    def test_is_path_within_accepts_root_and_descendants(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.assertTrue(is_path_within(root, root))
            self.assertTrue(is_path_within(root / "stage" / "mesh.msh", root))

    def test_is_path_within_rejects_sibling_and_parent(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "workspace"
            sibling = Path(temporary) / "workspace-other"
            self.assertFalse(is_path_within(sibling / "mesh.msh", root))
            self.assertFalse(is_path_within(root.parent, root))


    def test_relative_path_within_rejects_escape(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            with self.assertRaises(ValueError):
                relative_path_within(base / "outside" / "mesh.msh", base / "workspace")

    def test_is_path_within_rejects_symlink_escape(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            root = base / "workspace"
            outside = base / "outside"
            root.mkdir()
            outside.mkdir()
            link = root / "external"
            try:
                link.symlink_to(outside, target_is_directory=True)
            except (OSError, NotImplementedError) as error:
                self.skipTest(f"Directory symlinks are unavailable: {error}")
            self.assertFalse(is_path_within(link / "mesh.msh", root))


if __name__ == "__main__":
    unittest.main()
