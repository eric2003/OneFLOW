"""Tests for explicit source-editor discovery and process launch."""
from __future__ import annotations
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from tools.oneflow_gui.source_editor import open_source_files, resolve_editor_executable

class SourceEditorTest(unittest.TestCase):
    def test_vscode_launch_uses_reuse_window_and_explicit_file_arguments(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); source=root/"sample.py"; source.write_text("pass\n", encoding="ascii")
            executable=root/"Code.exe"
            with patch("tools.oneflow_gui.source_editor.resolve_editor_executable", return_value=executable), patch(
                "tools.oneflow_gui.source_editor.subprocess.Popen") as popen:
                open_source_files("vscode", (source,), working_directory=root)
            popen.assert_called_once_with([str(executable), "--reuse-window", str(source.resolve())], cwd=str(root.resolve()))

    def test_visual_studio_launches_multiple_files(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); files=(root/"a.cpp",root/"b.py")
            for path in files: path.write_text("", encoding="ascii")
            executable=root/"devenv.exe"
            with patch("tools.oneflow_gui.source_editor.resolve_editor_executable", return_value=executable), patch(
                "tools.oneflow_gui.source_editor.subprocess.Popen") as popen:
                open_source_files("vs2026", files, working_directory=root)
            self.assertEqual(popen.call_args.args[0], [str(executable), *(str(path.resolve()) for path in files)])

    def test_missing_editor_is_reported(self):
        with patch("tools.oneflow_gui.source_editor.shutil.which", return_value=None), patch.dict(
            "tools.oneflow_gui.source_editor.os.environ", {}, clear=True):
            with self.assertRaisesRegex(FileNotFoundError, "Visual Studio Code was not found"):
                resolve_editor_executable("vscode")

    def test_missing_file_is_rejected_before_launch(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(FileNotFoundError, "no longer exist"):
                open_source_files("vscode", (Path(directory)/"missing.py",), working_directory=Path(directory))

if __name__ == "__main__":
    unittest.main()
