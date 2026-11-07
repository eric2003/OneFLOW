"""Tests for real script-directory discovery and optional control manifests."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.script_file_discovery import discover_script_parameter_files
from tools.oneflow_gui.project_repository import ProjectType, open_project


class ScriptFileDiscoveryTest(unittest.TestCase):
    def test_real_files_are_loaded_in_filename_order_and_manifest_is_checked(self):
        with tempfile.TemporaryDirectory() as temporary:
            script = Path(temporary) / "script"
            script.mkdir()
            (script / "control.txt").write_text(
                '"b.txt";\n"a.txt";\n', encoding="ascii"
            )
            (script / "a.txt").write_text("int a = 1;\n", encoding="ascii")
            (script / "b.txt").write_text("int b = 2;\n", encoding="ascii")
            (script / "extra.txt").write_text("int extra = 3;\n", encoding="ascii")

            result = discover_script_parameter_files(script)

            self.assertEqual(result.mode, "manifest")
            self.assertEqual(
                [path.name for path in result.files],
                ["a.txt", "b.txt", "extra.txt"],
            )
            self.assertTrue(any("extra.txt" in item for item in result.diagnostics))
            self.assertTrue(any("filename order" in item for item in result.diagnostics))

    def test_missing_and_duplicate_manifest_entries_are_diagnosed(self):
        with tempfile.TemporaryDirectory() as temporary:
            script = Path(temporary) / "script"
            script.mkdir()
            (script / "control.txt").write_text(
                '"present.txt";\n"missing.txt";\n"present.txt";\n',
                encoding="ascii",
            )
            (script / "present.txt").write_text("int present = 1;\n", encoding="ascii")

            result = discover_script_parameter_files(script)

            self.assertEqual([path.name for path in result.files], ["present.txt"])
            self.assertTrue(any("missing.txt" in item for item in result.diagnostics))
            self.assertTrue(any("more than once" in item for item in result.diagnostics))
            self.assertFalse(any(item.startswith("ERROR:") for item in result.diagnostics))

    def test_parameter_style_control_is_loaded_as_an_ordinary_parameter_file(self):
        with tempfile.TemporaryDirectory() as temporary:
            script = Path(temporary) / "script"
            script.mkdir()
            (script / "control.txt").write_text(
                'string simutask = "Grid";\nint dimension = 2;\n',
                encoding="ascii",
            )
            (script / "grid.txt").write_text("int gridScale = 1;\n", encoding="ascii")

            result = discover_script_parameter_files(script)

            self.assertEqual(result.mode, "parameter-file")
            self.assertEqual(
                [path.name for path in result.files], ["control.txt", "grid.txt"]
            )
            self.assertTrue(any("ordinary parameter file" in item for item in result.diagnostics))

    def test_missing_control_discovers_all_real_parameter_files(self):
        with tempfile.TemporaryDirectory() as temporary:
            script = Path(temporary) / "script"
            script.mkdir()
            (script / "simu.txt").write_text(
                'string simutask = "Solve";\n', encoding="ascii"
            )
            (script / "cfd.txt").write_text("int maxSteps = 20;\n", encoding="ascii")
            (script / "solver.txt").write_text("NsSolver\n", encoding="ascii")

            result = discover_script_parameter_files(script)

            self.assertEqual(result.mode, "directory")
            self.assertEqual([path.name for path in result.files], ["cfd.txt", "simu.txt"])
            self.assertTrue(any("solver.txt" in item for item in result.diagnostics))

    def test_project_without_control_still_loads_parameters_and_detects_task(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            script = root / "script"
            script.mkdir()
            (script / "simu.txt").write_text(
                'string simutask = "Solve";\n', encoding="ascii"
            )
            (script / "cfd.txt").write_text("int maxSteps = 20;\n", encoding="ascii")

            project = open_project(root)

            self.assertEqual(project.project_type, ProjectType.FLOW)
            self.assertEqual(
                [path.name for path in project.parameter_files], ["cfd.txt", "simu.txt"]
            )
            self.assertIsNone(project.control_file)


    def test_project_loads_parameter_style_control_and_uses_its_task(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            script = root / "script"
            script.mkdir()
            (script / "control.txt").write_text(
                'string simutask = "Grid";\nint dimension = 2;\n',
                encoding="ascii",
            )
            (script / "grid.txt").write_text(
                "int gridScale = 1;\n", encoding="ascii"
            )

            project = open_project(root)

            self.assertEqual(project.project_type, ProjectType.GRID)
            self.assertEqual(project.parameter_file_mode, "parameter-file")
            self.assertEqual(
                [path.name for path in project.parameter_files],
                ["control.txt", "grid.txt"],
            )



if __name__ == "__main__":
    unittest.main()
