"""Tests for catalog-independent OneFLOW project discovery."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.project_repository import ProjectType, open_project


class ProjectRepositoryTest(unittest.TestCase):
    def test_opens_arbitrary_directory_without_catalog_registration(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory) / "user_project"
            root.mkdir()

            project = open_project(root)

            self.assertEqual(project.name, "user_project")
            self.assertEqual(project.root, root.resolve())
            self.assertEqual(project.project_type, ProjectType.UNKNOWN)
            self.assertIsNone(project.control_file)
            self.assertIsNone(project.runtime_file)

    def test_detects_grid_project_from_runtime_task(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            script = root / "script"
            script.mkdir()
            (script / "control.txt").write_text('"simu.txt";\n', encoding="utf-8")
            (script / "simu.txt").write_text(
                'string simutask = "Grid";\n', encoding="utf-8"
            )

            project = open_project(root)

            self.assertEqual(project.project_type, ProjectType.GRID)
            self.assertEqual(project.control_file, (script / "control.txt").resolve())
            self.assertEqual(project.runtime_file, (script / "simu.txt").resolve())

    def test_detects_flow_project_case_insensitively(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            (root / "control.txt").write_text('"flow.txt";\n', encoding="utf-8")
            (root / "flow.txt").write_text(
                'string simuTask = "Flow";\n', encoding="utf-8"
            )

            self.assertEqual(open_project(root).project_type, ProjectType.FLOW)

    def test_discovers_parameter_files_in_workflow_case_directories(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory) / "workflow"
            grid_script = root / "rae2822grid" / "script"
            flow_script = root / "rae2822_roe_sa" / "script"
            grid_script.mkdir(parents=True)
            flow_script.mkdir(parents=True)
            (root / "workflow.json").write_text("{}\n", encoding="utf-8")
            (grid_script / "control.txt").write_text(
                '"simu.txt";\n"dimension.txt";\n"grid.txt";\n',
                encoding="utf-8",
            )
            (grid_script / "simu.txt").write_text(
                'string simutask = "Grid";\n', encoding="utf-8"
            )
            (grid_script / "dimension.txt").write_text(
                "int dimension = 2;\n", encoding="utf-8"
            )
            (grid_script / "grid.txt").write_text(
                'string targetGridFileName = "results/grid.ofl";\n',
                encoding="utf-8",
            )
            (flow_script / "control.txt").write_text(
                '"simu.txt";\n"cfd.txt";\n', encoding="utf-8"
            )
            (flow_script / "simu.txt").write_text(
                'string simutask = "Solve";\n', encoding="utf-8"
            )
            (flow_script / "cfd.txt").write_text(
                'string gridFileName = "grid/grid.ofl";\n', encoding="utf-8"
            )

            project = open_project(root)

            self.assertEqual(project.control_file, (grid_script / "control.txt").resolve())
            self.assertEqual(project.runtime_file, (grid_script / "simu.txt").resolve())
            self.assertEqual(
                {path.name for path in project.parameter_files},
                {"simu.txt", "dimension.txt", "grid.txt", "cfd.txt"},
            )
            self.assertTrue(all(path.is_file() and path.stat().st_size for path in project.parameter_files))

    def test_rejects_a_missing_project_directory(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            missing = Path(temporary_directory) / "missing"
            with self.assertRaises(NotADirectoryError):
                open_project(missing)


if __name__ == "__main__":
    unittest.main()
