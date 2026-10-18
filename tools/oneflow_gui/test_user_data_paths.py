"""Tests for centralized user data path resolution."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.user_data_paths import UserDataPaths


class UserDataPathsTest(unittest.TestCase):
    def test_default_paths_use_dot_oneflow_under_home(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            home = Path(directory)
            paths = UserDataPaths.discover(home=home, environ={})

            self.assertEqual(paths.root, (home / ".oneflow").resolve())
            self.assertEqual(paths.workflow_data, paths.root / "workflow-data")
            self.assertEqual(
                paths.workflow_definitions,
                paths.workflow_data / "definitions",
            )
            self.assertEqual(paths.workflow_runs, paths.workflow_data / "runs")
            self.assertEqual(paths.executions, paths.workflow_data / "executions")
            self.assertEqual(paths.settings_file, paths.root / "settings.json")
            self.assertEqual(paths.logs, paths.root / "logs")
            self.assertEqual(
                paths.mesh_layout_workspace,
                paths.root / "mesh-layout-workspace",
            )
            self.assertEqual(
                paths.mesh_layout_specifications,
                paths.mesh_layout_workspace / "specifications",
            )
            self.assertEqual(
                paths.mesh_layout_outputs,
                paths.mesh_layout_workspace / "layouts",
            )
            self.assertFalse(paths.root.exists())

    def test_oneflow_home_overrides_default_home(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "custom-data"
            paths = UserDataPaths.discover(
                home=Path(directory) / "ignored-home",
                environ={"ONEFLOW_HOME": str(root)},
            )

            self.assertEqual(paths.root, root.resolve())
            self.assertEqual(paths.workflow_data, root.resolve() / "workflow-data")

    def test_empty_override_keeps_default_home(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            home = Path(directory)
            paths = UserDataPaths.discover(
                home=home,
                environ={"ONEFLOW_HOME": "  "},
            )
            self.assertEqual(paths.root, (home / ".oneflow").resolve())


if __name__ == "__main__":
    unittest.main()
