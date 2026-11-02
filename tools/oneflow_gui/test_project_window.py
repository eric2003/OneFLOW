"""Tests for parameter-file actions in the project window."""

from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication, QLabel, QListWidget, QPushButton

from tools.oneflow_gui.project_repository import open_project
from tools.oneflow_gui.project_window import ProjectWindow


class ProjectWindowParameterTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.application = QApplication.instance() or QApplication([])

    def test_runtime_file_is_not_repeated_and_parameter_files_are_selectable(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            script = root / "script"
            script.mkdir()
            (script / "control.txt").write_text('"simu.txt"\n', encoding="ascii")
            (script / "simu.txt").write_text('string simuTask = "flow";\n', encoding="ascii")

            window = ProjectWindow(
                open_project(root),
                root.parent,
            )
            try:
                labels = [label.text() for label in window.findChildren(QLabel)]
                self.assertNotIn("Runtime file", labels)
                parameter_list = window.findChild(QListWidget, "parameterFileList")
                self.assertIsNotNone(parameter_list)
                self.assertEqual(parameter_list.count(), 1)
                self.assertEqual(
                    Path(parameter_list.item(0).text()),
                    (script / "simu.txt").resolve(),
                )
                button = window.findChild(QPushButton, "modifyParametersButton")
                self.assertIsNotNone(button)
                self.assertTrue(button.isEnabled())
            finally:
                window.close()

    def test_modify_selected_parameter_uses_configured_editor(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            script = root / "script"
            script.mkdir()
            (script / "control.txt").write_text('"simu.txt"\n', encoding="ascii")
            (script / "simu.txt").write_text("string simuTask = \"flow\";\n", encoding="ascii")
            window = ProjectWindow(open_project(root), root.parent)
            try:
                with patch(
                    "tools.oneflow_gui.project_window.ApplicationSettings.load"
                ) as load_settings, patch(
                    "tools.oneflow_gui.project_window.open_source_files"
                ) as open_files:
                    load_settings.return_value.source_editor = "vscode"
                    window._modify_selected_parameter()
                open_files.assert_called_once_with(
                    "vscode",
                    ((script / "simu.txt").resolve(),),
                    working_directory=root.resolve(),
                )
            finally:
                window.close()


if __name__ == "__main__":
    unittest.main()
