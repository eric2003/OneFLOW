"""Tests for parameter-file actions in the project window."""

from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication, QLabel, QPushButton, QTreeWidget

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

            catalog_home = root.parent / "oneflow-user-data"
            with patch.dict(os.environ, {"ONEFLOW_HOME": str(catalog_home)}):
                window = ProjectWindow(
                    open_project(root),
                    root.parent,
                )
            self.assertTrue((catalog_home / "parameter-catalog.sqlite3").is_file())
            try:
                labels = [label.text() for label in window.findChildren(QLabel)]
                self.assertNotIn("Runtime file", labels)
                parameter_tree = window.findChild(QTreeWidget, "unifiedParameterTree")
                self.assertIsNotNone(parameter_tree)
                self.assertEqual(parameter_tree.topLevelItemCount(), 4)
                parameter = self._find_item(parameter_tree, "simuTask")
                self.assertIsNotNone(parameter)
                self.assertEqual(parameter.text(1), "flow")
                self.assertIsNone(self._find_item(parameter_tree, "control.txt"))
                button = window.findChild(QPushButton, "modifyParametersButton")
                self.assertIsNotNone(button)
                self.assertTrue(button.isEnabled())
            finally:
                window.close()

    @staticmethod
    def _find_item(tree: QTreeWidget, text: str):
        def visit(parent):
            for index in range(parent.childCount()):
                child = parent.child(index)
                if child.text(0) == text:
                    return child
                found = visit(child)
                if found is not None:
                    return found
            return None

        for index in range(tree.topLevelItemCount()):
            root = tree.topLevelItem(index)
            if root.text(0) == text:
                return root
            found = visit(root)
            if found is not None:
                return found
        return None

    def test_unified_view_merges_manifest_files_and_exposes_consumer_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            script = root / "script"
            script.mkdir()
            (script / "control.txt").write_text(
                '"simu.txt"\n"cfd.txt"\n"override.txt"\n', encoding="ascii"
            )
            (script / "simu.txt").write_text('string simuTask = "Solve";\nint maxSteps = 10;\n', encoding="ascii")
            (script / "cfd.txt").write_text("int maxSteps = 20;\n", encoding="ascii")
            (script / "override.txt").write_text("int MAXSTEPS = 30;\n", encoding="ascii")
            catalog_home = root.parent / "oneflow-user-data"
            with patch.dict(os.environ, {"ONEFLOW_HOME": str(catalog_home)}):
                window = ProjectWindow(open_project(root), root.parent)
            try:
                tree = window.findChild(QTreeWidget, "unifiedParameterTree")
                effective = self._find_item(tree, "MAXSTEPS")
                self.assertIsNotNone(effective)
                self.assertEqual(effective.text(1), "30")
                self.assertIn("file 3", effective.text(3))
                overridden = [
                    effective.child(index)
                    for index in range(effective.childCount())
                    if effective.child(index).text(0) == "[overridden] maxSteps"
                ]
                self.assertEqual([item.text(1) for item in overridden], ["10", "20"])
                self.assertIn("file 1", overridden[0].text(3))
                self.assertIn("file 2", overridden[1].text(3))
                consumer = next(
                    effective.child(index)
                    for index in range(effective.childCount())
                    if effective.child(index).text(0) == "Runtime consumer"
                )
                self.assertEqual(consumer.text(1), "Iteration::Init")
                self.assertIn("Iteration.cpp", consumer.text(3))
                self.assertIsNone(self._find_item(tree, "control.txt"))
            finally:
                window.close()

    def test_modify_selected_parameter_uses_configured_editor(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            script = root / "script"
            script.mkdir()
            (script / "control.txt").write_text('"simu.txt"\n', encoding="ascii")
            (script / "simu.txt").write_text("string simuTask = \"flow\";\n", encoding="ascii")
            catalog_home = root.parent / "oneflow-user-data"
            with patch.dict(os.environ, {"ONEFLOW_HOME": str(catalog_home)}):
                window = ProjectWindow(open_project(root), root.parent)
            try:
                tree = window.findChild(QTreeWidget, "unifiedParameterTree")
                tree.setCurrentItem(self._find_item(tree, "simuTask"))
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
