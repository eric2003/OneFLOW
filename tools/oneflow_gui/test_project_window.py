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

            catalog_home = root / "oneflow-user-data"
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
                self.assertEqual(parameter.text(4), "Both")
                runtime = next(
                    parameter.child(index)
                    for index in range(parameter.childCount())
                    if parameter.child(index).text(0) == "Runtime consumer"
                )
                self.assertEqual(runtime.text(1), "Not audited")
                self.assertIsNone(self._find_item(parameter_tree, "control.txt"))
                button = window.findChild(QPushButton, "modifyParametersButton")
                self.assertIsNotNone(button)
                self.assertTrue(button.isEnabled())
            finally:
                window.close()

    @staticmethod
    def _find_item(parent, text: str):
        def visit(item):
            if item.text(0) == text:
                return item
            for index in range(item.childCount()):
                found = visit(item.child(index))
                if found is not None:
                    return found
            return None

        if isinstance(parent, QTreeWidget):
            for index in range(parent.topLevelItemCount()):
                found = visit(parent.topLevelItem(index))
                if found is not None:
                    return found
            return None
        return visit(parent)

    def test_unified_view_merges_manifest_files_and_exposes_consumer_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            script = root / "script"
            script.mkdir()
            (script / "control.txt").write_text(
                '"simu.txt"\n"cfd.txt"\n"override.txt"\n', encoding="ascii"
            )
            (script / "simu.txt").write_text('string simuTask = "Solve";\nint maxSteps = 10;\nint time_integral = 1;\nint rk_stage = 3;\n', encoding="ascii")
            (script / "cfd.txt").write_text("int maxSteps = 20;\n", encoding="ascii")
            (script / "override.txt").write_text("int MAXSTEPS = 30;\nReal gridTrans[3] = {0.0, 0.0, 0.0};\n", encoding="ascii")
            (script / "extra.txt").write_text(
                "int extraManifestParameter = 9;\n", encoding="ascii"
            )
            catalog_home = root.parent / "oneflow-user-data"
            with patch.dict(os.environ, {"ONEFLOW_HOME": str(catalog_home)}):
                window = ProjectWindow(open_project(root), root.parent)
            try:
                tree = window.findChild(QTreeWidget, "unifiedParameterTree")
                effective = self._find_item(tree, "MAXSTEPS")
                self.assertIsNotNone(effective)
                self.assertEqual(effective.text(1), "30")
                self.assertEqual(effective.text(4), "Solve")
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
                extra = self._find_item(tree, "extraManifestParameter")
                self.assertIsNotNone(extra)
                self.assertEqual(extra.text(1), "9")
                self.assertIn("script/extra.txt", extra.text(3))
                for parameter_name in ("time_integral", "rk_stage"):
                    parameter_item = self._find_item(tree, parameter_name)
                    self.assertIsNotNone(parameter_item)
                    knowledge = self._find_item(parameter_item, "Formula & Knowledge")
                    self.assertIsNotNone(knowledge)
                    self.assertEqual(knowledge.text(1), "time_integral_rk_stage_weight")
                    if parameter_name == "time_integral":
                        tree.setCurrentItem(knowledge)
                        with patch(
                            "tools.oneflow_gui.verification_curve_dialog.VerificationCurveDialog"
                        ) as dialog_type:
                            window._modify_selected_parameter()
                        dialog_type.assert_called_once_with(root.parent.resolve(), window)
                        dialog_type.return_value.select_entry.assert_called_once_with(
                            "time_integral_rk_stage_weight"
                        )
                        dialog_type.return_value.exec.assert_called_once()
                unsupported = self._find_item(tree, "gridTrans")
                self.assertIsNotNone(unsupported)
                self.assertEqual(unsupported.text(4), "Unknown")
                self.assertIn("Array declarations", unsupported.toolTip(0))
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


    def test_project_without_control_loads_real_script_parameter_files(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            script = root / "script"
            script.mkdir()
            (script / "simu.txt").write_text(
                'string simutask = "Solve";\n', encoding="ascii"
            )
            (script / "cfd.txt").write_text(
                "int maxSteps = 40;\n", encoding="ascii"
            )
            catalog_home = root.parent / "oneflow-user-data"
            with patch.dict(os.environ, {"ONEFLOW_HOME": str(catalog_home)}):
                window = ProjectWindow(open_project(root), root.parent)
            try:
                tree = window.findChild(QTreeWidget, "unifiedParameterTree")
                self.assertIsNotNone(self._find_item(tree, "simutask"))
                maximum = self._find_item(tree, "maxSteps")
                self.assertIsNotNone(maximum)
                self.assertEqual(maximum.text(1), "40")
                self.assertIn("script/cfd.txt", maximum.text(3))
            finally:
                window.close()

    def test_parameter_style_control_is_visible_as_a_source_document(self) -> None:
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
            catalog_home = root.parent / "oneflow-user-data"
            with patch.dict(os.environ, {"ONEFLOW_HOME": str(catalog_home)}):
                window = ProjectWindow(open_project(root), root.parent)
            try:
                tree = window.findChild(QTreeWidget, "unifiedParameterTree")
                task = self._find_item(tree, "simutask")
                self.assertIsNotNone(task)
                self.assertIn("script/control.txt", task.text(3))
                self.assertEqual(task.text(1), "Grid")
                dimension = self._find_item(tree, "dimension")
                self.assertIsNotNone(dimension)
                self.assertEqual(dimension.text(1), "2")
            finally:
                window.close()



if __name__ == "__main__":
    unittest.main()
