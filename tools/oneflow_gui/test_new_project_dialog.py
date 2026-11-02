"""Regression tests for the new project template selector."""

from __future__ import annotations

import os
import unittest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication

from tools.oneflow_gui.new_project_dialog import NewProjectDialog
from tools.oneflow_gui.user_data_paths import UserDataPaths
from tools.oneflow_gui.project_templates import available_project_templates


class NewProjectDialogTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.application = QApplication.instance() or QApplication([])

    def test_default_destination_is_managed_under_oneflow_home(self) -> None:
        from pathlib import Path

        dialog = NewProjectDialog(Path.cwd())
        try:
            self.assertEqual(
                Path(dialog.destination_edit.text()),
                UserDataPaths.discover().projects,
            )
        finally:
            dialog.close()

    def test_laminar_template_uses_requested_default_project_name(self) -> None:
        from pathlib import Path

        dialog = NewProjectDialog(Path.cwd())
        try:
            index = dialog.template_combo.findData("laminar_flat_plate")
            dialog.template_combo.setCurrentIndex(index)
            self.assertEqual(dialog.name_edit.text(), "LaminarFlatPlate")

            dialog.name_edit.setText("MyCustomPlate")
            dialog.template_combo.setCurrentIndex(
                dialog.template_combo.findData("blank")
            )
            self.assertEqual(dialog.name_edit.text(), "MyCustomPlate")

            dialog.template_combo.setCurrentIndex(index)
            self.assertEqual(dialog.name_edit.text(), "MyCustomPlate")
        finally:
            dialog.close()

    def test_template_descriptions_do_not_replace_template_ids(self) -> None:
        from pathlib import Path

        dialog = NewProjectDialog(Path.cwd())
        try:
            templates = available_project_templates()
            self.assertEqual(dialog.template_combo.count(), len(templates))
            for index, template in enumerate(templates):
                self.assertEqual(
                    dialog.template_combo.itemData(index),
                    template.template_id,
                )
                self.assertEqual(
                    dialog.template_combo.itemData(index, Qt.ItemDataRole.ToolTipRole),
                    template.description,
                )
        finally:
            dialog.close()


if __name__ == "__main__":
    unittest.main()
