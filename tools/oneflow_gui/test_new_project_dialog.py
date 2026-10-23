"""Regression tests for the new project template selector."""

from __future__ import annotations

import os
import unittest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication

from tools.oneflow_gui.new_project_dialog import NewProjectDialog
from tools.oneflow_gui.project_templates import available_project_templates


class NewProjectDialogTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.application = QApplication.instance() or QApplication([])

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
