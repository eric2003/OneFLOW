"""UI layout tests for the development case selector."""

from __future__ import annotations

import os
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication, QHeaderView

from .models import CaseDefinition
from .test_dialog import CaseDatabaseSelectionDialog


class CaseDatabaseSelectionDialogLayoutTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls._application = QApplication.instance() or QApplication([])

    def test_tree_column_is_separate_from_run_checkbox_column(self) -> None:
        case = CaseDefinition(
            case_id="sample",
            title="Sample development case",
            case_directory=Path("."),
            runtime_config=None,
            profile=None,
            layout=None,
            group="Mesh Workflow",
            dimension=3,
        )
        dialog = CaseDatabaseSelectionDialog([case], [case])
        try:
            tree = dialog.case_tree
            header = tree.header()
            self.assertEqual(tree.treePosition(), 1)
            self.assertEqual(tree.columnWidth(0), 64)
            self.assertEqual(tree.columnWidth(1), 690)
            self.assertEqual(tree.columnWidth(2), 330)
            for column in range(3):
                self.assertEqual(
                    header.sectionResizeMode(column),
                    QHeaderView.ResizeMode.Interactive,
                )
            dimension = tree.topLevelItem(0)
            self.assertEqual(dimension.text(1), "3D")
            category = dimension.child(0)
            item = category.child(0)
            self.assertEqual(item.text(1), "Sample development case")
            self.assertEqual(item.checkState(0).value, 2)
        finally:
            dialog.close()


if __name__ == "__main__":
    unittest.main()
