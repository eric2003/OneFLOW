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

    def test_checkbox_is_in_case_column_and_columns_are_resizable(self) -> None:
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
            self.assertEqual(tree.columnCount(), 2)
            self.assertEqual(tree.treePosition(), 0)
            self.assertEqual(tree.headerItem().text(0), "Case / Group")
            self.assertEqual(tree.headerItem().text(1), "Configuration")
            self.assertEqual(tree.columnWidth(0), 760)
            self.assertEqual(tree.columnWidth(1), 330)
            for column in range(2):
                self.assertEqual(
                    header.sectionResizeMode(column),
                    QHeaderView.ResizeMode.Interactive,
                )
            dimension = tree.topLevelItem(0)
            self.assertEqual(dimension.text(0), "3D")
            category = dimension.child(0)
            item = category.child(0)
            self.assertEqual(item.text(0), "Sample development case")
            self.assertEqual(item.checkState(0).value, 2)
        finally:
            dialog.close()


if __name__ == "__main__":
    unittest.main()
