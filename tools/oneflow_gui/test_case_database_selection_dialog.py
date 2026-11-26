"""UI layout tests for the development case selector."""

from __future__ import annotations

import os
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication, QComboBox, QHeaderView

from .models import CaseDefinition, CaseImplementation
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
            self.assertEqual(tree.columnCount(), 3)
            self.assertEqual(tree.treePosition(), 0)
            self.assertEqual(tree.headerItem().text(0), "Case / Group")
            self.assertEqual(tree.headerItem().text(1), "Configuration")
            self.assertEqual(tree.headerItem().text(2), "Implementation")
            self.assertEqual(tree.columnWidth(0), 680)
            self.assertEqual(tree.columnWidth(1), 260)
            self.assertEqual(tree.columnWidth(2), 250)
            for column in range(3):
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



    def test_multiple_implementations_are_selectable_per_case(self) -> None:
        case = CaseDefinition(
            case_id="pipe-bend",
            title="90-Degree Pipe Bend",
            case_directory=Path("."),
            runtime_config=None,
            profile=None,
            layout=None,
            execution_kind="python",
            implementations=(
                CaseImplementation("python_reference", "python", Path("run_python.py")),
                CaseImplementation("native_cpp_geometry", "cpp", Path("run_cpp.py"), "c++17"),
            ),
            default_implementation="python_reference",
        )
        dialog = CaseDatabaseSelectionDialog(
            [case], [case], selected_implementations={"pipe-bend": "native_cpp_geometry"}
        )
        try:
            combo = dialog.findChild(QComboBox)
            self.assertIsNotNone(combo)
            self.assertEqual(combo.currentData(), "native_cpp_geometry")
            self.assertEqual(
                dialog.selected_implementations()["pipe-bend"], "native_cpp_geometry"
            )
            combo.setCurrentIndex(0)
            self.assertEqual(
                dialog.selected_implementations()["pipe-bend"], "python_reference"
            )
        finally:
            dialog.close()


if __name__ == "__main__":
    unittest.main()
