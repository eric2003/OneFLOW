"""UI layout and real item-view editor lifecycle tests for the development case selector."""

from __future__ import annotations

import os
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
# Turn Qt ownership/lifetime warnings into immediate test failures on both OSes.
os.environ.setdefault("QT_FATAL_WARNINGS", "1")

from PySide6.QtCore import QCoreApplication, QEvent, Qt
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QComboBox, QHeaderView

from .models import CaseDefinition, CaseImplementation
from .test_dialog import CaseDatabaseSelectionDialog


class CaseDatabaseSelectionDialogLayoutTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls._application = QApplication.instance() or QApplication([])

    def _flush_deferred_deletes(self) -> None:
        # processEvents() alone does not guarantee delivery of DeferredDelete
        # events when tests drive Qt without entering the main event loop.
        QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
        self._application.processEvents()

    def _select_implementation_through_view(
        self, dialog: CaseDatabaseSelectionDialog, implementation_index: int
    ) -> None:
        """Exercise QTreeWidget's real delegate/editor ownership protocol."""
        tree = dialog.case_tree
        item = tree.topLevelItem(0).child(0).child(0)
        dialog.show()
        self._application.processEvents()

        tree.editItem(item, 2)
        self._application.processEvents()
        editor = tree.findChild(QComboBox)
        self.assertIsNotNone(editor, "QTreeWidget did not create its delegate editor")
        assert editor is not None
        self.assertGreaterEqual(implementation_index, 0)
        self.assertLess(implementation_index, editor.count())
        editor.setCurrentIndex(implementation_index)
        # Let QAbstractItemView commit and close the editor. Tests must not call
        # createEditor(), destroyEditor(), or deleteLater() on the editor.
        QTest.keyClick(editor, Qt.Key.Key_Return)
        self._application.processEvents()

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
            dialog.deleteLater()
            self._flush_deferred_deletes()

    def _pipe_bend_case(self) -> CaseDefinition:
        return CaseDefinition(
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

    def test_delegate_editor_lifecycle_is_stable_across_repeated_dialogs(self) -> None:
        case = self._pipe_bend_case()
        for _ in range(25):
            dialog = CaseDatabaseSelectionDialog([case], [case])
            try:
                self._select_implementation_through_view(dialog, 1)
                self.assertEqual(
                    dialog.selected_implementations()["pipe-bend"],
                    "native_cpp_geometry",
                )
            finally:
                dialog.close()
                dialog.deleteLater()
                self._flush_deferred_deletes()

    def test_multiple_implementations_are_selectable_per_case(self) -> None:
        case = self._pipe_bend_case()
        dialog = CaseDatabaseSelectionDialog(
            [case], [case], selected_implementations={"pipe-bend": "native_cpp_geometry"}
        )
        try:
            item = dialog.case_tree.topLevelItem(0).child(0).child(0)
            self.assertEqual(item.text(2), "Native C++ geometry")
            self.assertEqual(
                dialog.selected_implementations()["pipe-bend"], "native_cpp_geometry"
            )
            self._select_implementation_through_view(dialog, 0)
            self.assertEqual(item.text(2), "Python reference")
            self.assertEqual(
                dialog.selected_implementations()["pipe-bend"], "python_reference"
            )
        finally:
            dialog.close()
            dialog.deleteLater()
            self._flush_deferred_deletes()


if __name__ == "__main__":
    unittest.main()
