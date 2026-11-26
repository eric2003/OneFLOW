"""UI layout and real item-view editor lifecycle tests for the development case selector."""

from __future__ import annotations

import os
import subprocess
import sys
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
if os.name == "nt":
    # PySide6 wheels no longer bundle Qt fonts. Point Qt at the fonts already
    # installed by Windows so the expected font-directory warning does not
    # become a process-aborting error under QT_FATAL_WARNINGS=1.
    windows_fonts = Path(os.environ.get("WINDIR", "C:/Windows")) / "Fonts"
    if windows_fonts.is_dir():
        os.environ.setdefault("QT_QPA_FONTDIR", str(windows_fonts))
# Turn Qt ownership/lifetime warnings into immediate test failures on both OSes.
os.environ.setdefault("QT_FATAL_WARNINGS", "1")

from PySide6.QtCore import QCoreApplication, QEvent, Qt, qInstallMessageHandler
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QComboBox, QHeaderView

from .models import CaseDefinition, CaseImplementation
from .test_dialog import CaseDatabaseSelectionDialog


class CaseDatabaseSelectionDialogLayoutTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if os.environ.get("ONEFLOW_QT_TEST_CHILD") == "1":
            import faulthandler

            faulthandler.enable()
            os.write(2, b"[qt-selector] setUpClass: before QApplication\n")

            def report_qt_message(mode, context, message):
                details = (
                    f"[qt-selector] Qt message mode={mode}: {message}\n"
                )
                os.write(2, details.encode("utf-8", errors="replace"))

            qInstallMessageHandler(report_qt_message)
        cls._application = QApplication.instance() or QApplication([])
        if os.environ.get("ONEFLOW_QT_TEST_CHILD") == "1":
            os.write(2, b"[qt-selector] setUpClass: QApplication ready\n")

    def _run_isolated_test(self, test_name: str) -> None:
        """Run one Qt test in a fresh process so native crashes are attributable."""
        environment = os.environ.copy()
        environment["ONEFLOW_QT_TEST_CHILD"] = "1"
        result = subprocess.run(
            [sys.executable, "-u", "-X", "faulthandler", "-m", "unittest", f"{__name__}.{self.__class__.__name__}.{test_name}"],
            cwd=Path(__file__).resolve().parents[2],
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        self.assertEqual(
            result.returncode,
            0,
            f"Isolated Qt test {test_name} failed with exit code "
            f"{result.returncode}:\n{result.stdout}",
        )

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
        if os.environ.get("ONEFLOW_QT_TEST_CHILD") != "1":
            self._run_isolated_test("test_checkbox_is_in_case_column_and_columns_are_resizable")
            return
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
        os.write(2, b"[qt-selector] layout: before dialog construction\n")
        dialog = CaseDatabaseSelectionDialog([case], [case])
        os.write(2, b"[qt-selector] layout: dialog constructed\n")
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
            os.write(2, b"[qt-selector] layout: before deferred-delete flush\n")
            self._flush_deferred_deletes()
            os.write(2, b"[qt-selector] layout: deferred-delete flush complete\n")

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
        if os.environ.get("ONEFLOW_QT_TEST_CHILD") != "1":
            self._run_isolated_test("test_delegate_editor_lifecycle_is_stable_across_repeated_dialogs")
            return
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
        if os.environ.get("ONEFLOW_QT_TEST_CHILD") != "1":
            self._run_isolated_test("test_multiple_implementations_are_selectable_per_case")
            return
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
