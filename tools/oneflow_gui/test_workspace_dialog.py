"""Tests for workspace case path editing and directory selection."""

from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication

from tools.oneflow_gui.workspace_dialog import WorkspaceDialog


class WorkspaceDialogTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls._application = QApplication.instance() or QApplication([])

    def test_browse_updates_case_path_and_uses_current_directory(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            current = root / "current"
            selected = root / "selected"
            current.mkdir()
            selected.mkdir()
            dialog = WorkspaceDialog()
            dialog._append_case("case-a", str(current))

            with patch(
                "tools.oneflow_gui.workspace_dialog.QFileDialog.getExistingDirectory",
                return_value=str(selected),
            ) as choose_directory:
                dialog._browse_case_directory(0)

            self.assertEqual(dialog.case_table.item(0, 1).text(), str(selected))
            self.assertEqual(choose_directory.call_args.args[2], str(current))
            dialog.close()

    def test_cancelled_browse_preserves_typed_case_path(self) -> None:
        dialog = WorkspaceDialog()
        typed_path = "D:/cases/my-case"
        dialog._append_case("case-a", typed_path)

        with patch(
            "tools.oneflow_gui.workspace_dialog.QFileDialog.getExistingDirectory",
            return_value="",
        ):
            dialog._browse_case_directory(0)

        self.assertEqual(dialog.case_table.item(0, 1).text(), typed_path)
        dialog.close()


if __name__ == "__main__":
    unittest.main()
