"""Tests for the unified Contract Workbench domain navigation."""

from __future__ import annotations

import os
import unittest
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication

from tools.oneflow_gui.contract_workbench_dialog import (
    ContractDomainDialog,
    ContractWorkbenchDialog,
)


class ContractWorkbenchDialogTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.application = QApplication.instance() or QApplication([])
        cls.repository_root = Path(__file__).resolve().parents[2]

    def test_domain_hub_exposes_current_contract_domains(self) -> None:
        dialog = ContractWorkbenchDialog(self.repository_root)
        try:
            self.assertEqual(
                set(dialog.domain_buttons),
                {"formula", "geometry", "case", "runtime"},
            )
            self.assertEqual(
                dialog.domain_buttons["formula"].text(),
                "Formula & Knowledge",
            )
            self.assertEqual(dialog.domain_buttons["case"].text(), "Case")
            self.assertEqual(dialog.domain_buttons["runtime"].text(), "Runtime")
        finally:
            dialog.close()

    def test_case_domain_displays_typed_case_setup_tab(self) -> None:
        dialog = ContractDomainDialog(self.repository_root, "case")
        try:
            self.assertGreaterEqual(dialog.contract_list.count(), 1)
            self.assertEqual(dialog.tabs.tabText(dialog.tabs.count() - 1), "Case Setup")
            self.assertIn("SLAU2", dialog.case_view.toPlainText())
            self.assertIn("Blasius", dialog.case_view.toPlainText())
        finally:
            dialog.close()

    def test_runtime_domain_uses_shared_contract_navigation(self) -> None:
        dialog = ContractDomainDialog(self.repository_root, "runtime")
        try:
            self.assertGreaterEqual(dialog.contract_list.count(), 1)
            tab_names = [
                dialog.tabs.tabText(index)
                for index in range(dialog.tabs.count())
            ]
            self.assertIn("Overview", tab_names)
            self.assertIn("Principle", tab_names)
            self.assertIn("Architecture", tab_names)
            self.assertIn("Implementation / Code", tab_names)
            self.assertIn("Verification", tab_names)
            self.assertNotIn("Case Setup", tab_names)
        finally:
            dialog.close()


if __name__ == "__main__":
    unittest.main()
