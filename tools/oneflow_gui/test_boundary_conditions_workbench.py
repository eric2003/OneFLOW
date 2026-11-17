"""Tests for the dedicated Boundary Conditions Workbench."""
from __future__ import annotations
import os
import unittest
from pathlib import Path
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from PySide6.QtWidgets import QApplication
from tools.oneflow_gui.boundary_conditions_workbench import BoundaryConditionsWorkbenchDialog, load_boundary_registry

class BoundaryConditionsWorkbenchTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])
        cls.root = Path(__file__).resolve().parents[2]

    def test_registry_has_unique_semantic_ids(self):
        data = load_boundary_registry(self.root)
        ids = [x["id"] for x in data["conditions"]]
        self.assertEqual(len(ids), len(set(ids)))
        self.assertEqual(len(ids), 5)

    def test_gui_exposes_domain_views_and_search(self):
        dialog = BoundaryConditionsWorkbenchDialog(self.root)
        try:
            self.assertEqual([dialog.tabs.tabText(i) for i in range(dialog.tabs.count())],
                ["Overview", "Formulation", "Parameters", "Mappings", "Implementation", "Verification", "References"])
            self.assertEqual(dialog.condition_list.count(), 5)
            self.assertIn("q(-1)=q(1)", dialog.views["Formulation"].toPlainText())
            dialog.search_edit.setText("2003")
            self.assertEqual(dialog.condition_list.count(), 1)
            self.assertIn("flow.inflow.total_conditions", dialog.condition_list.item(0).data(256))
        finally:
            dialog.close()

if __name__ == "__main__":
    unittest.main()
