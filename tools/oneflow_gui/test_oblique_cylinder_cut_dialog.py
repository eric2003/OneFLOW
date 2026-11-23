"""Tests for the blind-STL oblique cylinder cutting GUI handoff."""
from __future__ import annotations
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from PySide6.QtWidgets import QApplication, QPushButton
from tools.oneflow_gui.oblique_cylinder_cut_dialog import ObliqueCylinderCutDialog

class ObliqueCylinderCutDialogTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.application = QApplication.instance() or QApplication([])
        cls.repository_root = Path(__file__).resolve().parents[2]

    def test_conversion_and_oneflow_action_is_present_but_disabled_until_mesh_exists(self):
        dialog = ObliqueCylinderCutDialog(self.repository_root)
        try:
            button = next(item for item in dialog.findChildren(QPushButton)
                          if item.text() == "Convert + Run OneFLOW...")
            self.assertFalse(button.isEnabled())
            self.assertEqual(dialog._mesh_path().name, "cylinder_oblique_cut.msh")
        finally:
            dialog.close()

    def test_existing_mesh_is_handed_to_mesh_to_cfd_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            mesh = output / "results" / "cylinder_oblique_cut.msh"
            mesh.parent.mkdir()
            mesh.write_text("$MeshFormat\\n", encoding="ascii")
            dialog = ObliqueCylinderCutDialog(self.repository_root)
            try:
                dialog.output_directory.setText(str(output))
                with patch("tools.oneflow_gui.mesh_validation_dialog.MeshValidationDialog") as validation:
                    dialog._run_mesh_to_cfd()
                validation.assert_called_once_with(
                    self.repository_root.resolve(), mesh_file=mesh.resolve(),
                    project_directory=output.resolve(),
                    mesh_case_id="stl_to_gmsh_oblique_cylinder_cut", parent=dialog)
                validation.return_value.exec.assert_called_once()
            finally:
                dialog.close()

if __name__ == "__main__":
    unittest.main()
