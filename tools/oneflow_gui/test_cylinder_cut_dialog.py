"""Tests for the analytic cylinder cutting and CFD GUI."""
from __future__ import annotations
import os, tempfile, unittest
from pathlib import Path
from unittest.mock import patch
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from PySide6.QtWidgets import QApplication, QPushButton
from tools.oneflow_gui.cylinder_cut_dialog import CYLINDER_MESH_CASE_ID, CylinderCutWorkbenchDialog
from tools.oneflow_gui.mesh_generation_dialog import MeshGenerationDialog
from tools.oneflow_gui.user_data_paths import UserDataPaths

class CylinderCutWorkbenchDialogTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.application = QApplication.instance() or QApplication([])
        cls.repository_root = Path(__file__).resolve().parents[2]

    def test_mesh_generation_has_first_class_cylinder_cut_entry(self):
        dialog = MeshGenerationDialog(self.repository_root)
        try:
            button = next(item for item in dialog.findChildren(QPushButton) if item.text() == "Cylinder Cutting & CFD...")
            with patch("tools.oneflow_gui.mesh_generation_dialog.CylinderCutWorkbenchDialog") as workbench:
                dialog._open_cylinder_cut_workbench()
            workbench.assert_called_once_with(self.repository_root.resolve(), dialog)
            workbench.return_value.exec.assert_called_once()
            self.assertIsNotNone(button)
        finally:
            dialog.close()

    def test_default_command_contains_cut_and_mesh_parameters(self):
        with tempfile.TemporaryDirectory() as temporary:
            dialog = CylinderCutWorkbenchDialog(self.repository_root)
            try:
                command = dialog._command_text()
                self.assertIn("cylinder_volume_mesh.py", command)
                self.assertIn("--cut-start -1.5", command)
                self.assertIn("--cut-end 1.5", command)
                self.assertIn("--mesh-size 0.2", command)
                self.assertIn(str(UserDataPaths.discover().mesh_generation), command)
                dialog.output_directory.setText(str(Path(temporary) / "output"))
                self.assertFalse(dialog.run_oneflow_button.isEnabled())
            finally:
                dialog.close()

    def test_run_oneflow_hands_current_mesh_to_existing_mesh_validation(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            mesh = output / "cylinder_fluid_domain.msh"
            mesh.write_text("$MeshFormat\\n", encoding="ascii")
            dialog = CylinderCutWorkbenchDialog(self.repository_root)
            try:
                dialog.output_directory.setText(str(output))
                with patch("tools.oneflow_gui.mesh_validation_dialog.MeshValidationDialog") as validation:
                    dialog._run_oneflow()
                validation.assert_called_once_with(self.repository_root.resolve(), mesh_file=mesh.resolve(),
                    project_directory=output.resolve(), mesh_case_id=CYLINDER_MESH_CASE_ID, parent=dialog)
                validation.return_value.exec.assert_called_once()
            finally:
                dialog.close()

if __name__ == "__main__":
    unittest.main()
