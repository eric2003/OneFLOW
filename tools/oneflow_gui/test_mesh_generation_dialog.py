"""Tests for the Mesh Generation Workbench navigation and STL command builder."""

from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication, QPushButton
from tools.oneflow_gui.application_settings import ApplicationSettings
from tools.oneflow_gui.user_data_paths import UserDataPaths

from tools.oneflow_gui.mesh_generation_dialog import (
    MeshGenerationDialog,
    StlToGmshDialog,
)


class MeshGenerationDialogTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.application = QApplication.instance() or QApplication([])
        cls.repository_root = Path(__file__).resolve().parents[2]

    def test_home_window_exposes_mesh_generation_as_single_entry(self) -> None:
        from tools.oneflow_gui.home_window import WorkbenchHomeWindow

        settings = ApplicationSettings(
            repository_root=self.repository_root,
            build_configuration="Debug",
            executable_path="",
            run_mode="single_case",
        )
        with patch(
            "tools.oneflow_gui.home_window.ApplicationSettings.load",
            return_value=settings,
        ), patch("tools.oneflow_gui.home_window.QTimer.singleShot"):
            window = WorkbenchHomeWindow(self.repository_root)
        try:
            buttons = window.findChildren(QPushButton)
            mesh_buttons = [
                button.text() for button in buttons
                if "Mesh" in button.text()
            ]
            self.assertIn("Mesh Generation...", mesh_buttons)
            self.assertNotIn("Mesh Layout Editor...", mesh_buttons)
            with patch(
                "tools.oneflow_gui.home_window.MeshGenerationDialog"
            ) as dialog_type:
                window._open_mesh_generation()
            dialog_type.assert_called_once_with(
                self.repository_root.resolve(), window
            )
            dialog_type.return_value.exec.assert_called_once()
        finally:
            window.close()

    def test_mesh_generation_hub_exposes_layout_and_stl_tools(self) -> None:
        dialog = MeshGenerationDialog(self.repository_root)
        try:
            buttons = dialog.findChildren(QPushButton)
            labels = {button.text() for button in buttons}
            self.assertIn("Mesh Layout Editor", labels)
            self.assertIn("STL to Gmsh", labels)
            self.assertIn("Mesh Knowledge", labels)
            with patch(
                "tools.oneflow_gui.mesh_generation_dialog.MeshLayoutDialog"
            ) as layout_type:
                dialog._open_mesh_layout()
            layout_type.assert_called_once_with(self.repository_root.resolve(), dialog)
            layout_type.return_value.exec.assert_called_once()
            with patch(
                "tools.oneflow_gui.mesh_generation_dialog.StlToGmshDialog"
            ) as stl_type:
                dialog._open_stl_to_gmsh()
            stl_type.assert_called_once_with(self.repository_root.resolve(), dialog)
            stl_type.return_value.exec.assert_called_once()
            with patch(
                "tools.oneflow_gui.mesh_generation_dialog.ContractDomainDialog"
            ) as knowledge_type:
                dialog._open_mesh_knowledge()
            knowledge_type.assert_called_once_with(
                self.repository_root.resolve(), "mesh", dialog
            )
            knowledge_type.return_value.exec.assert_called_once()
        finally:
            dialog.close()

    def test_stl_dialog_starts_with_repository_unit_sphere_fixture(self) -> None:
        dialog = StlToGmshDialog(self.repository_root)
        try:
            self.assertEqual(
                Path(dialog.stl_path.text()),
                self.repository_root / "grid" / "stl_to_gmsh" / "3d" / "sphere"
                / "automatic_tetrahedra" / "python" / "input" / "model.stl",
            )
            self.assertTrue(Path(dialog.stl_path.text()).is_file())
            self.assertEqual(dialog.geometry_preset.currentData(), "sphere")
            self.assertEqual(
                Path(dialog.output_path.text()),
                UserDataPaths.discover().mesh_generation / "sphere.msh",
            )
        finally:
            dialog.close()

    def test_stl_dialog_builds_coarse_medium_and_fine_commands(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            stl = root / "unit_sphere.stl"
            stl.write_text("solid sphere\nendsolid sphere\n", encoding="ascii")
            output = root / "sphere.msh"
            dialog = StlToGmshDialog(self.repository_root)
            try:
                dialog.stl_path.setText(str(stl))
                dialog.output_path.setText(str(output))
                expected = {
                    "coarse": "0.40",
                    "medium": "0.20",
                    "fine": "0.10",
                }
                for preset, size in expected.items():
                    dialog.density.setCurrentIndex(
                        dialog.density.findData(preset)
                    )
                    command = dialog._command_text()
                    self.assertIn("--mesh-size", command)
                    self.assertIn(size, command)
                    self.assertIn("generate_volume_mesh.py", command)
                    self.assertIn(str(stl), command)
                    self.assertIn(str(output), command)
            finally:
                dialog.close()

    def test_stl_dialog_builds_custom_size_range_and_external_gmsh_command(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            stl = root / "surface.stl"
            stl.write_text("solid surface\nendsolid surface\n", encoding="ascii")
            dialog = StlToGmshDialog(self.repository_root)
            try:
                dialog.stl_path.setText(str(stl))
                dialog.output_path.setText(str(root / "surface.msh"))
                dialog.density.setCurrentIndex(dialog.density.findData("custom"))
                self.assertEqual(dialog.minimum_size.text(), "0.10")
                self.assertEqual(dialog.maximum_size.text(), "0.20")
                dialog.minimum_size.setText("0.05")
                dialog.maximum_size.setText("0.15")
                dialog.gmsh_path.setText(str(root / "gmsh"))
                command = dialog._command_text()
                self.assertIn("--mesh-size-min", command)
                self.assertIn("0.05", command)
                self.assertIn("--mesh-size-max", command)
                self.assertIn("0.15", command)
                self.assertIn("--gmsh", command)
                self.assertIn("fluid_domain", command)
                self.assertIn("wall", command)
            finally:
                dialog.close()

    def test_stl_dialog_rejects_invalid_custom_size_bounds(self) -> None:
        dialog = StlToGmshDialog(self.repository_root)
        try:
            dialog.stl_path.setText("surface.stl")
            dialog.density.setCurrentIndex(dialog.density.findData("custom"))
            dialog.minimum_size.setText("0.20")
            dialog.maximum_size.setText("0.10")
            with self.assertRaisesRegex(ValueError, "minimum must not exceed maximum"):
                dialog._arguments()
        finally:
            dialog.close()


if __name__ == "__main__":
    unittest.main()
