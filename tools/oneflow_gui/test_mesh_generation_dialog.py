"""Tests for the Mesh Generation Workbench navigation and STL command builder."""

from __future__ import annotations

import gc
import hashlib
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import QEvent
from PySide6.QtWidgets import QApplication, QPushButton, QScrollArea


_QT_APPLICATION = QApplication.instance() or QApplication([])
from tools.oneflow_gui.application_settings import ApplicationSettings
from tools.oneflow_gui.user_data_paths import UserDataPaths
from grid.stl_to_gmsh.generate_volume_mesh import (
    DEFAULT_MESH_SIZE,
    DEFAULT_SPHERE_STL,
    default_output_path,
    discover_gmsh,
)

from tools.oneflow_gui import mesh_validation_dialog
from tools.oneflow_gui.mesh_generation_dialog import (
    MeshGenerationDialog,
    StlGeometryDialog,
    StlToGmshDialog,
)
from tools.oneflow_gui.generic_surface_cut_dialog import GenericSurfaceCutDialog
from tools.oneflow_gui.case_repository import select_cases_for_selection_dialog
from tools.oneflow_gui.mesh_project import (
    create_mesh_project_from_stl,
    create_mesh_project_from_template,
    load_mesh_project,
)


class MeshGenerationDialogTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.application = _QT_APPLICATION
        cls.repository_root = Path(__file__).resolve().parents[2]

    def tearDown(self) -> None:
        # Close and destroy every top-level widget created by this test before
        # the next test starts. Closing alone leaves Qt-owned child objects and
        # deferred-delete events queued until interpreter shutdown on Linux.
        for widget in QApplication.topLevelWidgets():
            widget.close()
            widget.deleteLater()
        self.application.sendPostedEvents(None, QEvent.Type.DeferredDelete)
        self.application.processEvents()
        gc.collect()
        self.application.sendPostedEvents(None, QEvent.Type.DeferredDelete)
        self.application.processEvents()

    def test_geometry_library_is_available_from_mesh_generation(self) -> None:
        from tools.oneflow_gui.geometry_library_dialog import GeometryLibraryDialog

        dialog = MeshGenerationDialog(self.repository_root)
        try:
            button = next(
                item for item in dialog.findChildren(QPushButton)
                if item.text() == "Geometry Library..."
            )
            self.assertIsNotNone(button)
            edited_assets = []
            library = GeometryLibraryDialog(
                self.repository_root, dialog, on_edit_asset=edited_assets.append
            )
            try:
                self.assertGreaterEqual(library.asset_list.count(), 3)
                self.assertTrue(library.import_button.isEnabled())
                self.assertTrue(library.edit_button.isEnabled())
                self.assertIn("geometry.sphere.unit", library.details.toPlainText())
                library._edit_selected_asset()
                self.assertEqual(edited_assets[0]["id"], "geometry.sphere.unit")
                self.assertIn("\nDefinition\n", library.details.toPlainText())
                self.assertNotIn("\\nDefinition", library.details.toPlainText())
                displayed = library.details.toPlainText()
                self.assertIn("\nProvenance\n", displayed)
                self.assertIn("\nRepresentations\n", displayed)
                self.assertIn("\n  Kind: ", displayed)
                self.assertIn("\n\nDefinition\n", displayed)
                self.assertNotIn("\\n", displayed)
                self.assertIn("analytic-equation", library.details.toPlainText())
                library.asset_list.setCurrentRow(1)
                self.assertEqual(edited_assets[0]["id"], "geometry.sphere.unit")
                library._edit_selected_asset()
                self.assertEqual(edited_assets[1]["id"], "geometry.spheroid.prolate_6_to_1")
                self.assertIn("geometry.spheroid.prolate_6_to_1", library.details.toPlainText())
                self.assertIn("generate_stl_geometry.py", library.details.toPlainText())
            finally:
                library.close()
            with patch(
                "tools.oneflow_gui.mesh_generation_dialog.GeometryLibraryDialog"
            ) as library_type:
                dialog._open_geometry_library()
            library_type.assert_called_once()
            self.assertEqual(
                library_type.call_args.args,
                (self.repository_root.resolve(), dialog),
            )
            self.assertTrue(callable(library_type.call_args.kwargs["on_edit_asset"]))
            library_type.return_value.exec.assert_called_once()
        finally:
            dialog.close()

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
            self.assertIn("Mesh-to-CFD Validation...", mesh_buttons)
            self.assertNotIn("Mesh Layout Editor...", mesh_buttons)
            with patch(
                "tools.oneflow_gui.home_window.MeshGenerationDialog"
            ) as dialog_type:
                window._open_mesh_generation()
            dialog_type.assert_called_once_with(
                self.repository_root.resolve(), window
            )
            dialog_type.return_value.exec.assert_called_once()
            with patch(
                "tools.oneflow_gui.home_window.MeshValidationDialog"
            ) as validation_type:
                window._open_mesh_validation()
            validation_type.assert_called_once_with(
                self.repository_root.resolve(), parent=window
            )
            validation_type.return_value.exec.assert_called_once()
        finally:
            window.close()

    def test_mesh_generation_hub_exposes_layout_and_stl_tools(self) -> None:
        dialog = MeshGenerationDialog(self.repository_root)
        try:
            buttons = dialog.findChildren(QPushButton)
            labels = {button.text() for button in buttons}
            self.assertIn("Mesh Layout Editor", labels)
            self.assertIn("STL to Gmsh", labels)
            self.assertIn("Mesh-to-CFD Validation...", labels)
            self.assertIn("Geometry Library...", labels)
            self.assertNotIn("Generate STL Geometry", labels)
            self.assertIn("Mesh Knowledge", labels)
            self.assertIn("New from Template", labels)
            self.assertIn("Open Existing Project...", labels)
            self.assertIn("Import STL / Grid Database...", labels)
            self.assertIn(
                "stl_to_gmsh_python_prolate_spheroid_automatic_tetrahedra",
                {dialog.project_template.itemData(i) for i in range(dialog.project_template.count())},
            )
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
            with patch(
                "tools.oneflow_gui.mesh_generation_dialog.MeshValidationDialog"
            ) as validation_type:
                dialog._open_mesh_validation()
            validation_type.assert_called_once_with(
                self.repository_root.resolve(), parent=dialog
            )
            validation_type.return_value.exec.assert_called_once()
        finally:
            dialog.close()


    def test_geometry_library_asset_opens_matching_parameter_editor(self) -> None:
        dialog = StlGeometryDialog(
            self.repository_root, asset_id="geometry.cylinder.internal_flow"
        )
        try:
            self.assertEqual(dialog.geometry.currentData(), "cylinder_internal_flow")
            self.assertFalse(dialog.cylinder_group.isHidden())
            self.assertTrue(dialog.sphere_group.isHidden())
        finally:
            dialog.close()

    def test_generate_stl_geometry_dialog_creates_cylinder_and_manifest(self) -> None:
        from types import SimpleNamespace

        with tempfile.TemporaryDirectory() as directory:
            workspace = Path(directory) / "mesh-generation"
            with patch(
                "tools.oneflow_gui.mesh_generation_dialog.UserDataPaths.discover",
                return_value=SimpleNamespace(mesh_generation=workspace),
            ):
                dialog = StlGeometryDialog(self.repository_root)
                try:
                    self.assertEqual(dialog.geometry.currentData(), "sphere")
                    dialog.geometry.setCurrentIndex(
                        dialog.geometry.findData("cylinder_internal_flow")
                    )
                    self.assertTrue(dialog.sphere_group.isHidden())
                    self.assertTrue(dialog.spheroid_group.isHidden())
                    self.assertFalse(dialog.cylinder_group.isHidden())
                    self.assertEqual(dialog.cylinder_radius.text(), "1.0")
                    self.assertIn("not solid walls", dialog.patch_note.text())
                    output = Path(directory) / "cylinder.stl"
                    dialog.output.setText(str(output))
                    dialog._generate()
                    self.assertTrue(output.is_file())
                    manifest = output.with_suffix(".geometry.json")
                    self.assertTrue(manifest.is_file())
                    self.assertEqual(dialog.generated_stl, output.resolve())
                    self.assertIn("inflow", dialog.result_view.toPlainText())
                    self.assertIn("outflow", dialog.result_view.toPlainText())
                    self.assertIn("Manifest:", dialog.result_view.toPlainText())
                    self.assertIn("Patch triangles:", dialog.result_view.toPlainText())
                    geometry_manifest = json.loads(manifest.read_text(encoding="utf-8"))
                    self.assertEqual(
                        set(geometry_manifest["patch_triangle_indices"]),
                        {"wall", "inflow", "outflow"},
                    )
                    project = create_mesh_project_from_stl(
                        self.repository_root, output, "cylinder-ui"
                    )
                    project_geometry_manifest = (
                        project["stl_path"].with_suffix(".geometry.json")
                    )
                    self.assertTrue(project_geometry_manifest.is_file())
                    project_geometry = json.loads(
                        project_geometry_manifest.read_text(encoding="utf-8")
                    )
                    self.assertEqual(
                        project_geometry["asset_id"], "geometry.cylinder.internal_flow"
                    )
                    self.assertEqual(
                        project["geometry_lineage"]["provenance"]["kind"], "generated"
                    )
                    self.assertIsNone(project["boundary_regions"])
                    gmsh_dialog = StlToGmshDialog(
                        self.repository_root, project=project
                    )
                    try:
                        self.assertEqual(
                            gmsh_dialog.geometry_preset.currentData(),
                            "cylinder_internal_flow",
                        )
                        self.assertTrue(gmsh_dialog.geometry_preset.isEnabled())
                        self.assertFalse(gmsh_dialog.boundary_name.isEnabled())
                        self.assertFalse(gmsh_dialog.volume_name.isEnabled())
                        self.assertEqual(gmsh_dialog.stl_path.text(), str(project["stl_path"]))
                        self.assertEqual(
                            gmsh_dialog._mesh_case_id,
                            "stl_to_gmsh_existing_mesh_3d_cylinder_internal_flow",
                        )
                        self.assertEqual(
                            gmsh_dialog.validate_button.text(),
                            "Mesh-to-CFD Validation...",
                        )
                        lineage_text = gmsh_dialog.geometry_lineage_view.toPlainText()
                        self.assertIn(
                            "Boundary patches: inflow, outflow, wall",
                            lineage_text.splitlines(),
                        )
                        self.assertIn("\nGeometry asset:", lineage_text)
                        self.assertNotIn("\\n", lineage_text)
                    finally:
                        gmsh_dialog.close()
                finally:
                    dialog.close()

    def test_generated_cylinder_can_open_region_aware_stl_to_gmsh(self) -> None:
        from types import SimpleNamespace

        with tempfile.TemporaryDirectory() as directory:
            workspace = Path(directory) / "mesh-generation"
            with patch(
                "tools.oneflow_gui.mesh_generation_dialog.UserDataPaths.discover",
                return_value=SimpleNamespace(mesh_generation=workspace),
            ):
                dialog = StlGeometryDialog(self.repository_root)
                try:
                    dialog.geometry.setCurrentIndex(
                        dialog.geometry.findData("cylinder_internal_flow")
                    )
                    output = Path(directory) / "cylinder.stl"
                    dialog.output.setText(str(output))
                    dialog._generate()
                    self.assertTrue(dialog.open_gmsh_button.isEnabled())
                    with patch(
                        "tools.oneflow_gui.mesh_generation_dialog.create_mesh_project_from_stl",
                        return_value={
                            "directory": workspace / "cylinder",
                            "name": "cylinder",
                            "geometry": "cylinder_internal_flow",
                        },
                    ) as create_project, patch(
                        "tools.oneflow_gui.mesh_generation_dialog.StlToGmshDialog"
                    ) as stl_dialog:
                        dialog._open_in_gmsh()
                    create_project.assert_called_once_with(
                        self.repository_root.resolve(), output.resolve()
                    )
                    stl_dialog.assert_called_once()
                    self.assertEqual(
                        stl_dialog.call_args.kwargs["project"]["geometry"],
                        "cylinder_internal_flow",
                    )
                    stl_dialog.return_value.exec.assert_called_once()
                finally:
                    dialog.close()

    def test_geometry_dialog_settings_are_scrollable_on_short_screens(self) -> None:
        dialog = StlGeometryDialog(self.repository_root)
        try:
            self.assertEqual(dialog.minimumHeight(), 480)
            scroll = dialog.findChild(QScrollArea, "geometrySettingsScrollArea")
            self.assertIsNotNone(scroll)
            self.assertTrue(scroll.widgetResizable())
            self.assertIs(scroll.widget().layout().parentWidget(), scroll.widget())
            self.assertIs(dialog.sphere_group.parentWidget(), scroll.widget())
            self.assertIs(dialog.spheroid_group.parentWidget(), scroll.widget())
            self.assertIs(dialog.cylinder_group.parentWidget(), scroll.widget())
        finally:
            dialog.close()

    def test_geometry_parameter_panels_are_independent_and_retain_values(self) -> None:
        dialog = StlGeometryDialog(self.repository_root)
        try:
            self.assertFalse(dialog.sphere_group.isHidden())
            self.assertTrue(dialog.spheroid_group.isHidden())
            self.assertTrue(dialog.cylinder_group.isHidden())

            dialog.radius.setText("2.5")
            dialog.subdivisions.setText("4")
            dialog.geometry.setCurrentIndex(
                dialog.geometry.findData("cylinder_internal_flow")
            )
            self.assertTrue(dialog.sphere_group.isHidden())
            self.assertFalse(dialog.cylinder_group.isHidden())
            dialog.cylinder_radius.setText("0.75")
            dialog.length.setText("6.0")
            dialog.circumferential.setText("64")
            dialog.axial.setText("20")

            dialog.geometry.setCurrentIndex(
                dialog.geometry.findData("prolate_spheroid")
            )
            self.assertFalse(dialog.spheroid_group.isHidden())
            self.assertTrue(dialog.cylinder_group.isHidden())
            dialog.semi_axis_a.setText("4.0")
            dialog.semi_axis_b.setText("0.8")
            dialog.semi_axis_c.setText("0.6")
            dialog.spheroid_subdivisions.setText("5")

            dialog.geometry.setCurrentIndex(dialog.geometry.findData("sphere"))
            self.assertEqual(dialog.radius.text(), "2.5")
            self.assertEqual(dialog.subdivisions.text(), "4")
            self.assertEqual(dialog.cylinder_radius.text(), "0.75")
            self.assertEqual(dialog.length.text(), "6.0")
            self.assertEqual(dialog.circumferential.text(), "64")
            self.assertEqual(dialog.axial.text(), "20")
            self.assertEqual(dialog.semi_axis_a.text(), "4.0")
            self.assertEqual(dialog.semi_axis_b.text(), "0.8")
            self.assertEqual(dialog.semi_axis_c.text(), "0.6")
            self.assertEqual(dialog.spheroid_subdivisions.text(), "5")
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
                UserDataPaths.discover().mesh_generation / "sphere" / "sphere.msh",
            )
            self.assertTrue(dialog.advanced_settings.isHidden())
            self.assertEqual(dialog.gmsh_path.text(), discover_gmsh() or "")
            self.assertEqual(dialog.mesh_size.text(), "0.20")
            self.assertTrue(dialog.run_button.isEnabled())
            self.assertEqual(dialog.validate_button.isEnabled(),
                             Path(dialog.output_path.text()).is_file())
        finally:
            dialog.close()

    def test_geometry_preset_switches_stl_output_and_validation_case(self) -> None:
        dialog = StlToGmshDialog(self.repository_root)
        try:
            spheroid_stl = (
                self.repository_root / "grid" / "stl_to_gmsh" / "3d"
                / "prolate_spheroid" / "automatic_tetrahedra" / "python"
                / "input" / "model.stl"
            )
            dialog.geometry_preset.setCurrentIndex(
                dialog.geometry_preset.findData("prolate_spheroid")
            )
            self.assertEqual(Path(dialog.stl_path.text()), spheroid_stl)
            self.assertTrue(spheroid_stl.is_file())
            self.assertEqual(
                Path(dialog.output_path.text()),
                UserDataPaths.discover().mesh_generation
                / "prolate_spheroid" / "prolate_spheroid.msh",
            )
            self.assertEqual(
                dialog._mesh_case_id,
                "stl_to_gmsh_prolate_spheroid_automatic_tetrahedra",
            )
            self.assertIn(str(spheroid_stl), dialog._command_text())

            dialog.geometry_preset.setCurrentIndex(
                dialog.geometry_preset.findData("sphere")
            )
            self.assertEqual(Path(dialog.stl_path.text()), DEFAULT_SPHERE_STL)
            self.assertEqual(
                Path(dialog.output_path.text()),
                UserDataPaths.discover().mesh_generation / "sphere" / "sphere.msh",
            )
            self.assertEqual(
                dialog._mesh_case_id,
                "stl_to_gmsh_sphere_automatic_tetrahedra",
            )
        finally:
            dialog.close()

    def test_invalid_cylinder_asset_rolls_back_preset_and_keeps_input_consistent(self) -> None:
        dialog = StlToGmshDialog(self.repository_root)
        try:
            original_stl = Path(dialog.stl_path.text())
            with patch(
                "tools.oneflow_gui.mesh_generation_dialog.hashlib.sha256"
            ) as sha256:
                sha256.return_value.hexdigest.return_value = "invalid-checksum"
                with patch(
                    "tools.oneflow_gui.mesh_generation_dialog.QMessageBox.warning"
                ):
                    dialog.geometry_preset.setCurrentIndex(
                        dialog.geometry_preset.findData("cylinder_internal_flow")
                    )
            self.assertEqual(dialog.geometry_preset.currentData(), "sphere")
            self.assertEqual(Path(dialog.stl_path.text()), original_stl)
            self.assertEqual(
                dialog._mesh_case_id,
                "stl_to_gmsh_sphere_automatic_tetrahedra",
            )
        finally:
            dialog.close()

    def test_cylinder_geometry_preset_uses_bundled_standard_grid_asset(self) -> None:
        import hashlib
        import json

        dialog = StlToGmshDialog(self.repository_root)
        try:
            expected_stl = (
                self.repository_root / "grid" / "stl_to_gmsh" / "3d"
                / "cylinder_internal_flow" / "automatic_tetrahedra"
                / "cpp" / "input" / "model.stl"
            )
            dialog.geometry_preset.setCurrentIndex(
                dialog.geometry_preset.findData("cylinder_internal_flow")
            )
            self.assertEqual(Path(dialog.stl_path.text()), expected_stl)
            self.assertTrue(expected_stl.is_file())
            self.assertNotIn(".oneflow", str(expected_stl).lower())

            manifest_path = expected_stl.with_suffix(".geometry.json")
            self.assertTrue(manifest_path.is_file())
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(manifest["asset_id"], "geometry.cylinder.internal_flow")
            self.assertEqual(manifest["geometry"], "cylinder_internal_flow")
            self.assertEqual(
                manifest["stl_sha256"],
                hashlib.sha256(expected_stl.read_bytes()).hexdigest(),
            )
            self.assertTrue(manifest["surface"]["watertight"])
            self.assertEqual(
                set(manifest["patch_triangle_indices"]),
                {"wall", "inflow", "outflow"},
            )
            self.assertEqual(
                manifest["parameters"],
                {
                    "radius": 1,
                    "length": 4,
                    "circumferential_divisions": 48,
                    "axial_divisions": 16,
                },
            )
            self.assertEqual(
                dialog._mesh_case_id,
                "stl_to_gmsh_cylinder_internal_flow_automatic_tetrahedra",
            )
            self.assertEqual(
                Path(dialog.output_path.text()),
                UserDataPaths.discover().mesh_generation
                / "cylinder_internal_flow" / "cylinder_internal_flow.msh",
            )
            self.assertIn(str(expected_stl), dialog._command_text())
        finally:
            dialog.close()

    def test_mesh_validation_case_tree_uses_single_selection_column(self) -> None:
        from PySide6.QtCore import Qt
        from tools.oneflow_gui.mesh_validation_dialog import MeshValidationDialog

        dialog = MeshValidationDialog(self.repository_root)
        try:
            self.assertEqual(dialog.case_tree.columnCount(), 1)
            self.assertEqual(
                dialog.case_tree.headerItem().text(0), "STL / Mesh Case"
            )
            self.assertTrue(dialog.case_tree.isColumnHidden(0) is False)
            self.assertTrue(dialog._case_items)
            self.assertIn(
                "stl_to_gmsh_pipe_bend_90_degree",
                dialog._case_items,
            )

            case_item = next(iter(dialog._case_items.values()))
            category_item = case_item.parent()
            dimension_item = category_item.parent()
            case_item.setCheckState(0, Qt.CheckState.Unchecked)
            self.assertEqual(
                category_item.checkState(0), Qt.CheckState.PartiallyChecked
            )
            self.assertEqual(
                dimension_item.checkState(0), Qt.CheckState.PartiallyChecked
            )

            category_item.setCheckState(0, Qt.CheckState.Checked)
            for index in range(category_item.childCount()):
                self.assertEqual(
                    category_item.child(index).checkState(0), Qt.CheckState.Checked
                )
        finally:
            dialog.close()

    def test_mesh_validation_dialog_reconstructs_report_views(self) -> None:
        from tools.oneflow_gui.mesh_validation_dialog import MeshValidationDialog

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            console_path = root / "mesh-validation-console.log"
            console_path.write_text("restored console output\n", encoding="utf-8")
            report_path = root / "mesh-validation-report.json"
            report = {
                "schema": "oneflow.run-report",
                "version": 1,
                "workflow": "Mesh-to-CFD Validation",
                "status": "success",
                "exit_code": 0,
                "mesh_case_id": "test_mesh_case",
                "stages": {"gmsh_to_cgns": {"status": "PASS", "message": "Converted"}},
                "metrics": {"nodes": 12},
                "artifacts": {
                    "report": str(report_path),
                    "console_log": str(console_path),
                    "mesh": str(root / "mesh.msh"),
                },
                "execution": {"console_log": str(console_path)},
                "geometry_lineage": {},
                "error": None,
            }
            report_path.write_text(json.dumps(report), encoding="utf-8")
            dialog = MeshValidationDialog(self.repository_root)
            try:
                with patch(
                    "tools.oneflow_gui.mesh_validation_dialog.QFileDialog.getOpenFileName",
                    return_value=(str(report_path), "JSON files (*.json)"),
                ):
                    dialog._open_run_report()
                self.assertEqual(dialog.result_tabs.count(), 5)
                self.assertEqual(dialog.result_tabs.tabText(0), "Results")
                self.assertEqual(dialog.result_tabs.tabText(1), "Artifacts")
                self.assertEqual(dialog.result_tabs.tabText(4), "Run History")
                self.assertGreaterEqual(dialog.result_tree.topLevelItemCount(), 1)
                self.assertGreaterEqual(dialog.artifact_tree.topLevelItemCount(), 1)
                self.assertIn("restored console output", dialog.log_view.toPlainText())
                self.assertIn("test_mesh_case", dialog.raw_report_view.toPlainText())
            finally:
                dialog.close()
                dialog.deleteLater()
                self.application.sendPostedEvents(None, QEvent.Type.DeferredDelete)
                self.application.processEvents()

    def test_mesh_validation_dialog_recovers_run_directory_without_report(self) -> None:
        from tools.oneflow_gui.mesh_validation_dialog import MeshValidationDialog

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            stage_log = root / "stl_to_gmsh/results/mesh-validation.log"
            stage_log.parent.mkdir(parents=True)
            stage_log.write_text(
                '$ ["gmsh", "model.geo"]\n'
                "[stl_to_gmsh] succeeded: Recovered from stage log.\n",
                encoding="utf-8",
            )
            (root / "mesh-validation-console.log").write_text(
                "recovered console text\n", encoding="utf-8"
            )
            (root / "mesh-validation-run.json").write_text(
                json.dumps({
                    "schema": "oneflow.mesh-validation-run",
                    "version": 1,
                    "mesh_case_id": "stl_to_gmsh_test_case",
                    "mesh_implementation": "python",
                    "mesh_input_mode": "generated_mesh",
                }),
                encoding="utf-8",
            )
            dialog = MeshValidationDialog(self.repository_root)
            try:
                self.assertTrue(dialog.recover_run_button.isEnabled())
                with patch(
                    "tools.oneflow_gui.mesh_validation_dialog.QFileDialog.getExistingDirectory",
                    return_value=str(root),
                ):
                    dialog._open_run_directory()

                self.assertEqual(
                    dialog.result_tree.topLevelItem(0).text(1), "RECOVERED"
                )
                self.assertIn("stl_to_gmsh_test_case", dialog.raw_report_view.toPlainText())
                self.assertIn("recovered console text", dialog.log_view.toPlainText())
                self.assertGreaterEqual(dialog.artifact_tree.topLevelItemCount(), 1)
            finally:
                dialog.close()
                dialog.deleteLater()
                self.application.sendPostedEvents(None, QEvent.Type.DeferredDelete)
                self.application.processEvents()

    def test_mesh_validation_worker_persists_console_log_and_links_report(self) -> None:
        from tools.oneflow_gui.mesh_validation_dialog import MeshValidationWorker

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            report_path = root / "run" / "mesh-validation-report.json"
            report_path.parent.mkdir(parents=True)
            report_path.write_text(
                json.dumps({
                    "schema": "oneflow.run-report",
                    "version": 1,
                    "workflow": "Mesh-to-CFD Validation",
                    "status": "success",
                    "exit_code": 0,
                    "stages": {"gmsh_to_cgns": {"status": "PASS"}},
                    "metrics": {},
                    "artifacts": {"report": str(report_path)},
                    "error": None,
                }),
                encoding="utf-8",
            )

            class FakeProcess:
                stdout = iter([f"report: {report_path}\n", "conversion complete\n"])

                def __enter__(self):
                    return self

                def __exit__(self, *_args):
                    return False

                def wait(self):
                    return 0

            worker = MeshValidationWorker(
                [("Test Case", ("python", "-m", "fake_case"), root)],
                root / "session-logs",
            )
            with patch(
                "tools.oneflow_gui.mesh_validation_dialog.subprocess.Popen",
                return_value=FakeProcess(),
            ):
                worker.run()

            report = json.loads(report_path.read_text(encoding="utf-8"))
            console_reference = report["execution"]["console_log"]
            self.assertEqual(report["artifacts"]["console_log"], console_reference)
            console_log = Path(console_reference)
            self.assertFalse(console_log.is_absolute())
            console_log = report_path.parent / console_log
            self.assertTrue(console_log.is_file())
            log_text = console_log.read_text(encoding="utf-8")
            self.assertIn("Working directory:", log_text)
            self.assertIn("Command:", log_text)
            self.assertIn("conversion complete", log_text)
            self.assertIn("PASSED: Test Case", log_text)
            console_record = next(
                item for item in report["artifacts"]["files"]
                if item["path"] == "mesh-validation-console.log"
            )
            self.assertEqual(console_record["size_bytes"], console_log.stat().st_size)
            self.assertEqual(
                console_record["sha256"],
                hashlib.sha256(console_log.read_bytes()).hexdigest(),
            )

    def test_validate_current_mesh_opens_single_mesh_validation_mode(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            mesh = Path(directory) / "sphere.msh"
            mesh.write_text("mesh", encoding="ascii")
            dialog = StlToGmshDialog(self.repository_root)
            try:
                dialog.output_path.setText(str(mesh))
                with patch(
                    "tools.oneflow_gui.mesh_validation_dialog.MeshValidationDialog"
                ) as dialog_type:
                    dialog._validate_current_mesh()
                dialog_type.assert_called_once_with(
                    self.repository_root.resolve(),
                    mesh_file=mesh.resolve(),
                    project_directory=mesh.parent.resolve(),
                    mesh_case_id="stl_to_gmsh_sphere_automatic_tetrahedra",
                    parent=dialog,
                )
                dialog_type.return_value.exec.assert_called_once()
            finally:
                dialog.close()

    def test_generic_cut_region_artifacts_are_copied_into_mesh_project(self) -> None:
        import hashlib
        import json
        from types import SimpleNamespace
        from tools.oneflow_gui.mesh_project import create_mesh_project_from_stl

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "cut"
            source.mkdir()
            closed = source / "closed_surface.stl"
            closed.write_text("solid closed\nendsolid closed\n", encoding="ascii")
            regions = {}
            for key, filename in (
                ("wall", "wall.stl"),
                ("inflow", "inlet.stl"),
                ("outflow", "outlet.stl"),
            ):
                path = source / filename
                path.write_text(f"solid {key}\nendsolid {key}\n", encoding="ascii")
                regions[key] = {
                    "file": filename,
                    "triangle_count": 1,
                    "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                }
            (source / "region_manifest.json").write_text(
                json.dumps({
                    "schema": "oneflow.stl-surface-regions",
                    "version": 1,
                    "closed_surface_file": "closed_surface.stl",
                    "closed_surface_sha256": hashlib.sha256(closed.read_bytes()).hexdigest(),
                    "regions": regions,
                }),
                encoding="utf-8",
            )
            workspace = root / "workspace"
            with patch(
                "tools.oneflow_gui.mesh_project.UserDataPaths.discover",
                return_value=SimpleNamespace(mesh_generation=workspace),
            ):
                project = create_mesh_project_from_stl(
                    self.repository_root, closed, "region_project"
                )
            project_input = project["directory"] / "input"
            self.assertTrue((project_input / "region_manifest.json").is_file())
            self.assertTrue((project_input / "wall.stl").is_file())
            self.assertTrue((project_input / "inlet.stl").is_file())
            self.assertTrue((project_input / "outlet.stl").is_file())
            self.assertEqual(
                project["boundary_regions"], ["inflow", "outflow", "wall"]
            )
            self.assertIsNone(project["template_case_id"])
            self.assertEqual(
                project["mesh_case_id"], "stl_to_gmsh_generic_3d_mesh_handoff"
            )

    def test_open_mesh_project_displays_geometry_lineage(self) -> None:
        from types import SimpleNamespace

        with tempfile.TemporaryDirectory() as directory:
            workspace = Path(directory) / "mesh-generation"
            with patch(
                "tools.oneflow_gui.mesh_generation_dialog.UserDataPaths.discover",
                return_value=SimpleNamespace(mesh_generation=workspace),
            ):
                project = create_mesh_project_from_template(
                    self.repository_root,
                    "stl_to_gmsh_python_sphere_automatic_tetrahedra",
                    "lineage_display",
                )
            dialog = StlToGmshDialog(self.repository_root, project=project)
            try:
                displayed = dialog.geometry_lineage_view.toPlainText()
                self.assertIn("Provenance: UNKNOWN", displayed)
                self.assertIn("Source STL:", displayed)
            finally:
                dialog.close()

    def test_closed_shell_preset_resets_project_validation_to_generic_handoff(self) -> None:
        from types import SimpleNamespace

        with tempfile.TemporaryDirectory() as directory:
            workspace = Path(directory) / "mesh-generation"
            with patch(
                "tools.oneflow_gui.mesh_project.UserDataPaths.discover",
                return_value=SimpleNamespace(mesh_generation=workspace),
            ):
                project = create_mesh_project_from_template(
                    self.repository_root,
                    "stl_to_gmsh_python_prolate_spheroid_automatic_tetrahedra",
                    "switch_to_closed_shell",
                )
            dialog = StlToGmshDialog(self.repository_root, project=project)
            try:
                dialog.geometry_preset.setCurrentIndex(
                    dialog.geometry_preset.findData("closed_shell")
                )
                expected = "stl_to_gmsh_generic_3d_mesh_handoff"
                self.assertEqual(dialog._mesh_case_id, expected)
                saved = json.loads(
                    (project["directory"] / "project.json").read_text(encoding="utf-8")
                )
                self.assertIsNone(saved["template_case_id"])
                self.assertEqual(saved["mesh_case_id"], expected)
            finally:
                dialog.close()

    def test_project_settings_can_be_edited_and_saved_without_running_gmsh(self) -> None:
        from types import SimpleNamespace

        with tempfile.TemporaryDirectory() as directory:
            workspace = Path(directory) / "mesh-generation"
            with patch(
                "tools.oneflow_gui.mesh_project.UserDataPaths.discover",
                return_value=SimpleNamespace(mesh_generation=workspace),
            ):
                project = create_mesh_project_from_template(
                    self.repository_root,
                    "stl_to_gmsh_python_prolate_spheroid_automatic_tetrahedra",
                    "editable_spheroid",
                )
            dialog = StlToGmshDialog(self.repository_root, project=project)
            try:
                self.assertFalse(dialog.save_project_button.isHidden())
                self.assertEqual(dialog.geometry_preset.currentData(), "prolate_spheroid")
                dialog.minimum_size.setText("0.08")
                dialog.maximum_size.setText("0.16")
                dialog.boundary_name.setText("body_wall")
                dialog.volume_name.setText("spheroid_fluid")
                dialog.output_path.setText(
                    str(project["directory"] / "results" / "edited.msh")
                )
                with patch("tools.oneflow_gui.mesh_generation_dialog.QMessageBox.information"):
                    self.assertTrue(dialog._save_project_settings())
                saved = load_mesh_project(project["directory"])
                self.assertEqual(saved["mesh_size_min"], 0.08)
                self.assertEqual(saved["mesh_size_max"], 0.16)
                self.assertEqual(saved["boundary_name"], "body_wall")
                self.assertEqual(saved["volume_name"], "spheroid_fluid")
                self.assertEqual(saved["output_path"].name, "edited.msh")
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

    def test_generic_surface_cut_dialog_exposes_plane_controls_and_cli(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "generic-cut"
            dialog = GenericSurfaceCutDialog(self.repository_root)
            try:
                dialog.output_directory.setText(str(output))
                command = dialog._command()
                self.assertIn("grid.stl_import.cut_surface_cli", command)
                self.assertIn("--inlet-point", command)
                self.assertIn("--inlet-normal", command)
                self.assertIn("--outlet-point", command)
                self.assertIn("--outlet-normal", command)
                self.assertIn("--tolerance", command)
                self.assertNotIn("generate_cylinder_volume_mesh.py", command)
            finally:
                dialog.close()

    def test_volume_mesh_cli_has_beginner_friendly_defaults(self) -> None:
        self.assertTrue(DEFAULT_SPHERE_STL.is_file())
        self.assertEqual(DEFAULT_MESH_SIZE, 0.20)
        with patch.dict(os.environ, {"ONEFLOW_HOME": "C:/OneFLOWData"}):
            self.assertEqual(
                default_output_path(),
                Path("C:/OneFLOWData/mesh-generation/sphere/sphere.msh"),
            )

    def test_gmsh_can_be_discovered_from_environment_without_cli_argument(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            executable = Path(directory) / "gmsh.exe"
            executable.write_text("", encoding="ascii")
            with patch.dict(os.environ, {"GMSH_EXECUTABLE": str(executable)}):
                self.assertEqual(discover_gmsh(), str(executable.resolve()))


    def test_geometry_library_exposes_pipe_bend_parameter_editor(self) -> None:
        from tools.oneflow_gui.geometry_library_dialog import GeometryLibraryDialog

        dialog = MeshGenerationDialog(self.repository_root)
        edited_assets = []
        library = GeometryLibraryDialog(
            self.repository_root, dialog, on_edit_asset=edited_assets.append
        )
        try:
            row = next(
                index for index, asset in enumerate(library._assets)
                if asset["id"] == "geometry.pipe_bend.90_degree"
            )
            library.asset_list.setCurrentRow(row)
            self.assertTrue(library.edit_button.isEnabled())
            self.assertIn("90-Degree Circular Pipe Bend", library.details.toPlainText())
            library._edit_selected_asset()
            self.assertEqual(edited_assets[0]["id"], "geometry.pipe_bend.90_degree")
        finally:
            library.close()
            dialog.close()

    def test_stl_to_gmsh_geometry_presets_include_generated_pipe_bend(self) -> None:
        dialog = StlToGmshDialog(self.repository_root)
        try:
            index = dialog.geometry_preset.findData("pipe_bend_90_degree")
            self.assertGreaterEqual(index, 0)
            self.assertEqual(dialog.geometry_preset.itemText(index), "90-Degree Pipe Bend STL")
        finally:
            dialog.close()

    def test_stl_to_gmsh_pipe_bend_preset_generates_user_workspace_input(self) -> None:
        from types import SimpleNamespace

        with tempfile.TemporaryDirectory() as directory:
            workspace = Path(directory) / "mesh-generation"
            with patch(
                "tools.oneflow_gui.mesh_generation_dialog.UserDataPaths.discover",
                return_value=SimpleNamespace(mesh_generation=workspace),
            ):
                dialog = StlToGmshDialog(self.repository_root)
                try:
                    index = dialog.geometry_preset.findData("pipe_bend_90_degree")
                    self.assertGreaterEqual(index, 0)
                    dialog.geometry_preset.setCurrentIndex(index)
                    stl = Path(dialog.stl_path.text())
                    self.assertTrue(stl.is_file())
                    self.assertTrue(stl.with_suffix(".geometry.json").is_file())
                    self.assertTrue((stl.parent / "region_manifest.json").is_file())
                    self.assertEqual(
                        dialog._mesh_case_id, "stl_to_gmsh_pipe_bend_90_degree"
                    )
                    self.assertEqual(
                        dialog.output_path.text(),
                        str(workspace / "pipe_bend_90_degree" / "pipe_bend_90_degree.msh"),
                    )
                    self.assertFalse(dialog.boundary_name.isEnabled())
                    self.assertFalse(dialog.volume_name.isEnabled())
                    self.assertIn(str(stl), dialog.command_view.toPlainText())
                    self.assertNotIn(str(self.repository_root), str(stl))
                finally:
                    dialog.close()

    def test_development_case_selector_offers_pipe_bend_python_reference(self) -> None:
        from tools.oneflow_gui.test_dialog import TestRunnerDialog

        dialog = TestRunnerDialog(self.repository_root)
        try:
            self.assertFalse(dialog.python_development_check.isChecked())
            dialog.python_development_check.setChecked(True)
            self.assertIn("Mesh producers (", dialog.case_database_button.text())
            self.assertIn(
                "stl_to_gmsh_pipe_bend_90_degree",
                {case.case_id for case in dialog._available_database_cases_by_kind["python"]},
            )
            available = select_cases_for_selection_dialog(
                dialog._all_database_cases,
                run_native=False,
                run_python=True,
            )
            self.assertIn(
                "stl_to_gmsh_pipe_bend_90_degree",
                {case.case_id for case in available},
            )
        finally:
            dialog.close()

    def test_pipe_bend_asset_opens_matching_parameter_editor(self) -> None:
        dialog = StlGeometryDialog(
            self.repository_root, asset_id="geometry.pipe_bend.90_degree"
        )
        try:
            self.assertEqual(dialog.geometry.currentData(), "pipe_bend_90_degree")
            self.assertFalse(dialog.pipe_bend_group.isHidden())
            self.assertTrue(dialog.cylinder_group.isHidden())
            self.assertEqual(dialog.pipe_bend_radius.text(), "0.5")
            self.assertEqual(dialog.pipe_bend_bend_radius.text(), "1.0")
            self.assertEqual(dialog.pipe_bend_angle.text(), "90")
        finally:
            dialog.close()

    def test_generate_pipe_bend_stl_and_manifest_from_geometry_workbench(self) -> None:
        from types import SimpleNamespace

        with tempfile.TemporaryDirectory() as directory:
            workspace = Path(directory) / "mesh-generation"
            with patch(
                "tools.oneflow_gui.mesh_generation_dialog.UserDataPaths.discover",
                return_value=SimpleNamespace(mesh_generation=workspace),
            ):
                dialog = StlGeometryDialog(self.repository_root)
                try:
                    dialog.geometry.setCurrentIndex(
                        dialog.geometry.findData("pipe_bend_90_degree")
                    )
                    self.assertFalse(dialog.pipe_bend_group.isHidden())
                    self.assertIn("not solid walls", dialog.patch_note.text())
                    output = Path(directory) / "pipe_bend_90_degree.stl"
                    dialog.output.setText(str(output))
                    dialog._generate()
                    self.assertTrue(output.is_file())
                    manifest_path = output.with_suffix(".geometry.json")
                    self.assertTrue(manifest_path.is_file())
                    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
                    self.assertEqual(manifest["asset_id"], "geometry.pipe_bend.90_degree")
                    self.assertTrue(manifest["surface"]["watertight"])
                    self.assertLess(manifest["surface"]["relative_volume_error"], 0.01)
                    region_manifest = json.loads(
                        (output.parent / "region_manifest.json").read_text(encoding="utf-8")
                    )
                    self.assertEqual(
                        set(region_manifest["regions"]), {"wall", "inflow", "outflow"}
                    )
                    self.assertEqual(
                        region_manifest["closed_surface_sha256"],
                        manifest["stl_sha256"],
                    )
                    self.assertEqual(len(manifest["patch_triangle_indices"]["inflow"]), 48)
                    self.assertEqual(len(manifest["patch_triangle_indices"]["outflow"]), 48)
                    self.assertEqual(dialog.generated_geometry, "pipe_bend_90_degree")
                    self.assertEqual(dialog.generated_stl, output.resolve())
                    self.assertTrue(dialog.open_gmsh_button.isEnabled())
                    self.assertIn("inflow", dialog.result_view.toPlainText())
                    self.assertIn("outflow", dialog.result_view.toPlainText())
                    self.assertIn("Region manifest:", dialog.result_view.toPlainText())
                    self.assertIn("wall:", dialog.result_view.toPlainText())
                    project = create_mesh_project_from_stl(
                        self.repository_root, output, "pipe-bend-ui"
                    )
                    self.assertTrue(
                        (project["stl_path"].parent / "region_manifest.json").is_file()
                    )
                    gmsh_dialog = StlToGmshDialog(
                        self.repository_root, project=project
                    )
                    try:
                        self.assertEqual(
                            gmsh_dialog.geometry_preset.currentData(),
                            "pipe_bend_90_degree",
                        )
                        self.assertFalse(gmsh_dialog.geometry_preset.isEnabled())
                        self.assertFalse(gmsh_dialog.boundary_name.isEnabled())
                        self.assertFalse(gmsh_dialog.volume_name.isEnabled())
                        self.assertEqual(gmsh_dialog.stl_path.text(), str(project["stl_path"]))
                        self.assertEqual(
                            gmsh_dialog._mesh_case_id,
                            "stl_to_gmsh_pipe_bend_90_degree",
                        )
                        self.assertEqual(
                            gmsh_dialog.validate_button.text(),
                            "Mesh-to-CFD Validation...",
                        )
                    finally:
                        gmsh_dialog.close()
                finally:
                    dialog.close()

if __name__ == "__main__":
    unittest.main()
