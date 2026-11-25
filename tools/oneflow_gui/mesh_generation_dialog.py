"""Unified mesh-generation entry point and direct STL-to-Gmsh interface."""

from __future__ import annotations

import hashlib
import json
import math
import os
import shlex
import subprocess
import sys
from pathlib import Path

from PySide6.QtCore import QProcess
from PySide6.QtGui import QTextCursor
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPlainTextEdit,
    QInputDialog,
    QPushButton,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from .contract_workbench_dialog import ContractDomainDialog
from .mesh_layout_dialog import MeshLayoutDialog
from .mesh_project import (
    create_mesh_project_from_stl,
    create_mesh_project_from_template,
    import_mesh_project_directory,
    load_mesh_project,
    mesh_project_templates,
    replace_project_stl,
    save_mesh_project_settings,
)
from grid.stl_to_gmsh.generate_stl_geometry import generate_geometry
from grid.stl_to_gmsh.pipe_bend_geometry import generate as generate_pipe_bend
from grid.stl_to_gmsh.generate_volume_mesh import (
    DEFAULT_MESH_SIZE,
    DEFAULT_SPHERE_STL,
    default_output_path,
    discover_gmsh,
)
from .user_data_paths import UserDataPaths
from .geometry_library_dialog import GeometryLibraryDialog
from .cylinder_cut_dialog import CylinderCutWorkbenchDialog
from .generic_surface_cut_dialog import GenericSurfaceCutDialog
from .mesh_validation_dialog import MeshValidationDialog
from .oblique_cylinder_cut_dialog import ObliqueCylinderCutDialog


class MeshGenerationDialog(QDialog):
    """Navigation hub for mesh generation and conversion activities."""

    def __init__(self, repository_root: Path, parent=None) -> None:
        super().__init__(parent)
        self._repository_root = repository_root.resolve()
        self.setWindowTitle("Mesh Generation & Conversion")
        self.resize(900, 760)
        self.setMinimumSize(760, 620)

        outer = QVBoxLayout(self)
        scroll = QScrollArea(self)
        scroll.setObjectName("meshGenerationHubScrollArea")
        scroll.setWidgetResizable(True)
        content = QWidget()
        root = QVBoxLayout(content)
        root.setContentsMargins(8, 8, 8, 8)
        scroll.setWidget(content)
        outer.addWidget(scroll)
        title = QLabel("Mesh Generation & Conversion")
        title.setStyleSheet("font-size: 22px; font-weight: 650;")
        root.addWidget(title)
        intro = QLabel(
            "Generate meshes with the existing OneFLOW layout tools or external "
            "mesh engines. Each operation exposes its parameters, command and results."
        )
        intro.setWordWrap(True)
        root.addWidget(intro)

        geometry_library_button = QPushButton("Geometry Library...")
        geometry_library_button.setMinimumHeight(64)
        geometry_library_button.setToolTip(
            "Select a registered geometry asset, inspect its provenance and open "
            "the matching parameter editor or representation workflow."
        )
        geometry_library_button.clicked.connect(self._open_geometry_library)
        root.addWidget(geometry_library_button)

        cylinder_cut_button = QPushButton("Cylinder Cutting & CFD...")
        cylinder_cut_button.setMinimumHeight(64)
        cylinder_cut_button.setToolTip(
            "Cut the analytic internal-flow cylinder, validate its open surface, "
            "generate a Gmsh fluid-volume mesh and run the OneFLOW ingestion smoke test."
        )
        cylinder_cut_button.clicked.connect(self._open_cylinder_cut_workbench)
        root.addWidget(cylinder_cut_button)

        generic_cut_button = QPushButton("Generic STL Surface Cutting...")
        generic_cut_button.setMinimumHeight(64)
        generic_cut_button.setToolTip(
            "Treat any closed STL shell as an unknown surface and cut it with two "
            "user-defined planes; no analytic cylinder or pipe model is used."
        )
        generic_cut_button.clicked.connect(self._open_generic_surface_cut)
        root.addWidget(generic_cut_button)

        oblique_cut_button = QPushButton("Oblique-cut Cylinder STL to Gmsh...")
        oblique_cut_button.setMinimumHeight(64)
        oblique_cut_button.setToolTip(
            "Run the blind-STL oblique inlet cutting reference case, generate the "
            "tetrahedral Gmsh volume mesh and verify all Physical Groups."
        )
        oblique_cut_button.clicked.connect(self._open_oblique_cylinder_cut)
        root.addWidget(oblique_cut_button)

        project_heading = QLabel("Mesh Projects")
        project_heading.setStyleSheet("font-size: 16px; font-weight: 600;")
        root.addWidget(project_heading)
        project_row = QHBoxLayout()
        self.project_template = QComboBox()
        try:
            for case_id, title in mesh_project_templates(self._repository_root):
                self.project_template.addItem(title, case_id)
        except (OSError, ValueError, KeyError) as error:
            self.project_template.addItem(f"Templates unavailable: {error}", "")
        self.create_project_button = QPushButton("New from Template")
        self.create_project_button.clicked.connect(self._create_mesh_project)
        project_row.addWidget(self.project_template, 1)
        project_row.addWidget(self.create_project_button)
        root.addLayout(project_row)

        project_actions = QHBoxLayout()
        open_project_button = QPushButton("Open Existing Project...")
        open_project_button.clicked.connect(self._open_mesh_project)
        import_stl_button = QPushButton("Import STL / Grid Database...")
        import_stl_button.clicked.connect(self._import_stl_project)
        project_actions.addWidget(open_project_button)
        project_actions.addWidget(import_stl_button)
        root.addLayout(project_actions)

        layout_button = QPushButton("Mesh Layout Editor")
        layout_button.setMinimumHeight(64)
        layout_button.setToolTip("Edit structured mesh layouts and compile OneFLOW layout files.")
        layout_button.clicked.connect(self._open_mesh_layout)
        root.addWidget(layout_button)

        stl_button = QPushButton("STL to Gmsh")
        stl_button.setMinimumHeight(64)
        stl_button.setToolTip("Generate a 3-D tetrahedral volume mesh from a closed STL shell.")
        stl_button.clicked.connect(self._open_stl_to_gmsh)
        root.addWidget(stl_button)

        validation_button = QPushButton("Mesh-to-CFD Validation...")
        validation_button.setMinimumHeight(54)
        validation_button.setToolTip(
            "Run registered mesh Cases through Gmsh, CGNS conversion, OneFLOW "
            "ingestion and a short CFD smoke test. A generated project can also "
            "validate its current MSH directly from the STL-to-Gmsh page."
        )
        validation_button.clicked.connect(self._open_mesh_validation)
        root.addWidget(validation_button)

        knowledge_button = QPushButton("Mesh Knowledge")
        knowledge_button.setMinimumHeight(54)
        knowledge_button.setToolTip(
            "Browse mesh principles, formats, algorithms, implementation and evidence."
        )
        knowledge_button.clicked.connect(self._open_mesh_knowledge)
        root.addWidget(knowledge_button)

        future = QLabel(
            "More tools can be added here: STL inspection and repair, mesh conversion "
            "and mesh quality checks."
        )
        future.setWordWrap(True)
        future.setStyleSheet("color: #64748b;")
        root.addWidget(future)
        root.addStretch(1)

    def _open_geometry_library(self) -> None:
        def edit_asset(asset: dict) -> None:
            StlGeometryDialog(
                self._repository_root, self, asset_id=str(asset["id"])
            ).exec()

        GeometryLibraryDialog(
            self._repository_root, self, on_edit_asset=edit_asset
        ).exec()

    def _open_cylinder_cut_workbench(self) -> None:
        CylinderCutWorkbenchDialog(self._repository_root, self).exec()

    def _open_generic_surface_cut(self) -> None:
        GenericSurfaceCutDialog(self._repository_root, self).exec()

    def _open_oblique_cylinder_cut(self) -> None:
        ObliqueCylinderCutDialog(self._repository_root, self).exec()

    def _create_mesh_project(self) -> None:
        case_id = self.project_template.currentData()
        if not case_id:
            QMessageBox.warning(self, "No mesh template", "No supported mesh template is available.")
            return
        default_name = str(case_id).replace("stl_to_gmsh_python_", "")
        default_name = default_name.replace("_automatic_tetrahedra", "")
        name, accepted = QInputDialog.getText(
            self, "Create Mesh Project", "Project name:", text=default_name,
        )
        if not accepted:
            return
        try:
            project = create_mesh_project_from_template(
                self._repository_root, str(case_id), name
            )
        except (OSError, ValueError, KeyError) as error:
            QMessageBox.warning(self, "Cannot create mesh project", str(error))
            return
        StlToGmshDialog(self._repository_root, self, project=project).exec()

    def _import_stl_project(self) -> None:
        start = self._repository_root / "grid"
        path, _ = QFileDialog.getOpenFileName(
            self, "Import STL from Grid database or disk", str(start),
            "STL files (*.stl);;All files (*)",
        )
        if not path:
            return
        try:
            project = create_mesh_project_from_stl(
                self._repository_root, Path(path)
            )
        except (OSError, ValueError, KeyError) as error:
            QMessageBox.warning(self, "Cannot import STL", str(error))
            return
        StlToGmshDialog(self._repository_root, self, project=project).exec()

    def _open_mesh_project(self) -> None:
        project_root = UserDataPaths.discover().mesh_generation
        project_root.mkdir(parents=True, exist_ok=True)
        path = QFileDialog.getExistingDirectory(
            self, "Open OneFLOW Mesh Project", str(project_root),
        )
        if not path:
            return
        try:
            selected_directory = Path(path)
            if (selected_directory / "project.json").is_file():
                project = load_mesh_project(selected_directory)
            else:
                project = import_mesh_project_directory(
                    self._repository_root, selected_directory
                )
        except (OSError, ValueError, KeyError) as error:
            QMessageBox.warning(self, "Cannot open mesh project", str(error))
            return
        StlToGmshDialog(self._repository_root, self, project=project).exec()

    def _open_mesh_layout(self) -> None:
        MeshLayoutDialog(self._repository_root, self).exec()

    def _open_stl_geometry(self) -> None:
        StlGeometryDialog(self._repository_root, self).exec()

    def _open_stl_to_gmsh(self) -> None:
        StlToGmshDialog(self._repository_root, self).exec()

    def _open_mesh_validation(self) -> None:
        MeshValidationDialog(self._repository_root, parent=self).exec()

    def _open_mesh_knowledge(self) -> None:
        ContractDomainDialog(self._repository_root, "mesh", self).exec()


class StlGeometryDialog(QDialog):
    """Generate deterministic STL assets from the built-in geometry library."""

    def __init__(self, repository_root: Path, parent=None, *, asset_id: str | None = None) -> None:
        super().__init__(parent)
        self._repository_root = repository_root.resolve()
        self.generated_stl: Path | None = None
        self.generated_geometry: str | None = None
        self._active_geometry: str | None = None
        self.setWindowTitle("Generate STL Geometry")
        self.resize(720, 650)
        self.setMinimumSize(560, 480)
        root = QVBoxLayout(self)
        heading = QLabel("Generate STL Geometry")
        heading.setStyleSheet("font-size: 22px; font-weight: 650;")
        root.addWidget(heading)
        description = QLabel(
            "Choose a geometry from the built-in library, adjust its analytic parameters "
            "and surface discretization, then generate an STL plus a checksum-bound "
            "geometry manifest. Geometry resolution is independent of Gmsh volume-mesh size."
        )
        description.setWordWrap(True)
        root.addWidget(description)

        # Keep configuration controls in a scroll area. On Windows, the
        # combined size hints of all geometry groups and the result pane can
        # otherwise force a minimum height larger than the available screen.
        settings_scroll = QScrollArea(self)
        settings_scroll.setObjectName("geometrySettingsScrollArea")
        settings_scroll.setWidgetResizable(True)
        settings_content = QWidget()
        settings_layout = QVBoxLayout(settings_content)
        settings_layout.setContentsMargins(0, 0, 0, 0)
        settings_scroll.setWidget(settings_content)
        root.addWidget(settings_scroll, 1)

        selector_form = QFormLayout()
        self.geometry = QComboBox()
        self.geometry.addItem("Sphere", "sphere")
        self.geometry.addItem("Prolate Spheroid", "prolate_spheroid")
        self.geometry.addItem("Cylinder — Internal Flow", "cylinder_internal_flow")
        self.geometry.addItem("90° Pipe Bend — Internal Flow", "pipe_bend_90_degree")
        selector_form.addRow("Geometry asset", self.geometry)
        settings_layout.addLayout(selector_form)
        asset_to_geometry = {
            "geometry.sphere.unit": "sphere",
            "geometry.spheroid.prolate_6_to_1": "prolate_spheroid",
            "geometry.cylinder.internal_flow": "cylinder_internal_flow",
            "geometry.pipe_bend.90_degree": "pipe_bend_90_degree",
        }
        if asset_id is not None:
            geometry_id = asset_to_geometry.get(asset_id)
            if geometry_id is None:
                raise ValueError(f"No parameter editor is registered for geometry asset {asset_id!r}.")
            self.geometry.setCurrentIndex(self.geometry.findData(geometry_id))

        # Keep one widget set per geometry. Hiding inactive groups preserves
        # their independent values when the user switches between shapes.
        self.sphere_group = QGroupBox("Sphere Parameters")
        sphere_form = QFormLayout(self.sphere_group)
        self.radius = QLineEdit("1.0")
        sphere_form.addRow("Radius", self.radius)
        self.subdivisions = QLineEdit("3")
        sphere_form.addRow("Icosphere subdivision level", self.subdivisions)
        settings_layout.addWidget(self.sphere_group)

        self.spheroid_group = QGroupBox("Prolate Spheroid Parameters")
        spheroid_form = QFormLayout(self.spheroid_group)
        self.semi_axis_a = QLineEdit("3.0")
        self.semi_axis_b = QLineEdit("0.5")
        self.semi_axis_c = QLineEdit("0.5")
        axes = QHBoxLayout()
        for label, widget in (("a", self.semi_axis_a), ("b", self.semi_axis_b), ("c", self.semi_axis_c)):
            axes.addWidget(QLabel(label))
            axes.addWidget(widget)
        spheroid_form.addRow("Semi-axes", axes)
        self.spheroid_subdivisions = QLineEdit("3")
        spheroid_form.addRow("Surface subdivision level", self.spheroid_subdivisions)
        settings_layout.addWidget(self.spheroid_group)

        self.cylinder_group = QGroupBox("Cylinder Parameters")
        cylinder_form = QFormLayout(self.cylinder_group)
        self.cylinder_radius = QLineEdit("1.0")
        cylinder_form.addRow("Radius", self.cylinder_radius)
        self.length = QLineEdit("4.0")
        cylinder_form.addRow("Length", self.length)
        self.circumferential = QLineEdit("48")
        cylinder_form.addRow("Circumferential divisions", self.circumferential)
        self.axial = QLineEdit("16")
        cylinder_form.addRow("Axial divisions", self.axial)
        settings_layout.addWidget(self.cylinder_group)

        self.pipe_bend_group = QGroupBox("90-Degree Pipe-Bend Parameters")
        pipe_bend_form = QFormLayout(self.pipe_bend_group)
        self.pipe_bend_radius = QLineEdit("0.5")
        self.pipe_bend_bend_radius = QLineEdit("1.0")
        self.pipe_bend_angle = QLineEdit("90")
        self.pipe_bend_inlet_length = QLineEdit("1.5")
        self.pipe_bend_outlet_length = QLineEdit("1.5")
        self.pipe_bend_circumferential = QLineEdit("48")
        self.pipe_bend_bend_divisions = QLineEdit("24")
        self.pipe_bend_inlet_divisions = QLineEdit("12")
        self.pipe_bend_outlet_divisions = QLineEdit("12")
        pipe_bend_form.addRow("Pipe radius", self.pipe_bend_radius)
        pipe_bend_form.addRow("Centreline bend radius", self.pipe_bend_bend_radius)
        pipe_bend_form.addRow("Bend angle (degrees)", self.pipe_bend_angle)
        pipe_bend_form.addRow("Inlet straight length", self.pipe_bend_inlet_length)
        pipe_bend_form.addRow("Outlet straight length", self.pipe_bend_outlet_length)
        pipe_bend_form.addRow("Circumferential divisions", self.pipe_bend_circumferential)
        pipe_bend_form.addRow("Bend divisions", self.pipe_bend_bend_divisions)
        pipe_bend_form.addRow("Inlet / outlet divisions", self.pipe_bend_inlet_divisions)
        pipe_bend_form.addRow("Outlet divisions", self.pipe_bend_outlet_divisions)
        settings_layout.addWidget(self.pipe_bend_group)

        common_group = QGroupBox("Common Output Settings")
        common_form = QFormLayout(common_group)
        self.units = QComboBox()
        self.units.addItems(["m", "mm", "cm", "in"])
        common_form.addRow("Coordinate units", self.units)

        self.output = QLineEdit()
        output_row = QHBoxLayout()
        output_row.addWidget(self.output, 1)
        browse = QPushButton("Browse...")
        browse.clicked.connect(self._browse_output)
        output_row.addWidget(browse)
        common_form.addRow("Output STL", output_row)
        settings_layout.addWidget(common_group)

        self.patch_note = QLabel()
        self.patch_note.setWordWrap(True)
        self.patch_note.setStyleSheet("color: #475569;")
        settings_layout.addWidget(self.patch_note)

        self.generate_button = QPushButton("Generate STL")
        self.generate_button.setMinimumHeight(42)
        self.generate_button.clicked.connect(self._generate)
        root.addWidget(self.generate_button)
        self.open_gmsh_button = QPushButton("Open Generated STL in STL to Gmsh")
        self.open_gmsh_button.setEnabled(False)
        self.open_gmsh_button.clicked.connect(self._open_in_gmsh)
        root.addWidget(self.open_gmsh_button)
        self.result_view = QPlainTextEdit()
        self.result_view.setReadOnly(True)
        self.result_view.setMinimumHeight(110)
        root.addWidget(self.result_view, 1)
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.accept)
        root.addWidget(close_button)
        self.geometry.currentIndexChanged.connect(self._geometry_changed)
        self._geometry_changed()

    def _default_output(self, geometry: str) -> Path:
        return UserDataPaths.discover().mesh_generation / "generated_stl" / geometry / f"{geometry}.stl"

    def _geometry_changed(self, *_args) -> None:
        geometry = str(self.geometry.currentData())
        previous_default = self._default_output(self._active_geometry) if self._active_geometry else None
        current = Path(self.output.text()).expanduser() if self.output.text().strip() else None
        if current is None or current == previous_default:
            self.output.setText(str(self._default_output(geometry)))
        is_sphere = geometry == "sphere"
        is_spheroid = geometry == "prolate_spheroid"
        is_cylinder = geometry == "cylinder_internal_flow"
        is_pipe_bend = geometry == "pipe_bend_90_degree"
        self.sphere_group.setVisible(is_sphere)
        self.spheroid_group.setVisible(is_spheroid)
        self.cylinder_group.setVisible(is_cylinder)
        self.pipe_bend_group.setVisible(is_pipe_bend)
        self.open_gmsh_button.setEnabled(
            self.generated_stl is not None and self.generated_geometry == geometry
        )
        if is_pipe_bend:
            self.patch_note.setText(
                "The closed swept shell contains separate wall, inflow and outflow STL "
                "patches plus a checksum-bound region manifest. The generic region-aware "
                "STL-to-Gmsh route preserves these Physical Surfaces without using "
                "pipe-bend equations. Inlet/outlet disks close the fluid volume; they "
                "are not solid walls."
            )
        elif is_cylinder:
            self.patch_note.setText(
                "The closed STL contains wall, inflow and outflow patches. End caps close "
                "the shell geometrically but are not solid walls. Patch identities are "
                "stored in the companion .geometry.json manifest. The generic STL-to-Gmsh "
                "path does not yet reconstruct these separate Physical Surfaces."
            )
        else:
            self.patch_note.setText("One closed wall surface; the STL is a piecewise-planar approximation.")
        self._active_geometry = geometry

    def _browse_output(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self, "Save generated STL", self.output.text(), "STL files (*.stl);;All files (*)"
        )
        if path:
            self.output.setText(path if path.lower().endswith(".stl") else path + ".stl")

    def _generate(self) -> None:
        try:
            geometry = str(self.geometry.currentData())
            # Pass only the active geometry's parameters. Inactive panels
            # retain their values but cannot accidentally affect generation.
            if geometry == "sphere":
                kwargs = {
                    "radius": float(self.radius.text()),
                    "subdivisions": int(self.subdivisions.text()),
                }
            elif geometry == "prolate_spheroid":
                kwargs = {
                    "semi_axes": (
                        float(self.semi_axis_a.text()),
                        float(self.semi_axis_b.text()),
                        float(self.semi_axis_c.text()),
                    ),
                    "subdivisions": int(self.spheroid_subdivisions.text()),
                }
            elif geometry == "pipe_bend_90_degree":
                kwargs = {
                    "radius": float(self.pipe_bend_radius.text()),
                    "bend_radius": float(self.pipe_bend_bend_radius.text()),
                    "angle_degrees": float(self.pipe_bend_angle.text()),
                    "inlet_length": float(self.pipe_bend_inlet_length.text()),
                    "outlet_length": float(self.pipe_bend_outlet_length.text()),
                    "circumferential_divisions": int(self.pipe_bend_circumferential.text()),
                    "bend_divisions": int(self.pipe_bend_bend_divisions.text()),
                    "inlet_divisions": int(self.pipe_bend_inlet_divisions.text()),
                    "outlet_divisions": int(self.pipe_bend_outlet_divisions.text()),
                }
            else:
                kwargs = {
                    "radius": float(self.cylinder_radius.text()),
                    "length": float(self.length.text()),
                    "circumferential_divisions": int(self.circumferential.text()),
                    "axial_divisions": int(self.axial.text()),
                }
            kwargs["units"] = str(self.units.currentText())
            output = Path(self.output.text()).expanduser()
            if not self.output.text().strip():
                raise ValueError("Choose an output STL path.")
            report = (
                generate_pipe_bend(output, **kwargs)
                if geometry == "pipe_bend_90_degree"
                else generate_geometry(geometry, output, **kwargs)
            )
        except (OSError, ValueError) as error:
            QMessageBox.warning(self, "STL generation failed", str(error))
            return
        self.generated_stl = Path(report["stl"])
        self.generated_geometry = geometry
        self.open_gmsh_button.setEnabled(True)
        surface = report["surface"]
        patch_counts = {name: len(indices) for name, indices in report["patch_triangle_indices"].items()}
        result_lines = [
            "STL generated successfully.",
            f"Geometry: {geometry}",
            f"STL: {report['stl']}",
            f"Manifest: {report['manifest']}",
            f"Units: {report['units']}",
            f"Vertices: {surface['vertex_count']}",
            f"Triangles: {surface['triangle_count']}",
            f"Surface area: {surface['surface_area']:.9g}",
            f"Signed volume: {surface['signed_volume']:.9g}",
            f"Watertight: {surface['watertight']}",
            f"Patch triangles: {patch_counts}",
            f"SHA-256: {report['stl_sha256']}",
        ]
        if report.get("region_manifest_path"):
            result_lines.append(f"Region manifest: {report['region_manifest_path']}")
            for name, region_path in report.get("region_files", {}).items():
                result_lines.append(f"  {name}: {region_path}")
        self.result_view.setPlainText("\n".join(result_lines))

    def _open_in_gmsh(self) -> None:
        if self.generated_stl is None or not self.generated_stl.is_file():
            QMessageBox.warning(self, "No generated STL", "Generate an STL before opening it in Gmsh.")
            return
        try:
            project = create_mesh_project_from_stl(self._repository_root, self.generated_stl)
        except (OSError, ValueError, KeyError) as error:
            QMessageBox.warning(self, "Cannot create mesh project", str(error))
            return
        StlToGmshDialog(self._repository_root, self, project=project).exec()


class StlToGmshDialog(QDialog):
    """Configure, preview and execute the direct closed-STL volume mesher."""

    def __init__(self, repository_root: Path, parent=None, *, project: dict | None = None) -> None:
        super().__init__(parent)
        self._repository_root = repository_root.resolve()
        self._project = project
        self._project_directory = Path(project["directory"]) if project else None
        self._mesh_case_id = (
            str(project.get("mesh_case_id"))
            if project else "stl_to_gmsh_sphere_automatic_tetrahedra"
        )
        self._process = QProcess(self)
        self._process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self._process.readyReadStandardOutput.connect(self._read_output)
        self._process.finished.connect(self._process_finished)
        self._process.errorOccurred.connect(self._process_error)

        self.setWindowTitle("STL to Gmsh - 3D Volume Mesh")
        self.resize(900, 760)
        root = QVBoxLayout(self)

        heading = QLabel("STL to Gmsh")
        heading.setStyleSheet("font-size: 22px; font-weight: 650;")
        root.addWidget(heading)
        description = QLabel(
            "Create a tetrahedral volume mesh from a closed, connected and consistently "
            "oriented STL shell. Gmsh is the meshing engine; this page configures and "
            "runs the OneFLOW command-line adapter."
        )
        description.setWordWrap(True)
        root.addWidget(description)
        if project:
            project_label = QLabel(
                f"Project: {project['name']}  |  Directory: {project['directory']}"
            )
            project_label.setWordWrap(True)
            project_label.setStyleSheet("color: #475569;")
            root.addWidget(project_label)
            root.addWidget(QLabel("Geometry provenance and derived artifacts"))
            self.geometry_lineage_view = QPlainTextEdit()
            self.geometry_lineage_view.setReadOnly(True)
            self.geometry_lineage_view.setMaximumHeight(130)
            root.addWidget(self.geometry_lineage_view)
            self._refresh_geometry_lineage()

        form = QFormLayout()
        self.geometry_preset = QComboBox()
        self.geometry_preset.addItem("Unit Sphere STL", "sphere")
        self.geometry_preset.addItem("6:1 Prolate Spheroid STL", "prolate_spheroid")
        self.geometry_preset.addItem("Internal-flow cylinder STL", "cylinder_internal_flow")
        self.geometry_preset.addItem("90-Degree Pipe Bend STL", "pipe_bend_90_degree")
        self.geometry_preset.addItem("Closed STL shell (genus zero)", "closed_shell")
        self.geometry_preset.setToolTip(
            "Selects a bundled or generated geometry, its STL input, validation Case and default output location."
        )
        form.addRow("Geometry preset", self.geometry_preset)

        stl_row = QHBoxLayout()
        self.stl_path = QLineEdit()
        default_sphere_stl = DEFAULT_SPHERE_STL
        if project:
            self.stl_path.setText(str(project["stl_path"]))
        elif default_sphere_stl.is_file():
            self.stl_path.setText(str(default_sphere_stl))
        self.stl_path.setPlaceholderText("Select a closed STL surface...")
        stl_browse = QPushButton("Browse...")
        stl_browse.clicked.connect(self._browse_stl)
        stl_row.addWidget(self.stl_path, 1)
        stl_row.addWidget(stl_browse)
        form.addRow("Input STL", stl_row)

        self.density = QComboBox()
        self.density.addItem("Coarse (0.40)", "coarse")
        self.density.addItem("Medium (0.20)", "medium")
        self.density.addItem("Fine (0.10)", "fine")
        self.density.addItem("Custom size range", "custom")
        self.density.setCurrentIndex(self.density.findData("medium"))
        self.density.currentIndexChanged.connect(self._density_changed)
        form.addRow("Mesh density", self.density)

        self.mesh_size = QLineEdit(f"{DEFAULT_MESH_SIZE:.2f}")
        self.mesh_size.setToolTip("One target size; maps to both MeshSizeMin and MeshSizeMax.")
        form.addRow("Target size", self.mesh_size)

        self.minimum_size = QLineEdit("0.10")
        self.maximum_size = QLineEdit("0.20")
        size_row = QHBoxLayout()
        size_row.addWidget(QLabel("Minimum"))
        size_row.addWidget(self.minimum_size)
        size_row.addWidget(QLabel("Maximum"))
        size_row.addWidget(self.maximum_size)

        self.output_path = QLineEdit(
            str(project["output_path"]) if project else str(default_output_path())
        )
        output_row = QHBoxLayout()
        output_row.addWidget(self.output_path, 1)
        output_browse = QPushButton("Browse...")
        output_browse.clicked.connect(self._browse_output)
        output_row.addWidget(output_browse)
        form.addRow("Output mesh", output_row)
        root.addLayout(form)

        advanced_toggle = QPushButton("Advanced settings")
        advanced_toggle.setCheckable(True)
        root.addWidget(advanced_toggle)
        self.advanced_settings = QWidget()
        advanced_form = QFormLayout(self.advanced_settings)
        advanced_form.addRow("Size bounds", size_row)
        self.gmsh_path = QLineEdit(discover_gmsh() or "")
        self.gmsh_path.setPlaceholderText("Auto-detected; set only if needed")
        gmsh_row = QHBoxLayout()
        gmsh_row.addWidget(self.gmsh_path, 1)
        gmsh_browse = QPushButton("Browse...")
        gmsh_browse.clicked.connect(self._browse_gmsh)
        gmsh_row.addWidget(gmsh_browse)
        advanced_form.addRow("Gmsh executable", gmsh_row)
        self.boundary_name = QLineEdit("wall")
        self.volume_name = QLineEdit("fluid_domain")
        groups_row = QHBoxLayout()
        groups_row.addWidget(QLabel("Surface"))
        groups_row.addWidget(self.boundary_name)
        groups_row.addWidget(QLabel("Volume"))
        groups_row.addWidget(self.volume_name)
        advanced_form.addRow("Physical Groups", groups_row)
        self.advanced_settings.setVisible(False)
        root.addWidget(self.advanced_settings)
        advanced_toggle.toggled.connect(self.advanced_settings.setVisible)
        command_heading = QLabel("Command Preview")
        command_heading.setStyleSheet("font-weight: 650;")
        root.addWidget(command_heading)
        self.command_view = QPlainTextEdit()
        self.command_view.setReadOnly(True)
        self.command_view.setMaximumBlockCount(100)
        self.command_view.setMinimumHeight(110)
        root.addWidget(self.command_view)

        actions = QHBoxLayout()
        self.refresh_button = QPushButton("Refresh Command")
        self.refresh_button.clicked.connect(self._refresh_command)
        actions.addWidget(self.refresh_button)
        self.copy_button = QPushButton("Copy Command")
        self.copy_button.clicked.connect(self._copy_command)
        actions.addWidget(self.copy_button)
        actions.addStretch(1)
        self.save_project_button = QPushButton("Save Project Settings")
        self.save_project_button.setVisible(project is not None)
        self.save_project_button.clicked.connect(self._save_project_settings)
        actions.addWidget(self.save_project_button)
        self.run_button = QPushButton("Generate Mesh")
        self.run_button.setMinimumWidth(150)
        self.run_button.clicked.connect(self._run)
        actions.addWidget(self.run_button)
        self.validate_button = QPushButton("Mesh-to-CFD Validation...")
        self.validate_button.setEnabled(Path(self.output_path.text()).is_file())
        self.validate_button.setToolTip(
            "Validate the current MSH file without rerunning Gmsh."
        )
        self.validate_button.clicked.connect(self._validate_current_mesh)
        actions.addWidget(self.validate_button)
        root.addLayout(actions)

        root.addWidget(QLabel("Execution Output"))
        self.output_view = QPlainTextEdit()
        self.output_view.setReadOnly(True)
        self.output_view.setMinimumHeight(140)
        root.addWidget(self.output_view, 1)

        self.status_label = QLabel("Ready. Generated files are written outside the repository.")
        self.status_label.setWordWrap(True)
        self.status_label.setStyleSheet("color: #475569;")
        root.addWidget(self.status_label)

        for widget in (
            self.stl_path, self.mesh_size, self.minimum_size, self.maximum_size,
            self.output_path, self.gmsh_path, self.boundary_name, self.volume_name,
        ):
            if isinstance(widget, QLineEdit):
                widget.textChanged.connect(self._refresh_command)
        if project:
            self.boundary_name.setText(str(project.get("boundary_name", "wall")))
            self.volume_name.setText(str(project.get("volume_name", "fluid_domain")))
            template_id = project.get("template_case_id")
            if project.get("geometry") == "cylinder_internal_flow":
                self._mesh_case_id = "stl_to_gmsh_existing_mesh_3d_cylinder_internal_flow"
                self.geometry_preset.setCurrentIndex(
                    self.geometry_preset.findData("cylinder_internal_flow")
                )
                self.boundary_name.setEnabled(False)
                self.volume_name.setEnabled(False)
            elif project.get("geometry") == "pipe_bend_90_degree":
                self._mesh_case_id = "stl_to_gmsh_pipe_bend_90_degree"
                self.geometry_preset.setCurrentIndex(
                    self.geometry_preset.findData("pipe_bend_90_degree")
                )
                self.geometry_preset.setEnabled(False)
                self.boundary_name.setEnabled(False)
                self.volume_name.setEnabled(False)
            elif template_id == "stl_to_gmsh_python_prolate_spheroid_automatic_tetrahedra":
                self.geometry_preset.setCurrentIndex(self.geometry_preset.findData("prolate_spheroid"))
            elif template_id == "stl_to_gmsh_python_sphere_automatic_tetrahedra":
                self.geometry_preset.setCurrentIndex(self.geometry_preset.findData("sphere"))
            self.density.setCurrentIndex(self.density.findData("custom"))
            target = float(project.get("mesh_size", DEFAULT_MESH_SIZE))
            minimum = float(project.get("mesh_size_min", target * 0.5))
            maximum = float(project.get("mesh_size_max", target))
            self.minimum_size.setText(f"{minimum:.6g}")
            self.maximum_size.setText(f"{maximum:.6g}")
        self._active_geometry_preset = str(self.geometry_preset.currentData() or "sphere")
        self.geometry_preset.currentIndexChanged.connect(self._geometry_preset_changed)
        self._density_changed()
        self._refresh_command()

    def _geometry_preset_paths(self, preset: str) -> tuple[Path, str, str] | None:
        """Return the bundled STL, native validation Case and Python template Case."""
        if preset == "cylinder_internal_flow":
            case_directory = (
                self._repository_root / "grid" / "stl_to_gmsh" / "3d"
                / "cylinder_internal_flow" / "automatic_tetrahedra"
            )
            return (
                case_directory / "cpp" / "input" / "model.stl",
                "stl_to_gmsh_cylinder_internal_flow_automatic_tetrahedra",
                "stl_to_gmsh_python_cylinder_internal_flow_automatic_tetrahedra",
            )
        if preset not in {"sphere", "prolate_spheroid"}:
            return None
        case_directory = (
            self._repository_root / "grid" / "stl_to_gmsh" / "3d" / preset
            / "automatic_tetrahedra"
        )
        stl = case_directory / "python" / "input" / "model.stl"
        return (
            stl,
            f"stl_to_gmsh_{preset}_automatic_tetrahedra",
            f"stl_to_gmsh_python_{preset}_automatic_tetrahedra",
        )

    def _geometry_default_output(self, preset: str) -> Path:
        return UserDataPaths.discover().mesh_generation / preset / f"{preset}.msh"

    def _geometry_preset_changed(self, *_args) -> None:
        """Switch every geometry-dependent input when a preset is selected."""
        preset = str(self.geometry_preset.currentData() or "closed_shell")
        previous = self._active_geometry_preset
        previous_output = self._geometry_default_output(previous)
        current_output = Path(self.output_path.text()).expanduser()
        # The dialog's initial sphere output comes from generate_volume_mesh's
        # default_output_path(), while subsequent preset defaults are resolved
        # through UserDataPaths. Treat both as untouched defaults so switching
        # presets also works when the user-data directory is overridden (e.g.
        # in tests or by a configured workspace).
        initial_default_output = Path(default_output_path()).expanduser()
        output_is_default = (
            current_output.resolve() == previous_output.resolve()
            or current_output.resolve() == initial_default_output.resolve()
        )
        self.boundary_name.setEnabled(True)
        self.volume_name.setEnabled(True)
        definition = self._geometry_preset_paths(preset)

        if preset == "pipe_bend_90_degree":
            stl = (
                UserDataPaths.discover().mesh_generation
                / "generated_stl" / preset / f"{preset}.stl"
            )
            try:
                stl.parent.mkdir(parents=True, exist_ok=True)
                generate_pipe_bend(stl, units="m")
                if self._project is not None:
                    self._project = replace_project_stl(self._project, stl)
                    manifest_path = Path(self._project["directory"]) / "project.json"
                    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
                    manifest["template_case_id"] = None
                    manifest["mesh_case_id"] = "stl_to_gmsh_pipe_bend_90_degree"
                    manifest_path.write_text(
                        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8",
                        newline="\n",
                    )
                    self._project = load_mesh_project(Path(self._project["directory"]))
                    self._project_directory = Path(self._project["directory"])
                    self.stl_path.setText(str(self._project["stl_path"]))
                else:
                    self.stl_path.setText(str(stl))
                self._mesh_case_id = "stl_to_gmsh_pipe_bend_90_degree"
                self.boundary_name.setEnabled(False)
                self.volume_name.setEnabled(False)
                if output_is_default:
                    self.output_path.setText(str(self._geometry_default_output(preset)))
            except (OSError, ValueError, KeyError) as error:
                QMessageBox.warning(self, "Cannot prepare pipe-bend STL", str(error))
                self.geometry_preset.blockSignals(True)
                self.geometry_preset.setCurrentIndex(
                    self.geometry_preset.findData(previous)
                )
                self.geometry_preset.blockSignals(False)
                return
            self._active_geometry_preset = preset
            self._refresh_command()
            return

        if definition is None:
            generic_case_id = "stl_to_gmsh_generic_3d_mesh_handoff"
            if self._project is not None:
                try:
                    manifest_path = Path(self._project["directory"]) / "project.json"
                    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
                    manifest["template_case_id"] = None
                    manifest["mesh_case_id"] = generic_case_id
                    manifest_path.write_text(
                        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8",
                        newline="\n",
                    )
                    self._project = load_mesh_project(Path(self._project["directory"]))
                except (OSError, ValueError, KeyError) as error:
                    QMessageBox.warning(
                        self, "Cannot select generic STL geometry", str(error)
                    )
                    self.geometry_preset.blockSignals(True)
                    self.geometry_preset.setCurrentIndex(
                        self.geometry_preset.findData(previous)
                    )
                    self.geometry_preset.blockSignals(False)
                    return
            self._mesh_case_id = generic_case_id
            self._active_geometry_preset = preset
            self._refresh_command()
            return

        stl, native_case_id, project_template_id = definition
        if not stl.is_file():
            QMessageBox.warning(
                self, "Geometry preset unavailable",
                f"The bundled STL for '{preset}' does not exist: {stl}",
            )
            self.geometry_preset.blockSignals(True)
            self.geometry_preset.setCurrentIndex(
                self.geometry_preset.findData(previous)
            )
            self.geometry_preset.blockSignals(False)
            return

        if preset == "cylinder_internal_flow":
            manifest_path = stl.with_suffix(".geometry.json")
            try:
                metadata = json.loads(manifest_path.read_text(encoding="utf-8"))
                valid_asset = (
                    metadata.get("schema") == "oneflow.stl-geometry-manifest"
                    and metadata.get("version") == 1
                    and metadata.get("asset_id") == "geometry.cylinder.internal_flow"
                    and metadata.get("geometry") == "cylinder_internal_flow"
                    and metadata.get("stl_sha256")
                    == hashlib.sha256(stl.read_bytes()).hexdigest()
                    and metadata.get("surface", {}).get("watertight") is True
                    and set(metadata.get("patch_triangle_indices", {}))
                    == {"wall", "inflow", "outflow"}
                )
                if not valid_asset:
                    raise ValueError("Cylinder STL provenance or patch contract is invalid.")
            except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as error:
                QMessageBox.warning(
                    self, "Invalid cylinder geometry asset",
                    f"The standard cylinder asset or its manifest is invalid: {error}",
                )
                # Keep the selector and the dependent input path in sync. A
                # failed validation must not leave the new preset displayed
                # while the previous geometry remains selected as input.
                self.geometry_preset.blockSignals(True)
                self.geometry_preset.setCurrentIndex(
                    self.geometry_preset.findData(previous)
                )
                self.geometry_preset.blockSignals(False)
                return

        try:
            if self._project is not None:
                self._project = replace_project_stl(self._project, stl)
                manifest_path = Path(self._project["directory"]) / "project.json"
                manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
                manifest["template_case_id"] = project_template_id
                manifest["mesh_case_id"] = project_template_id
                manifest_path.write_text(
                    json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8",
                    newline="\n",
                )
                self._project["template_case_id"] = project_template_id
                self._project["mesh_case_id"] = project_template_id
                self._project = load_mesh_project(Path(self._project["directory"]))
                self.stl_path.setText(str(self._project["stl_path"]))
                self._mesh_case_id = project_template_id
            else:
                self.stl_path.setText(str(stl))
                self._mesh_case_id = native_case_id
                if output_is_default:
                    self.output_path.setText(str(self._geometry_default_output(preset)))
        except (OSError, ValueError, KeyError) as error:
            QMessageBox.warning(self, "Cannot switch geometry preset", str(error))
            return

        self._active_geometry_preset = preset
        self._refresh_command()

    def _browse_stl(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Select STL surface", str(Path.home()), "STL files (*.stl);;All files (*)"
        )
        if path:
            selected = Path(path).expanduser().resolve()
            previous = Path(self.stl_path.text()).expanduser().resolve() if self.stl_path.text().strip() else None
            same_project_input = self._project is not None and selected == previous
            if self._project is not None and not same_project_input:
                try:
                    self._project = replace_project_stl(self._project, selected)
                except (OSError, ValueError, KeyError) as error:
                    QMessageBox.warning(self, "Cannot replace project STL", str(error))
                    return
                self._project_directory = Path(self._project["directory"])
                self._mesh_case_id = str(self._project["mesh_case_id"])
            self.stl_path.setText(str(self._project["stl_path"]) if self._project else str(selected))
            template_id = self._project.get("template_case_id") if self._project else None
            if selected == DEFAULT_SPHERE_STL.resolve():
                self.geometry_preset.setCurrentIndex(self.geometry_preset.findData("sphere"))
                if self._project is None:
                    self._mesh_case_id = "stl_to_gmsh_sphere_automatic_tetrahedra"
            elif same_project_input and template_id == "stl_to_gmsh_python_prolate_spheroid_automatic_tetrahedra":
                self.geometry_preset.setCurrentIndex(self.geometry_preset.findData("prolate_spheroid"))
            elif same_project_input and template_id == "stl_to_gmsh_python_sphere_automatic_tetrahedra":
                self.geometry_preset.setCurrentIndex(self.geometry_preset.findData("sphere"))
            else:
                self.geometry_preset.setCurrentIndex(self.geometry_preset.findData("closed_shell"))
                if self._project is None:
                    self._mesh_case_id = "stl_to_gmsh_generic_3d_mesh_handoff"

    def _browse_output(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Gmsh mesh", self.output_path.text(), "Gmsh mesh (*.msh)"
        )
        if path:
            if not path.lower().endswith(".msh"):
                path += ".msh"
            self.output_path.setText(path)

    def _browse_gmsh(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Select Gmsh executable", str(Path.home()),
            "Gmsh executable (gmsh.exe gmsh);;All files (*)"
        )
        if path:
            self.gmsh_path.setText(path)

    def _density_changed(self, *_args) -> None:
        preset = self.density.currentData()
        sizes = {"coarse": "0.40", "medium": "0.20", "fine": "0.10"}
        custom = preset == "custom"
        self.mesh_size.setEnabled(not custom)
        self.minimum_size.setEnabled(custom)
        self.maximum_size.setEnabled(custom)
        if not custom:
            self.mesh_size.setText(sizes[preset])
        self._refresh_command()

    def _arguments(self) -> list[str]:
        stl = self.stl_path.text().strip()
        output = self.output_path.text().strip()
        if not stl:
            raise ValueError("Select an input STL file.")
        if not output:
            raise ValueError("Choose an output .msh path.")
        if Path(output).suffix.lower() != ".msh":
            raise ValueError("Output mesh path must use the .msh extension.")
        args = [
            str(self._repository_root / "grid" / "stl_to_gmsh" / "generate_volume_mesh.py"),
            "--stl", str(Path(stl).expanduser()),
            "--output", str(Path(output).expanduser()),
        ]
        if self.density.currentData() == "custom":
            minimum = float(self.minimum_size.text())
            maximum = float(self.maximum_size.text())
            if (
                not math.isfinite(minimum) or not math.isfinite(maximum)
                or minimum <= 0 or maximum <= 0 or minimum > maximum
            ):
                raise ValueError("Mesh size bounds must be positive and minimum must not exceed maximum.")
            args.extend(["--mesh-size-min", self.minimum_size.text().strip(),
                         "--mesh-size-max", self.maximum_size.text().strip()])
        else:
            size = float(self.mesh_size.text())
            if not math.isfinite(size) or size <= 0:
                raise ValueError("Mesh size must be positive.")
            args.extend(["--mesh-size", self.mesh_size.text().strip()])
        gmsh = self.gmsh_path.text().strip()
        if gmsh:
            args.extend(["--gmsh", gmsh])
        boundary = self.boundary_name.text().strip()
        volume = self.volume_name.text().strip()
        if not boundary or not volume:
            raise ValueError("Physical Group names must not be empty.")
        args.extend(["--boundary-name", boundary, "--volume-name", volume])
        return args

    def _command_text(self) -> str:
        args = [sys.executable, *self._arguments()]
        if os.name == "nt":
            return subprocess.list2cmdline(args)
        return shlex.join(args)

    def _refresh_command(self, *_args) -> None:
        try:
            self.command_view.setPlainText(self._command_text())
            idle = self._process.state() == QProcess.ProcessState.NotRunning
            self.run_button.setEnabled(idle)
            self.validate_button.setEnabled(
                idle and Path(self.output_path.text()).expanduser().is_file()
            )
            self.status_label.setText("Command is valid.")
        except (ValueError, OSError) as error:
            self.command_view.setPlainText(f"Command cannot be prepared:\n{error}")
            self.run_button.setEnabled(False)
            self.status_label.setText("Check the input parameters.")

    def _copy_command(self) -> None:
        if self.run_button.isEnabled():
            from PySide6.QtWidgets import QApplication
            QApplication.clipboard().setText(self.command_view.toPlainText())
            self.status_label.setText("Command copied to clipboard.")

    def _save_project_settings(
        self, _checked: bool = False, *, show_message: bool = True,
    ) -> bool:
        if self._project is None:
            return True
        try:
            self._arguments()
            output = Path(self.output_path.text()).expanduser().resolve()
            settings = {
                "output_path": output,
                "boundary_name": self.boundary_name.text().strip(),
                "volume_name": self.volume_name.text().strip(),
            }
            if self.density.currentData() == "custom":
                settings["mesh_size_min"] = float(self.minimum_size.text())
                settings["mesh_size_max"] = float(self.maximum_size.text())
            else:
                settings["mesh_size"] = float(self.mesh_size.text())
            self._project = save_mesh_project_settings(self._project, settings)
        except (OSError, ValueError, KeyError) as error:
            QMessageBox.warning(self, "Cannot save mesh project", str(error))
            return False
        self.status_label.setText("Project settings saved.")
        if show_message:
            QMessageBox.information(
                self, "Mesh project saved",
                f"Settings were saved to {self._project['directory'] / 'project.json'}",
            )
        return True

    def _run(self) -> None:
        try:
            arguments = self._arguments()
            stl = Path(self.stl_path.text()).expanduser().resolve()
            if not stl.is_file():
                raise ValueError(f"Input STL does not exist: {stl}")
            output = Path(self.output_path.text()).expanduser().resolve()
            if output == stl:
                raise ValueError("Output mesh must not overwrite the input STL.")
            output.parent.mkdir(parents=True, exist_ok=True)
            if not self._save_project_settings(show_message=False):
                return
        except (OSError, ValueError, KeyError) as error:
            QMessageBox.warning(self, "Cannot generate mesh", str(error))
            return

        script, *script_arguments = arguments
        self.output_view.clear()
        self.output_view.appendPlainText(f"Working directory: {self._repository_root}")
        self.output_view.appendPlainText(f"Program: {sys.executable}")
        self.output_view.appendPlainText(f"Script: {script}")
        rendered_arguments = (
            subprocess.list2cmdline(script_arguments)
            if os.name == "nt" else shlex.join(script_arguments)
        )
        self.output_view.appendPlainText(f"Arguments: {rendered_arguments}")
        self._process.setWorkingDirectory(str(self._repository_root))
        self.run_button.setEnabled(False)
        self.validate_button.setEnabled(False)
        self.status_label.setText("Starting mesh generation...")
        self._process.start(sys.executable, [script, *script_arguments])

    def _read_output(self) -> None:
        data = bytes(self._process.readAllStandardOutput()).decode("utf-8", errors="replace")
        if data:
            cursor = self.output_view.textCursor()
            cursor.movePosition(QTextCursor.MoveOperation.End)
            cursor.insertText(data)
            self.output_view.setTextCursor(cursor)
            self.output_view.ensureCursorVisible()

    def _refresh_geometry_lineage(self) -> None:
        if self._project is None:
            return
        lineage = self._project.get("geometry_lineage") or {}
        provenance = lineage.get("provenance") or {}
        lines = [
            f"Provenance: {provenance.get('kind', 'unknown').upper()}",
            str(provenance.get("notes", "No verified source record is available.")),
        ]
        boundary_regions = (
            self._project.get("boundary_regions")
            or lineage.get("boundary_patches")
            or []
        )
        if boundary_regions:
            lines.append(f"Boundary patches: {', '.join(boundary_regions)}")
        asset_id = lineage.get("asset_id")
        if asset_id:
            lines.append(f"Geometry asset: {asset_id}")
        if lineage.get("manifest_status"):
            lines.append(f"Geometry manifest: UNTRUSTED ({lineage['manifest_status']})")
        source = lineage.get("source_representation") or {}
        if source:
            lines.append(
                f"Source STL: {source.get('path', 'unknown')} | "
                f"SHA-256: {source.get('sha256', 'unknown')}"
            )
        definition_hash = lineage.get("definition_sha256")
        if definition_hash:
            lines.append(f"Geometry definition SHA-256: {definition_hash}")
        generator = lineage.get("generator") or {}
        if generator:
            lines.append(
                f"Generator: {generator.get('name', 'unknown')} "
                f"{generator.get('version', '')}"
            )
            if generator.get("source_sha256"):
                lines.append(f"Generator source SHA-256: {generator['source_sha256']}")
        mesh = lineage.get("derived_mesh") or {}
        if mesh:
            lines.append(
                f"Gmsh MSH: {mesh.get('path', 'unknown')} | "
                f"SHA-256: {mesh.get('sha256', 'unknown')}"
            )
        for artifact in lineage.get("downstream_artifacts", []):
            lines.append(
                f"{artifact.get('role', 'Derived artifact')}: "
                f"{artifact.get('path', 'unknown')} | "
                f"SHA-256: {artifact.get('sha256', 'unknown')}"
            )
        self.geometry_lineage_view.setPlainText("\n".join(lines))

    def _validate_current_mesh(self) -> None:
        mesh_path = Path(self.output_path.text()).expanduser().resolve()
        if not mesh_path.is_file():
            QMessageBox.warning(
                self, "Mesh file not found",
                f"Generate or select an existing MSH file first: {mesh_path}",
            )
            return
        from .mesh_validation_dialog import MeshValidationDialog

        dialog = MeshValidationDialog(
            self._repository_root,
            mesh_file=mesh_path,
            project_directory=self._project_directory or mesh_path.parent,
            mesh_case_id=self._mesh_case_id,
            parent=self,
        )
        dialog.exec()
        if self._project_directory is not None:
            try:
                self._project = load_mesh_project(self._project_directory)
                self._refresh_geometry_lineage()
            except (OSError, ValueError, KeyError) as error:
                self.status_label.setText(f"Cannot refresh project geometry lineage: {error}")

    def _process_finished(self, exit_code: int, _exit_status) -> None:
        self._read_output()
        self.run_button.setEnabled(True)
        if exit_code == 0:
            mesh_path = Path(self.output_path.text()).expanduser()
            self.validate_button.setEnabled(mesh_path.is_file())
            refresh_error = None
            if self._project_directory is not None:
                try:
                    self._project = load_mesh_project(self._project_directory)
                    self._refresh_geometry_lineage()
                except (OSError, ValueError, KeyError) as error:
                    refresh_error = str(error)
            message = (
                "Completed and validated. Outputs: "
                f"{mesh_path}, {mesh_path.with_suffix('.geo')}, "
                f"{mesh_path.with_suffix('.report.json')}"
            )
            if refresh_error:
                message += f" Project lineage refresh failed: {refresh_error}"
            self.status_label.setText(message)
        else:
            self.status_label.setText(f"Generation failed (exit code {exit_code}). See execution output.")

    def _process_error(self, _error) -> None:
        self._read_output()
        self.run_button.setEnabled(True)
        self.validate_button.setEnabled(False)
        self.status_label.setText("Mesh generation process error. See execution output.")

    def closeEvent(self, event) -> None:
        if self._process.state() != QProcess.ProcessState.NotRunning:
            QMessageBox.information(
                self,
                "Mesh generation is running",
                "Wait for the mesh-generation process to finish before closing this page.",
            )
            event.ignore()
            return
        event.accept()


__all__ = ["MeshGenerationDialog", "StlToGmshDialog"]
