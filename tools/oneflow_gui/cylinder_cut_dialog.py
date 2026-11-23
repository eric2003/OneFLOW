"""GUI for analytic cylinder cutting, Gmsh volume meshing and OneFLOW validation."""
from __future__ import annotations
import math
import os
import shlex
import subprocess
import sys
from pathlib import Path
from PySide6.QtCore import QProcess
from PySide6.QtGui import QTextCursor
from PySide6.QtWidgets import (
    QDialog, QFileDialog, QFormLayout, QGroupBox, QHBoxLayout, QLabel,
    QLineEdit, QMessageBox, QPlainTextEdit, QPushButton, QVBoxLayout,
)
from grid.stl_to_gmsh.generate_volume_mesh import discover_gmsh
from .user_data_paths import UserDataPaths

CYLINDER_MESH_CASE_ID = "stl_to_gmsh_existing_mesh_3d_cylinder_internal_flow"


class CylinderCutWorkbenchDialog(QDialog):
    """Generate a cut-cylinder surface and volume mesh, then hand it to CFD validation."""

    def __init__(self, repository_root: Path, parent=None) -> None:
        super().__init__(parent)
        self._repository_root = Path(repository_root).resolve()
        self._process = QProcess(self)
        self._process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self._process.readyReadStandardOutput.connect(self._read_output)
        self._process.finished.connect(self._process_finished)
        self._process.errorOccurred.connect(self._process_error)
        self.setWindowTitle("Cylinder Cutting & CFD")
        self.resize(820, 760)
        self.setMinimumSize(640, 560)
        root = QVBoxLayout(self)
        heading = QLabel("Analytic Cylinder Cutting & CFD")
        heading.setStyleSheet("font-size: 22px; font-weight: 650;")
        root.addWidget(heading)
        description = QLabel(
            "Start from the Geometry Library cylinder definition. Set two axial cutting "
            "planes, generate the open lateral STL, independently verify its two boundary "
            "loops, and create a tetrahedral fluid-volume mesh with Gmsh."
        )
        description.setWordWrap(True)
        root.addWidget(description)

        asset_group = QGroupBox("Geometry Asset")
        asset_form = QFormLayout(asset_group)
        asset_form.addRow("Asset", QLabel("Cylinder — Internal Flow"))
        asset_form.addRow("Asset ID", QLabel("geometry.cylinder.internal_flow"))
        asset_form.addRow("Axis / center", QLabel("X / (0, 0, 0)"))
        root.addWidget(asset_group)

        parameter_group = QGroupBox("Cutting & Surface Parameters")
        form = QFormLayout(parameter_group)
        self.radius = QLineEdit("1.0")
        self.axis_start = QLineEdit("-2.0")
        self.axis_end = QLineEdit("2.0")
        self.cut_start = QLineEdit("-1.5")
        self.cut_end = QLineEdit("1.5")
        self.circumferential = QLineEdit("48")
        self.axial = QLineEdit("12")
        form.addRow("Radius", self.radius)
        form.addRow("Original axial start", self.axis_start)
        form.addRow("Original axial end", self.axis_end)
        form.addRow("Cut plane 1 (x)", self.cut_start)
        form.addRow("Cut plane 2 (x)", self.cut_end)
        form.addRow("Circumferential divisions", self.circumferential)
        form.addRow("Axial divisions", self.axial)
        root.addWidget(parameter_group)

        mesh_group = QGroupBox("Gmsh Volume Mesh")
        mesh_form = QFormLayout(mesh_group)
        self.mesh_size = QLineEdit("0.2")
        mesh_form.addRow("Target mesh size", self.mesh_size)
        self.gmsh_path = QLineEdit(discover_gmsh() or "")
        self.gmsh_path.setPlaceholderText("Auto-detected from PATH; browse if needed")
        gmsh_row = QHBoxLayout()
        gmsh_row.addWidget(self.gmsh_path, 1)
        gmsh_browse = QPushButton("Browse...")
        gmsh_browse.clicked.connect(self._browse_gmsh)
        gmsh_row.addWidget(gmsh_browse)
        mesh_form.addRow("Gmsh executable", gmsh_row)
        self.output_directory = QLineEdit(str(
            UserDataPaths.discover().mesh_generation / "cylinder_internal_flow" / "cut_workflow"
        ))
        output_row = QHBoxLayout()
        output_row.addWidget(self.output_directory, 1)
        output_browse = QPushButton("Browse...")
        output_browse.clicked.connect(self._browse_output)
        output_row.addWidget(output_browse)
        mesh_form.addRow("Output directory", output_row)
        root.addWidget(mesh_group)

        command_heading = QLabel("Command Preview")
        command_heading.setStyleSheet("font-weight: 650;")
        root.addWidget(command_heading)
        self.command_view = QPlainTextEdit()
        self.command_view.setReadOnly(True)
        self.command_view.setMaximumHeight(100)
        root.addWidget(self.command_view)
        actions = QHBoxLayout()
        self.generate_button = QPushButton("Generate Cut Surface + Volume Mesh")
        self.generate_button.clicked.connect(self._generate)
        actions.addWidget(self.generate_button)
        self.run_oneflow_button = QPushButton("Run OneFLOW...")
        self.run_oneflow_button.setEnabled(self._mesh_path().is_file())
        self.run_oneflow_button.setToolTip(
            "Open Mesh-to-CFD Validation for the generated MSH. The registered CFD consumer "
            "is a short smoke-test template, not a cylinder accuracy case."
        )
        self.run_oneflow_button.clicked.connect(self._run_oneflow)
        actions.addWidget(self.run_oneflow_button)
        root.addLayout(actions)
        self.status_label = QLabel(
            "Ready. Generated STL, MSH, reports and solver validation data stay outside the source tree."
        )
        self.status_label.setWordWrap(True)
        root.addWidget(self.status_label)
        self.output_view = QPlainTextEdit()
        self.output_view.setReadOnly(True)
        self.output_view.setMaximumBlockCount(12000)
        self.output_view.setPlaceholderText("Cutting, Stage A, Gmsh and OneFLOW output will appear here.")
        root.addWidget(self.output_view, 1)
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.accept)
        root.addWidget(close_button)
        for widget in (
            self.radius, self.axis_start, self.axis_end, self.cut_start, self.cut_end,
            self.circumferential, self.axial, self.mesh_size, self.gmsh_path, self.output_directory,
        ):
            widget.textChanged.connect(self._refresh_command)
        self._refresh_command()

    def _mesh_path(self) -> Path:
        return Path(self.output_directory.text()).expanduser() / "cylinder_fluid_domain.msh"

    def _arguments(self) -> list[str]:
        output = Path(self.output_directory.text().strip()).expanduser()
        if not self.output_directory.text().strip():
            raise ValueError("Choose an output directory.")
        radius = float(self.radius.text())
        axis_start, axis_end = float(self.axis_start.text()), float(self.axis_end.text())
        cut_start, cut_end = float(self.cut_start.text()), float(self.cut_end.text())
        mesh_size = float(self.mesh_size.text())
        circumferential, axial = int(self.circumferential.text()), int(self.axial.text())
        numeric = (radius, axis_start, axis_end, cut_start, cut_end, mesh_size)
        if not all(math.isfinite(value) for value in numeric):
            raise ValueError("Cutting and mesh parameters must be finite numbers.")
        if radius <= 0 or not axis_start < cut_start < cut_end < axis_end:
            raise ValueError("Require radius > 0 and axis_start < cut_start < cut_end < axis_end.")
        try:
            output.resolve().relative_to(self._repository_root)
        except ValueError:
            pass
        else:
            raise ValueError("Generated files must stay outside the source repository.")
        if circumferential < 3 or axial < 1:
            raise ValueError("Use at least 3 circumferential and 1 axial division.")
        if mesh_size <= 0:
            raise ValueError("Gmsh target mesh size must be positive.")
        args = [
            str(self._repository_root / "grid" / "stl_to_gmsh" / "generate_cylinder_volume_mesh.py"),
            "--output-directory", str(output),
            "--radius", self.radius.text().strip(),
            "--axis-start", self.axis_start.text().strip(),
            "--axis-end", self.axis_end.text().strip(),
            "--cut-start", self.cut_start.text().strip(),
            "--cut-end", self.cut_end.text().strip(),
            "--circumferential-divisions", str(circumferential),
            "--axial-divisions", str(axial),
            "--mesh-size", self.mesh_size.text().strip(),
        ]
        gmsh = self.gmsh_path.text().strip()
        if gmsh:
            args.extend(["--gmsh", gmsh])
        return args

    def _command_text(self) -> str:
        command = [sys.executable, *self._arguments()]
        return subprocess.list2cmdline(command) if os.name == "nt" else shlex.join(command)

    def _refresh_command(self, *_args) -> None:
        try:
            self.command_view.setPlainText(self._command_text())
            idle = self._process.state() == QProcess.ProcessState.NotRunning
            self.generate_button.setEnabled(idle)
            self.status_label.setText("Parameters are valid. Ready to generate.")
        except (ValueError, OSError) as error:
            self.command_view.setPlainText(f"Command cannot be prepared:\n{error}")
            self.generate_button.setEnabled(False)
            self.status_label.setText("Check the cutting and mesh parameters.")

    def _browse_gmsh(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Select Gmsh executable", str(Path.home()),
            "Gmsh executable (gmsh.exe gmsh);;All files (*)",
        )
        if path:
            self.gmsh_path.setText(path)

    def _browse_output(self) -> None:
        path = QFileDialog.getExistingDirectory(
            self, "Select output directory", self.output_directory.text()
        )
        if path:
            self.output_directory.setText(path)

    def _generate(self) -> None:
        try:
            args = self._arguments()
        except (ValueError, OSError) as error:
            QMessageBox.warning(self, "Invalid cylinder cutting parameters", str(error))
            return
        self.output_view.clear()
        self.run_oneflow_button.setEnabled(False)
        self.generate_button.setEnabled(False)
        self.status_label.setText(
            "Generating analytic cut surface, checking boundary loops and running Gmsh..."
        )
        self._process.start(sys.executable, args)

    def _read_output(self) -> None:
        data = bytes(self._process.readAllStandardOutput()).decode("utf-8", errors="replace")
        if data:
            self.output_view.moveCursor(QTextCursor.MoveOperation.End)
            self.output_view.insertPlainText(data)
            self.output_view.ensureCursorVisible()

    def _process_finished(self, exit_code: int, _exit_status) -> None:
        self._read_output()
        self.generate_button.setEnabled(True)
        if exit_code == 0 and self._mesh_path().is_file():
            self.run_oneflow_button.setEnabled(True)
            self.status_label.setText(
                "Cut surface and Gmsh volume mesh generated. Stage A boundary-loop and Physical Group checks passed."
            )
        else:
            self.run_oneflow_button.setEnabled(False)
            self.status_label.setText(f"Generation failed (exit code {exit_code}). See the execution output.")

    def _process_error(self, _error) -> None:
        self._read_output()
        self.generate_button.setEnabled(True)
        self.run_oneflow_button.setEnabled(False)
        self.status_label.setText("Could not start the mesh-generation process. Check Python and Gmsh paths.")

    def _run_oneflow(self) -> None:
        mesh_path = self._mesh_path().expanduser().resolve()
        if not mesh_path.is_file():
            QMessageBox.warning(self, "Volume mesh required", "Generate the Gmsh volume mesh first.")
            return
        from .mesh_validation_dialog import MeshValidationDialog
        dialog = MeshValidationDialog(
            self._repository_root,
            mesh_file=mesh_path,
            project_directory=mesh_path.parent,
            mesh_case_id=CYLINDER_MESH_CASE_ID,
            parent=self,
        )
        dialog.exec()
