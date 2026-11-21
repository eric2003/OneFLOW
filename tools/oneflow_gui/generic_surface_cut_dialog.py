"""Workbench page for cutting arbitrary closed STL shells by two planes."""
from __future__ import annotations

import os
import shlex
import subprocess
import sys
from pathlib import Path

from PySide6.QtCore import QProcess
from PySide6.QtGui import QTextCursor
from PySide6.QtWidgets import (
    QDialog, QFileDialog, QFormLayout, QHBoxLayout, QLabel, QLineEdit,
    QMessageBox, QPlainTextEdit, QPushButton, QVBoxLayout,
)
from grid.stl_to_gmsh.generate_volume_mesh import DEFAULT_SPHERE_STL
from .user_data_paths import UserDataPaths


class GenericSurfaceCutDialog(QDialog):
    """Cut a user-selected STL without consulting analytic geometry metadata."""

    def __init__(self, repository_root: Path, parent=None) -> None:
        super().__init__(parent)
        self._repository_root = Path(repository_root).resolve()
        self._process = QProcess(self)
        self._process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self._process.setWorkingDirectory(str(self._repository_root))
        self._process.readyReadStandardOutput.connect(self._read_output)
        self._process.finished.connect(self._finished)
        self._process.errorOccurred.connect(self._process_error)
        self.setWindowTitle("Generic STL Surface Cutting")
        self.resize(900, 760)
        root = QVBoxLayout(self)
        heading = QLabel("Generic STL Surface Cutting")
        heading.setStyleSheet("font-size: 22px; font-weight: 650;")
        root.addWidget(heading)
        description = QLabel(
            "Treat the STL as an unknown triangulated surface. No cylinder, pipe, "
            "radius, centerline or analytic geometry definition is used. Define two "
            "cutting planes by a point and normal; the retained shell is capped at "
            "both ends. The first version requires one simple closed loop per plane."
        )
        description.setWordWrap(True)
        root.addWidget(description)

        form = QFormLayout()
        self.input_stl = QLineEdit()
        cylinder = (
            self._repository_root / "grid" / "stl_to_gmsh" / "3d"
            / "cylinder_internal_flow" / "automatic_tetrahedra" / "python"
            / "input" / "model.stl"
        )
        self.input_stl.setText(str(cylinder if cylinder.is_file() else DEFAULT_SPHERE_STL))
        input_row = QHBoxLayout()
        input_row.addWidget(self.input_stl, 1)
        browse = QPushButton("Browse...")
        browse.clicked.connect(self._browse_input)
        input_row.addWidget(browse)
        form.addRow("Input STL", input_row)

        self.inlet_point = self._line("-1.5, 0, 0")
        self.inlet_normal = self._line("1, 0, 0")
        self.outlet_point = self._line("1.5, 0, 0")
        self.outlet_normal = self._line("1, 0, 0")
        self.tolerance = self._line("1e-8")
        form.addRow("Inlet plane point (x, y, z)", self.inlet_point)
        form.addRow("Inlet plane normal", self.inlet_normal)
        form.addRow("Outlet plane point (x, y, z)", self.outlet_point)
        form.addRow("Outlet plane normal", self.outlet_normal)
        form.addRow("Geometric tolerance", self.tolerance)

        default_output = UserDataPaths.discover().mesh_generation / "generic_surface_cut" / Path(self.input_stl.text()).stem
        self.output_directory = self._line(str(default_output))
        output_row = QHBoxLayout()
        output_row.addWidget(self.output_directory, 1)
        output_browse = QPushButton("Browse...")
        output_browse.clicked.connect(self._browse_output)
        output_row.addWidget(output_browse)
        form.addRow("Output directory", output_row)
        root.addLayout(form)

        root.addWidget(QLabel(
            "Outputs: wall.stl, inlet.stl, outlet.stl, closed_surface.stl and cut_report.json"
        ))
        root.addWidget(QLabel("Command preview"))
        self.command_view = QPlainTextEdit()
        self.command_view.setReadOnly(True)
        self.command_view.setMaximumHeight(105)
        root.addWidget(self.command_view)
        actions = QHBoxLayout()
        self.cut_button = QPushButton("Cut STL Surface")
        self.cut_button.clicked.connect(self._run_cut)
        actions.addWidget(self.cut_button)
        self.open_gmsh_button = QPushButton("Open Closed Surface in STL to Gmsh...")
        self.open_gmsh_button.setEnabled(False)
        self.open_gmsh_button.clicked.connect(self._open_gmsh)
        actions.addWidget(self.open_gmsh_button)
        root.addLayout(actions)
        self.status = QLabel("Ready. Output will be written outside the source tree.")
        self.status.setWordWrap(True)
        root.addWidget(self.status)
        self.output_view = QPlainTextEdit()
        self.output_view.setReadOnly(True)
        self.output_view.setMaximumBlockCount(10000)
        root.addWidget(self.output_view, 1)
        close = QPushButton("Close")
        close.clicked.connect(self.accept)
        root.addWidget(close)
        for widget in (
            self.input_stl, self.inlet_point, self.inlet_normal,
            self.outlet_point, self.outlet_normal, self.tolerance,
            self.output_directory,
        ):
            widget.textChanged.connect(self._refresh)
        self._refresh()

    @staticmethod
    def _line(value: str) -> QLineEdit:
        return QLineEdit(value)

    @staticmethod
    def _parse_vector(value: str, label: str) -> tuple[float, float, float]:
        try:
            parts = [float(item.strip()) for item in value.replace(";", ",").split(",")]
        except ValueError as error:
            raise ValueError(f"{label} must contain three comma-separated numbers.") from error
        if len(parts) != 3:
            raise ValueError(f"{label} must contain exactly three coordinates.")
        return tuple(parts)

    def _arguments(self) -> list[str]:
        source = Path(self.input_stl.text().strip()).expanduser()
        output = Path(self.output_directory.text().strip()).expanduser()
        if not source.is_file():
            raise ValueError("Select an existing STL file.")
        if not self.output_directory.text().strip():
            raise ValueError("Choose an output directory.")
        try:
            output.resolve().relative_to(self._repository_root)
        except ValueError:
            pass
        else:
            raise ValueError("Generated files must stay outside the source repository.")
        inlet_point = self._parse_vector(self.inlet_point.text(), "Inlet point")
        inlet_normal = self._parse_vector(self.inlet_normal.text(), "Inlet normal")
        outlet_point = self._parse_vector(self.outlet_point.text(), "Outlet point")
        outlet_normal = self._parse_vector(self.outlet_normal.text(), "Outlet normal")
        tolerance = float(self.tolerance.text())
        if tolerance <= 0:
            raise ValueError("Tolerance must be positive.")
        script_args = [
            "-m", "grid.stl_import.cut_surface_cli",
            "--input", str(source.resolve()),
            "--output-directory", str(output.resolve()),
            "--inlet-point", *map(str, inlet_point),
            "--inlet-normal", *map(str, inlet_normal),
            "--outlet-point", *map(str, outlet_point),
            "--outlet-normal", *map(str, outlet_normal),
            "--tolerance", str(tolerance),
        ]
        return script_args

    def _command(self) -> str:
        command = [sys.executable, *self._arguments()]
        return subprocess.list2cmdline(command) if os.name == "nt" else shlex.join(command)

    def _refresh(self, *_args) -> None:
        try:
            self.command_view.setPlainText(self._command())
            self.cut_button.setEnabled(self._process.state() == QProcess.ProcessState.NotRunning)
            self.status.setText("Parameters are valid. Ready to cut.")
        except (OSError, ValueError) as error:
            self.command_view.setPlainText(f"Command cannot be prepared:\n{error}")
            self.cut_button.setEnabled(False)
            self.status.setText("Check the STL path and cutting parameters.")

    def _browse_input(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Select input STL", str(Path.home()), "STL files (*.stl);;All files (*)"
        )
        if path:
            self.input_stl.setText(path)

    def _browse_output(self) -> None:
        path = QFileDialog.getExistingDirectory(
            self, "Select output directory", self.output_directory.text()
        )
        if path:
            self.output_directory.setText(path)

    def _run_cut(self) -> None:
        try:
            args = self._arguments()
        except (OSError, ValueError) as error:
            QMessageBox.warning(self, "Invalid generic cutting parameters", str(error))
            return
        self.output_view.clear()
        self.open_gmsh_button.setEnabled(False)
        self.cut_button.setEnabled(False)
        self.status.setText("Reading STL triangles, clipping both planes and capping the loops...")
        self._process.start(sys.executable, args)

    def _read_output(self) -> None:
        data = bytes(self._process.readAllStandardOutput()).decode("utf-8", errors="replace")
        if data:
            self.output_view.moveCursor(QTextCursor.MoveOperation.End)
            self.output_view.insertPlainText(data)
            self.output_view.ensureCursorVisible()

    def _finished(self, exit_code: int, _exit_status) -> None:
        self._read_output()
        self.cut_button.setEnabled(True)
        closed = Path(self.output_directory.text()).expanduser() / "closed_surface.stl"
        if exit_code == 0 and closed.is_file():
            self.open_gmsh_button.setEnabled(True)
            self.status.setText("Generic STL cutting succeeded. Four STL artifacts and the report were written.")
        else:
            self.open_gmsh_button.setEnabled(False)
            self.status.setText(f"Generic STL cutting failed (exit code {exit_code}).")

    def _process_error(self, _error) -> None:
        self._read_output()
        self.cut_button.setEnabled(True)
        self.open_gmsh_button.setEnabled(False)
        self.status.setText("Could not start the Python cutting process.")

    def _open_gmsh(self) -> None:
        closed = Path(self.output_directory.text()).expanduser() / "closed_surface.stl"
        if not closed.is_file():
            QMessageBox.warning(self, "Cut surface required", "Run generic STL cutting first.")
            return
        from .mesh_project import create_mesh_project_from_stl
        from .mesh_generation_dialog import StlToGmshDialog
        try:
            project = create_mesh_project_from_stl(self._repository_root, closed)
        except (OSError, ValueError, KeyError) as error:
            QMessageBox.warning(self, "Cannot create mesh project", str(error))
            return
        StlToGmshDialog(self._repository_root, self, project=project).exec()
