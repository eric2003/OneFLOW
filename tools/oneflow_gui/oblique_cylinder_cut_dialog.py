"""GUI runner for the oblique-cut cylinder STL reference case."""
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

from grid.stl_to_gmsh.generate_volume_mesh import discover_gmsh
from .user_data_paths import UserDataPaths


class ObliqueCylinderCutDialog(QDialog):
    """Run the complete blind-STL oblique cutting and Gmsh meshing pipeline."""

    def __init__(self, repository_root: Path, parent=None) -> None:
        super().__init__(parent)
        self._repository_root = Path(repository_root).resolve()
        self._runner = (self._repository_root / "grid" / "stl_to_gmsh" / "3d" / "cylinder_internal_flow" / "oblique_cut" / "python" / "run_cylinder_oblique_cut.py")
        self._case = (
            self._repository_root / "grid" / "stl_to_gmsh" / "3d"
            / "cylinder_internal_flow" / "oblique_cut" / "python"
        )
        self._process = QProcess(self)
        self._process.setWorkingDirectory(str(self._repository_root))
        self._process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self._process.readyReadStandardOutput.connect(self._read_output)
        self._process.finished.connect(self._finished)
        self._process.errorOccurred.connect(self._process_error)

        self.setWindowTitle("Oblique-cut Cylinder — STL to Gmsh")
        self.resize(900, 700)
        self.setMinimumSize(680, 560)
        root = QVBoxLayout(self)
        heading = QLabel("Oblique-cut Cylinder — Complete STL Workflow")
        heading.setStyleSheet("font-size: 22px; font-weight: 650;")
        root.addWidget(heading)
        description = QLabel(
            "Run the registered STL reference case end-to-end: read the triangulated "
            "cylinder STL as an unknown surface, cut it with the configured oblique inlet "
            "and normal outlet planes, cap both loops, generate a Gmsh tetrahedral volume "
            "mesh, and verify the wall/inflow/outflow/fluid_domain Physical Groups. "
            "The reference case uses a 30° inlet cut; the cutter does not use an analytic cylinder."
        )
        description.setWordWrap(True)
        root.addWidget(description)

        form = QFormLayout()
        self.case_directory = QLineEdit(str(self._case))
        case_row = QHBoxLayout()
        case_row.addWidget(self.case_directory, 1)
        case_browse = QPushButton("Browse...")
        case_browse.clicked.connect(self._browse_case)
        case_row.addWidget(case_browse)
        form.addRow("Reference case", case_row)

        default_output = UserDataPaths.discover().mesh_generation / "cylinder_oblique_cut"
        self.output_directory = QLineEdit(str(default_output))
        output_row = QHBoxLayout()
        output_row.addWidget(self.output_directory, 1)
        output_browse = QPushButton("Browse...")
        output_browse.clicked.connect(self._browse_output)
        output_row.addWidget(output_browse)
        form.addRow("Output directory", output_row)

        self.gmsh_path = QLineEdit(discover_gmsh() or "")
        gmsh_row = QHBoxLayout()
        gmsh_row.addWidget(self.gmsh_path, 1)
        gmsh_browse = QPushButton("Browse...")
        gmsh_browse.clicked.connect(self._browse_gmsh)
        gmsh_row.addWidget(gmsh_browse)
        form.addRow("Gmsh executable", gmsh_row)
        root.addLayout(form)

        root.addWidget(QLabel("Command Preview"))
        self.command_view = QPlainTextEdit()
        self.command_view.setReadOnly(True)
        self.command_view.setMaximumHeight(90)
        root.addWidget(self.command_view)

        actions = QHBoxLayout()
        self.run_button = QPushButton("Run Cutting + Gmsh")
        self.run_button.clicked.connect(self._run)
        actions.addWidget(self.run_button)
        self.open_output_button = QPushButton("Open Output Folder")
        self.open_output_button.setEnabled(False)
        self.open_output_button.clicked.connect(self._open_output)
        actions.addWidget(self.open_output_button)
        self.run_cfd_button = QPushButton("Convert + Run OneFLOW...")
        self.run_cfd_button.setEnabled(False)
        self.run_cfd_button.setToolTip(
            "Reuse the generated MSH and run Gmsh → CGNS → native OneFLOW grid → CFD smoke validation."
        )
        self.run_cfd_button.clicked.connect(self._run_mesh_to_cfd)
        actions.addWidget(self.run_cfd_button)
        root.addLayout(actions)

        self.status = QLabel("Ready. All generated files stay outside the source tree.")
        self.status.setWordWrap(True)
        root.addWidget(self.status)
        self.output_view = QPlainTextEdit()
        self.output_view.setReadOnly(True)
        self.output_view.setMaximumBlockCount(12000)
        self.output_view.setPlaceholderText("STL cutting, Gmsh and Physical Group verification output...")
        root.addWidget(self.output_view, 1)
        close = QPushButton("Close")
        close.clicked.connect(self.accept)
        root.addWidget(close)

        for widget in (self.case_directory, self.output_directory, self.gmsh_path):
            widget.textChanged.connect(self._refresh_command)
        self._refresh_command()

    def _arguments(self) -> list[str]:
        runner = self._runner
        case = Path(self.case_directory.text().strip()).expanduser()
        output_text = self.output_directory.text().strip()
        if not runner.is_file():
            raise ValueError(f"Runner script was not found: {runner}")
        if not case.is_dir():
            raise ValueError("Select an existing oblique-cut reference case directory.")
        if not (case.parent / "script" / "cutting.json").is_file() or not (case / "input" / "model.stl").is_file():
            raise ValueError("The selected case must contain script/cutting.json and input/model.stl.")
        if not output_text:
            raise ValueError("Choose an output directory.")
        output = Path(output_text).expanduser().resolve()
        try:
            output.relative_to(self._repository_root)
        except ValueError:
            pass
        else:
            raise ValueError("Generated files must stay outside the source repository.")
        gmsh = self.gmsh_path.text().strip()
        if not gmsh:
            raise ValueError("Gmsh was not detected. Enter its executable path.")
        return [
            str(runner), "--case-directory", str(case.resolve()),
            "--output-directory", str(output), "--gmsh", gmsh,
        ]

    def _command(self) -> str:
        command = [sys.executable, *self._arguments()]
        return subprocess.list2cmdline(command) if os.name == "nt" else shlex.join(command)

    def _refresh_command(self, *_args) -> None:
        try:
            self.command_view.setPlainText(self._command())
            idle = self._process.state() == QProcess.ProcessState.NotRunning
            self.run_button.setEnabled(idle)
            self.status.setText("Parameters are valid. Ready to run the complete workflow.")
        except (OSError, ValueError) as error:
            self.command_view.setPlainText(f"Command cannot be prepared:\n{error}")
            self.run_button.setEnabled(False)
            self.status.setText("Check the reference case, output directory and Gmsh executable.")

    def _browse_case(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Select oblique-cut reference case", self.case_directory.text())
        if path:
            self.case_directory.setText(path)

    def _browse_output(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Select output directory", self.output_directory.text())
        if path:
            self.output_directory.setText(path)

    def _browse_gmsh(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Select Gmsh executable", str(Path.home()),
            "Gmsh executable (gmsh.exe gmsh);;All files (*)",
        )
        if path:
            self.gmsh_path.setText(path)

    def _run(self) -> None:
        try:
            args = self._arguments()
        except (OSError, ValueError) as error:
            QMessageBox.warning(self, "Invalid oblique-cut workflow settings", str(error))
            return
        self.output_view.clear()
        self.open_output_button.setEnabled(False)
        self.run_cfd_button.setEnabled(False)
        self.run_button.setEnabled(False)
        self.status.setText("Running STL cutting, Gmsh volume meshing and Physical Group checks...")
        self._process.start(sys.executable, args)

    def _read_output(self) -> None:
        data = bytes(self._process.readAllStandardOutput()).decode("utf-8", errors="replace")
        if data:
            self.output_view.moveCursor(QTextCursor.MoveOperation.End)
            self.output_view.insertPlainText(data)
            self.output_view.ensureCursorVisible()

    def _finished(self, exit_code: int, _exit_status) -> None:
        self._read_output()
        self.run_button.setEnabled(True)
        mesh_ready = exit_code == 0 and self._mesh_path().is_file()
        self.open_output_button.setEnabled(mesh_ready)
        self.run_cfd_button.setEnabled(mesh_ready)
        if mesh_ready:
            self.status.setText("PASS: oblique STL cutting and Gmsh mesh complete. You can now convert and run OneFLOW.")
        else:
            self.status.setText(f"FAIL: workflow exited with code {exit_code} or the expected MSH is missing. See execution output.")

    def _process_error(self, _error) -> None:
        self._read_output()
        self.run_button.setEnabled(True)
        self.run_cfd_button.setEnabled(False)
        self.status.setText("Could not start Python. Check the active Python environment.")

    def _mesh_path(self) -> Path:
        return Path(self.output_directory.text()).expanduser() / "results" / "cylinder_oblique_cut.msh"

    def _run_mesh_to_cfd(self) -> None:
        mesh_path = self._mesh_path().expanduser().resolve()
        if not mesh_path.is_file():
            QMessageBox.warning(self, "Volume mesh required", "Generate the oblique-cut Gmsh volume mesh first.")
            return
        from .mesh_validation_dialog import MeshValidationDialog
        dialog = MeshValidationDialog(
            self._repository_root,
            mesh_file=mesh_path,
            project_directory=mesh_path.parent.parent,
            mesh_case_id="stl_to_gmsh_oblique_cylinder_cut",
            parent=self,
        )
        dialog.exec()

    def _open_output(self) -> None:
        output = Path(self.output_directory.text()).expanduser().resolve()
        if output.is_dir():
            QProcess.startDetached("explorer" if os.name == "nt" else "xdg-open", [str(output)])
