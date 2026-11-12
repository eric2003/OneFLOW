"""Unified mesh-generation entry point and direct STL-to-Gmsh interface."""

from __future__ import annotations

import math
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

from PySide6.QtCore import QProcess, QTextCursor
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QVBoxLayout,
)

from .contract_workbench_dialog import ContractDomainDialog
from .mesh_layout_dialog import MeshLayoutDialog
from .user_data_paths import UserDataPaths


class MeshGenerationDialog(QDialog):
    """Navigation hub for mesh generation and conversion activities."""

    def __init__(self, repository_root: Path, parent=None) -> None:
        super().__init__(parent)
        self._repository_root = repository_root.resolve()
        self.setWindowTitle("Mesh Generation & Conversion")
        self.resize(620, 420)

        root = QVBoxLayout(self)
        title = QLabel("Mesh Generation & Conversion")
        title.setStyleSheet("font-size: 22px; font-weight: 650;")
        root.addWidget(title)
        intro = QLabel(
            "Generate meshes with the existing OneFLOW layout tools or external "
            "mesh engines. Each operation exposes its parameters, command and results."
        )
        intro.setWordWrap(True)
        root.addWidget(intro)

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

    def _open_mesh_layout(self) -> None:
        MeshLayoutDialog(self._repository_root, self).exec()

    def _open_stl_to_gmsh(self) -> None:
        StlToGmshDialog(self._repository_root, self).exec()

    def _open_mesh_knowledge(self) -> None:
        ContractDomainDialog(self._repository_root, "mesh", self).exec()


class StlToGmshDialog(QDialog):
    """Configure, preview and execute the direct closed-STL volume mesher."""

    def __init__(self, repository_root: Path, parent=None) -> None:
        super().__init__(parent)
        self._repository_root = repository_root.resolve()
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

        form = QFormLayout()
        self.geometry_preset = QComboBox()
        self.geometry_preset.addItem("Unit Sphere STL", "sphere")
        self.geometry_preset.addItem("General STL", "general")
        self.geometry_preset.setToolTip(
            "The preset describes the geometry workflow; select the actual STL file below."
        )
        form.addRow("Geometry preset", self.geometry_preset)

        stl_row = QHBoxLayout()
        self.stl_path = QLineEdit()
        default_sphere_stl = (
            self._repository_root / "grid" / "stl_to_gmsh" / "3d" / "sphere"
            / "automatic_tetrahedra" / "python" / "input" / "model.stl"
        )
        if default_sphere_stl.is_file():
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

        self.mesh_size = QLineEdit("0.20")
        self.mesh_size.setToolTip("One target size; maps to both MeshSizeMin and MeshSizeMax.")
        form.addRow("Target size", self.mesh_size)

        self.minimum_size = QLineEdit("0.10")
        self.maximum_size = QLineEdit("0.20")
        size_row = QHBoxLayout()
        size_row.addWidget(QLabel("Minimum"))
        size_row.addWidget(self.minimum_size)
        size_row.addWidget(QLabel("Maximum"))
        size_row.addWidget(self.maximum_size)
        form.addRow("Advanced size bounds", size_row)

        self.output_path = QLineEdit(
            str(UserDataPaths.discover().mesh_generation / "sphere.msh")
        )
        output_row = QHBoxLayout()
        output_row.addWidget(self.output_path, 1)
        output_browse = QPushButton("Browse...")
        output_browse.clicked.connect(self._browse_output)
        output_row.addWidget(output_browse)
        form.addRow("Output mesh", output_row)

        self.gmsh_path = QLineEdit(shutil.which("gmsh") or "")
        self.gmsh_path.setPlaceholderText("Optional: leave blank to search PATH")
        gmsh_row = QHBoxLayout()
        gmsh_row.addWidget(self.gmsh_path, 1)
        gmsh_browse = QPushButton("Browse...")
        gmsh_browse.clicked.connect(self._browse_gmsh)
        gmsh_row.addWidget(gmsh_browse)
        form.addRow("Gmsh executable", gmsh_row)

        self.boundary_name = QLineEdit("wall")
        self.volume_name = QLineEdit("fluid_domain")
        groups_row = QHBoxLayout()
        groups_row.addWidget(QLabel("Surface"))
        groups_row.addWidget(self.boundary_name)
        groups_row.addWidget(QLabel("Volume"))
        groups_row.addWidget(self.volume_name)
        form.addRow("Physical Groups", groups_row)
        root.addLayout(form)

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
        self.run_button = QPushButton("Generate Mesh")
        self.run_button.setMinimumWidth(150)
        self.run_button.clicked.connect(self._run)
        actions.addWidget(self.run_button)
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
        self._density_changed()
        self._refresh_command()

    def _browse_stl(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Select STL surface", str(Path.home()), "STL files (*.stl);;All files (*)"
        )
        if path:
            self.stl_path.setText(path)
            if "sphere" in Path(path).name.lower():
                self.geometry_preset.setCurrentIndex(0)
            else:
                self.geometry_preset.setCurrentIndex(1)

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

    def _density_changed(self) -> None:
        preset = self.density.currentData()
        sizes = {"coarse": "0.40", "medium": "0.20", "fine": "0.10"}
        custom = preset == "custom"
        self.mesh_size.setEnabled(not custom)
        self.minimum_size.setEnabled(custom)
        self.maximum_size.setEnabled(custom)
        if not custom:
            self.mesh_size.setText(sizes[preset])
            self.minimum_size.setText(sizes[preset])
            self.maximum_size.setText(sizes[preset])
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
            self.run_button.setEnabled(self._process.state() == QProcess.ProcessState.NotRunning)
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
        except (OSError, ValueError) as error:
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

    def _process_finished(self, exit_code: int, _exit_status) -> None:
        self._read_output()
        self.run_button.setEnabled(True)
        if exit_code == 0:
            mesh_path = Path(self.output_path.text()).expanduser()
            self.status_label.setText(
                "Completed and validated. Outputs: "
                f"{mesh_path}, {mesh_path.with_suffix('.geo')}, "
                f"{mesh_path.with_suffix('.report.json')}"
            )
        else:
            self.status_label.setText(f"Generation failed (exit code {exit_code}). See execution output.")

    def _process_error(self, _error) -> None:
        self._read_output()
        self.run_button.setEnabled(True)
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
