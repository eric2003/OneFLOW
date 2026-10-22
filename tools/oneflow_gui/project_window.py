"""Open and run a user-selected OneFLOW project directory."""

from __future__ import annotations

import shutil
import threading
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from PySide6.QtCore import QThread, Signal
from PySide6.QtWidgets import (
    QComboBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from .executable_discovery import executable_candidates
from .project_repository import ProjectDefinition, ProjectType
from .user_data_paths import UserDataPaths
from .workflow_process import run_process
from .workflow_repository import WorkflowRepository, WorkflowRunRecord
from .runner import create_run_id


def _project_workflow_id(name: str) -> str:
    """Return a history identifier that satisfies the workflow ID schema."""
    slug = "".join(
        character.lower() if character.isascii() and character.isalnum() else "-"
        for character in name
    )
    slug = "-".join(part for part in slug.split("-") if part) or "project"
    return f"project-{slug[:48]}"


@dataclass(frozen=True, slots=True)
class ProjectRunOutcome:
    """Serializable outcome of one project run."""

    run_id: str
    status: str
    message: str
    project_directory: Path
    return_code: int | None


class ProjectRunWorker(QThread):
    """Stage a project when requested, then execute it away from the UI thread."""

    log_line = Signal(str)
    completed = Signal(object)

    def __init__(
        self,
        project: ProjectDefinition,
        executable: Path,
        run_id: str,
        execution_mode: str,
        custom_root: Path | None,
        repository_root: Path,
    ) -> None:
        super().__init__()
        self._project = project
        self._executable = executable
        self._repository_root = repository_root.resolve()
        self._run_id = run_id
        self._execution_mode = execution_mode
        self._custom_root = custom_root
        self._cancel = threading.Event()

    def request_stop(self) -> None:
        self._cancel.set()

    def run(self) -> None:
        project_directory = self._project.root
        try:
            source_directory = self._project.root.resolve()
            repository_root = self._repository_root.resolve()
            is_repository_content = (
                source_directory == repository_root
                or repository_root in source_directory.parents
            )
            if self._execution_mode == "source" and is_repository_content:
                raise ValueError(
                    "Projects inside the OneFLOW source repository cannot run in place. "
                    "Copy the project to a workspace first."
                )
            if self._execution_mode != "source":
                if self._execution_mode == "default_workspace":
                    destination_root = UserDataPaths.discover().executions / self._run_id
                else:
                    if self._custom_root is None:
                        raise ValueError("Choose a target directory for the project copy.")
                    destination_root = self._custom_root
                destination_root = destination_root.resolve()
                if (
                    destination_root == self._repository_root
                    or self._repository_root in destination_root.parents
                ):
                    raise ValueError(
                        "Project copies must be placed outside the OneFLOW source repository."
                    )
                project_directory = (destination_root / self._project.name).resolve()
                source_directory = self._project.root.resolve()
                if (
                    project_directory == source_directory
                    or source_directory in project_directory.parents
                ):
                    raise ValueError(
                        "The project copy destination must be outside the source project."
                    )
                destination_root.mkdir(parents=True, exist_ok=True)
                if project_directory.exists():
                    raise FileExistsError(
                        f"Project copy already exists: {project_directory}"
                    )
                shutil.copytree(
                    self._project.root,
                    project_directory,
                    ignore=shutil.ignore_patterns(".git", "__pycache__"),
                )
                self.log_line.emit(f"Project copied to: {project_directory}")

            command = [str(self._executable.resolve()), "r", str(project_directory)]
            self.log_line.emit("Command arguments:")
            for index, argument in enumerate(command):
                self.log_line.emit(f"  [{index}] {argument}")
            result = run_process(
                command,
                cwd=self._executable.resolve().parent,
                on_output=self.log_line.emit,
                is_cancelled=self._cancel.is_set,
            )
            self.completed.emit(
                ProjectRunOutcome(
                    self._run_id,
                    result.status,
                    result.message,
                    project_directory,
                    result.return_code,
                )
            )
        except (OSError, shutil.Error, ValueError) as error:
            self.completed.emit(
                ProjectRunOutcome(
                    self._run_id, "failed", str(error), project_directory, None
                )
            )


class ProjectWindow(QMainWindow):
    """Project-oriented entry point independent of catalog registration."""

    def __init__(
        self,
        project: ProjectDefinition,
        repository_root: Path,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self._project = project
        self._repository_root = repository_root.resolve()
        self._worker: ProjectRunWorker | None = None
        self._active_record: WorkflowRunRecord | None = None
        self._log_path: Path | None = None
        paths = UserDataPaths.discover()
        self._history = WorkflowRepository(paths.workflow_data)
        self.setWindowTitle(f"OneFLOW Project - {project.name}")
        self.resize(900, 650)
        self._build_ui()

    def _build_ui(self) -> None:
        central = QWidget(self)
        root = QVBoxLayout(central)
        root.setContentsMargins(22, 20, 22, 20)
        root.setSpacing(12)

        title = QLabel(self._project.name)
        title.setStyleSheet("font-size: 24px; font-weight: 650;")
        root.addWidget(title)
        root.addWidget(QLabel(str(self._project.root)))

        project_type = self._project.project_type
        type_label = {
            ProjectType.GRID: "Grid project",
            ProjectType.FLOW: "CFD / flow project",
            ProjectType.UNKNOWN: "Unknown project type (directory opened without catalog)",
        }[project_type]
        details = QFormLayout()
        details.addRow("Project type", QLabel(type_label))
        details.addRow(
            "Control file",
            QLabel(str(self._project.control_file) if self._project.control_file else "Not detected"),
        )
        details.addRow(
            "Runtime file",
            QLabel(str(self._project.runtime_file) if self._project.runtime_file else "Not detected"),
        )
        root.addLayout(details)

        executable_row = QHBoxLayout()
        executable_row.addWidget(QLabel("OneFLOW executable"))
        self.executable_combo = QComboBox()
        self.executable_combo.setEditable(True)
        self.executable_combo.setMinimumContentsLength(42)
        self.executable_combo.addItems(
            [str(path) for path in executable_candidates(self._repository_root)]
        )
        if self.executable_combo.count():
            self.executable_combo.setCurrentIndex(0)
        executable_row.addWidget(self.executable_combo, 1)
        browse_executable = QPushButton("Browse...")
        browse_executable.clicked.connect(self._choose_executable)
        executable_row.addWidget(browse_executable)
        root.addLayout(executable_row)

        location_row = QHBoxLayout()
        location_row.addWidget(QLabel("Execution location"))
        self.location_combo = QComboBox()
        source_is_repository_content = (
            self._project.root.resolve() == self._repository_root
            or self._repository_root in self._project.root.resolve().parents
        )
        if not source_is_repository_content:
            self.location_combo.addItem("Run in original project directory", "source")
        self.location_combo.addItem("Copy to OneFLOW default workspace", "default_workspace")
        self.location_combo.addItem("Copy to specified directory", "custom_workspace")
        self.location_combo.currentIndexChanged.connect(self._location_changed)
        location_row.addWidget(self.location_combo, 1)
        root.addLayout(location_row)

        target_row = QHBoxLayout()
        target_row.addWidget(QLabel("Copy destination"))
        self.target_edit = QLineEdit()
        self.target_edit.setPlaceholderText("Choose a parent directory for the project copy")
        target_row.addWidget(self.target_edit, 1)
        browse_target = QPushButton("Browse...")
        browse_target.clicked.connect(self._choose_target)
        target_row.addWidget(browse_target)
        root.addLayout(target_row)
        self._browse_target_button = browse_target
        self._location_changed()

        self.log_view = QPlainTextEdit()
        self.log_view.setReadOnly(True)
        root.addWidget(self.log_view, 1)

        actions = QHBoxLayout()
        self.run_button = QPushButton("Run Project")
        self.run_button.clicked.connect(self._start_run)
        self.stop_button = QPushButton("Stop")
        self.stop_button.setEnabled(False)
        self.stop_button.clicked.connect(self._stop_run)
        actions.addWidget(self.run_button)
        actions.addWidget(self.stop_button)
        actions.addStretch(1)
        root.addLayout(actions)
        self.setCentralWidget(central)

    def _choose_executable(self) -> None:
        current = self.executable_combo.currentText().strip()
        current_path = Path(current).expanduser() if current else None
        start = (
            str(current_path.parent)
            if current_path is not None and current_path.parent.is_dir()
            else str(self._repository_root)
        )
        path, _ = QFileDialog.getOpenFileName(
            self, "Select OneFLOW executable", start,
            "OneFLOW executable (OneFLOW.exe);;Executable (*.exe);;All files (*)",
        )
        if path:
            self.executable_combo.setEditText(path)

    def _choose_target(self) -> None:
        path = QFileDialog.getExistingDirectory(
            self, "Select project-copy destination", self.target_edit.text()
        )
        if path:
            self.target_edit.setText(path)

    def _location_changed(self, *_args) -> None:
        custom = self.location_combo.currentData() == "custom_workspace"
        self.target_edit.setEnabled(custom)
        self._browse_target_button.setEnabled(custom)

    def _start_run(self) -> None:
        if self._worker is not None and self._worker.isRunning():
            return
        executable = Path(self.executable_combo.currentText().strip()).expanduser()
        if not executable.is_file():
            QMessageBox.warning(self, "Executable required", "Select a valid OneFLOW executable.")
            return

        mode = self.location_combo.currentData()
        custom_root = None
        if mode == "source":
            answer = QMessageBox.question(
                self,
                "Run in original project directory",
                "OneFLOW may create or overwrite files inside this project. Continue?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No,
            )
            if answer != QMessageBox.StandardButton.Yes:
                return
        elif mode == "custom_workspace":
            target = self.target_edit.text().strip()
            if not target:
                QMessageBox.warning(self, "Destination required", "Choose a copy destination.")
                return
            custom_root = Path(target).expanduser().resolve()

        run_id = create_run_id()
        paths = UserDataPaths.discover()
        log_path = (paths.logs / f"{run_id}.log").resolve()
        try:
            log_path.parent.mkdir(parents=True, exist_ok=True)
            log_path.write_text(
                f"Project run {run_id}\nProject: {self._project.root}\nExecutable: {executable.resolve()}\n",
                encoding="utf-8",
            )
            record = WorkflowRunRecord(
                run_id=run_id,
                workflow_id=_project_workflow_id(self._project.name),
                status="running",
                started_at=datetime.now(timezone.utc).isoformat(timespec="seconds"),
                executable=str(executable.resolve()),
                steps=({
                    "type": "project",
                    "project_name": self._project.name,
                    "project_directory": str(self._project.root),
                    "project_type": self._project.project_type.value,
                    "execution_mode": mode,
                },),
                message="Project run started.",
                log_path=str(log_path),
            )
            self._history.save_run(record)
        except (OSError, ValueError) as error:
            QMessageBox.critical(self, "Cannot record run", str(error))
            return

        self._active_record = record
        self._log_path = log_path
        self.log_view.appendPlainText(f"Run {run_id}\nProject: {self._project.root}")
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self._worker = ProjectRunWorker(
            self._project, executable, run_id, mode, custom_root,
            self._repository_root,
        )
        self._worker.log_line.connect(self._append_log)
        self._worker.completed.connect(self._run_completed)
        self._worker.start()

    def _append_log(self, line: str) -> None:
        self.log_view.appendPlainText(line)
        if self._log_path is not None:
            try:
                with self._log_path.open("a", encoding="utf-8") as stream:
                    stream.write(line + "\n")
            except OSError as error:
                self.statusBar().showMessage(f"Could not write run log: {error}")

    def _run_completed(self, outcome: ProjectRunOutcome) -> None:
        self._append_log(f"{outcome.message}\nExecution directory: {outcome.project_directory}")
        record = self._active_record
        if record is not None:
            updated = WorkflowRunRecord(
                run_id=record.run_id,
                workflow_id=record.workflow_id,
                status=outcome.status,
                started_at=record.started_at,
                finished_at=datetime.now(timezone.utc).isoformat(timespec="seconds"),
                executable=record.executable,
                steps=({
                    "type": "project",
                    "project_name": self._project.name,
                    "project_directory": str(self._project.root),
                    "execution_directory": str(outcome.project_directory),
                    "status": outcome.status,
                    "return_code": outcome.return_code,
                    "message": outcome.message,
                },),
                message=outcome.message,
                log_path=record.log_path,
            )
            try:
                self._history.save_run(updated)
            except (OSError, ValueError) as error:
                self._append_log(f"Could not save run history: {error}")
        self._active_record = None
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.statusBar().showMessage(f"Project run {outcome.status}")
        self._worker = None

    def _stop_run(self) -> None:
        if self._worker is not None:
            self._worker.request_stop()
            self.stop_button.setEnabled(False)
            self.statusBar().showMessage("Stop requested; waiting for the active process to exit.")

    def closeEvent(self, event) -> None:
        if self._worker is not None and self._worker.isRunning():
            self._worker.request_stop()
            event.ignore()
            return
        event.accept()
