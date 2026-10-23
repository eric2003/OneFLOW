"""Open and run a user-selected OneFLOW project directory."""

from __future__ import annotations

import json
import shutil
import threading
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from PySide6.QtCore import QThread, Signal, Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QDialog,
    QComboBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QSplitter,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QTabWidget,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from .executable_discovery import executable_candidates
from .project_inventory import (
    FileChange,
    ProjectFileEntry,
    compare_inventory,
    scan_project,
    scan_workflow_run_outputs,
)
from .project_repository import ProjectDefinition, ProjectType
from .user_data_paths import UserDataPaths
from .workflow_graph_runner import run_case_graph
from .workflow_plan import build_case_graph_plan
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

    def _run_case_graph(self, project_directory: Path):
        workflow_file = project_directory / "workflow.json"
        try:
            workflow = json.loads(workflow_file.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError) as error:
            raise ValueError(f"Cannot read workflow definition {workflow_file}: {error}") from error

        steps = workflow.get("steps")
        if not isinstance(steps, list):
            raise ValueError("Workflow definition must contain a steps list.")
        graph_steps = [step for step in steps if isinstance(step, dict) and step.get("type") == "case_graph"]
        if len(graph_steps) != 1 or len(steps) != 1:
            raise ValueError("Project workflows must contain exactly one case_graph step.")

        plan = build_case_graph_plan(graph_steps[0].get("cases", []), project_directory)
        run_root = project_directory / ".oneflow-runs" / self._run_id
        self.log_line.emit(f"Workflow definition: {workflow_file}")
        self.log_line.emit(f"Workflow run directory: {run_root}")

        def run_step(step):
            command = [str(self._executable.resolve()), "r", str(step.directory)]
            self.log_line.emit(f"Workflow case: {step.case_id}")
            self.log_line.emit("Command arguments:")
            for index, argument in enumerate(command):
                self.log_line.emit(f"  [{index}] {argument}")
            result = run_process(
                command,
                cwd=self._executable.resolve().parent,
                on_output=self.log_line.emit,
                is_cancelled=self._cancel.is_set,
            )
            self.log_line.emit(f"Case {step.case_id}: {result.status} - {result.message}")
            return result.status, result.message

        in_place = self._execution_mode == "source"
        if in_place:
            self.log_line.emit("Execution mode: run workflow cases in their project directories")
        else:
            self.log_line.emit("Execution mode: run isolated workflow case copies")
        graph_result = run_case_graph(
            plan,
            run_step,
            run_root=run_root,
            is_cancelled=self._cancel.is_set,
            in_place=in_place,
        )
        summary = "; ".join(
            f"{case.case_id}={case.status}" for case in graph_result.cases
        )
        self.log_line.emit(f"Workflow result: {summary}")
        return ProjectRunOutcome(
            self._run_id,
            graph_result.status,
            f"Workflow {graph_result.status}: {summary}",
            project_directory,
            0 if graph_result.status == "succeeded" else (
                None if graph_result.status == "cancelled" else 1
            ),
        )

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

            workflow_file = project_directory / "workflow.json"
            if workflow_file.is_file():
                result = self._run_case_graph(project_directory)
            else:
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
        self._inventory_baseline = scan_project(project.root)
        self._run_changes: list[FileChange] = []
        self._run_output_count = 0
        self._run_change_counts = {"unchanged": 0, "added": 0, "modified": 0, "deleted": 0}
        self._view_windows: dict[str, QMainWindow] = {}
        self._detached_log: QPlainTextEdit | None = None
        self._detached_project_tree: QTreeWidget | None = None
        self._detached_changes_tree: QTreeWidget | None = None
        self._detached_changes_summary: QLabel | None = None
        self._detached_show_unchanged: QCheckBox | None = None
        paths = UserDataPaths.discover()
        self._history = WorkflowRepository(paths.workflow_data)
        self.setWindowTitle(f"OneFLOW Project - {project.name}")
        self.resize(1180, 820)
        self.setMinimumSize(900, 650)
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
        parameter_text = "\n".join(str(path) for path in self._project.parameter_files)
        details.addRow(
            "Parameter files",
            QLabel(parameter_text or "None detected"),
        )

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

        self.inventory_tabs = QTabWidget()
        self.project_tree = QTreeWidget()
        self.project_tree.setHeaderLabels(["Path", "Kind", "Category", "Size", "Purpose"])
        self.project_tree.setColumnWidth(0, 300)
        self.project_tree.itemSelectionChanged.connect(self._show_selected_file)
        self.inventory_tabs.addTab(self.project_tree, "Project files")

        self.changes_tree = QTreeWidget()
        self.changes_tree.setHeaderLabels(["Path", "Change", "Category", "Size", "Purpose"])
        self.changes_tree.setColumnWidth(0, 300)
        self.changes_summary = QLabel("Run the project to compare files.")
        self.show_unchanged = QCheckBox("Show unchanged files")
        self.show_unchanged.setChecked(False)
        self.show_unchanged.toggled.connect(self._render_run_changes)
        changes_page = QWidget()
        changes_layout = QVBoxLayout(changes_page)
        changes_layout.setContentsMargins(10, 10, 10, 10)
        changes_layout.setSpacing(8)
        changes_layout.addWidget(self.changes_summary)
        changes_layout.addWidget(self.show_unchanged)
        changes_layout.addWidget(self.changes_tree, 1)
        self.inventory_tabs.addTab(changes_page, "Run changes")

        self.log_view = QPlainTextEdit()
        self.log_view.setReadOnly(True)
        log_page = QWidget()
        log_layout = QVBoxLayout(log_page)
        log_layout.setContentsMargins(10, 10, 10, 10)
        log_layout.addWidget(self.log_view)
        self.inventory_tabs.addTab(log_page, "Run log")

        self.file_details = QPlainTextEdit()
        self.file_details.setReadOnly(True)
        self.file_details.setMaximumHeight(85)
        self.file_details.setPlaceholderText("Select a file to see its classification and purpose.")
        self.inventory_tabs.setDocumentMode(True)
        self.inventory_tabs.setMovable(True)
        self.inventory_tabs.setMinimumHeight(280)
        self.inventory_tabs.setTabPosition(QTabWidget.TabPosition.North)
        self.inventory_tabs.setTabText(2, "Run information")
        root.addWidget(self.inventory_tabs, 1)
        self.project_status = QLabel("Ready")
        self.project_status.setWordWrap(True)
        root.addWidget(self.project_status)

        self._refresh_project_inventory()

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

    def _remember_view(self, key: str, window: QMainWindow) -> None:
        existing = self._view_windows.get(key)
        if existing is not None:
            existing.show()
            existing.raise_()
            existing.activateWindow()
            return
        self._view_windows[key] = window
        window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        window.destroyed.connect(lambda *_args, view_key=key: self._forget_view(view_key))
        window.show()
        window.raise_()
        window.activateWindow()

    def _forget_view(self, key: str) -> None:
        self._view_windows.pop(key, None)
        if key == "files":
            self._detached_project_tree = None
        elif key == "changes":
            self._detached_changes_tree = None
            self._detached_changes_summary = None
            self._detached_show_unchanged = None
        elif key == "log":
            self._detached_log = None

    def _activate_existing_view(self, key: str) -> bool:
        window = self._view_windows.get(key)
        if window is None:
            return False
        window.show()
        window.raise_()
        window.activateWindow()
        return True

    @staticmethod
    def _copy_tree_item(source: QTreeWidgetItem, target_parent) -> QTreeWidgetItem:
        item = QTreeWidgetItem([source.text(column) for column in range(source.columnCount())])
        item.setData(0, Qt.ItemDataRole.UserRole, source.data(0, Qt.ItemDataRole.UserRole))
        for column in range(source.columnCount()):
            item.setToolTip(column, source.toolTip(column))
        if target_parent is None:
            return item
        target_parent.addChild(item)
        return item

    def _clone_tree(self, source: QTreeWidget, target: QTreeWidget) -> None:
        def clone_children(source_parent, target_parent):
            for index in range(source_parent.childCount()):
                source_item = source_parent.child(index)
                target_item = self._copy_tree_item(source_item, target_parent)
                clone_children(source_item, target_item)

        for index in range(source.topLevelItemCount()):
            source_item = source.topLevelItem(index)
            target_item = self._copy_tree_item(source_item, None)
            target.addTopLevelItem(target_item)
            clone_children(source_item, target_item)
        target.setSortingEnabled(True)
        target.sortItems(0, Qt.SortOrder.AscendingOrder)
        target.expandToDepth(1)

    def _open_project_files(self) -> None:
        if self._activate_existing_view("files"):
            return
        window = QMainWindow(self)
        window.setWindowTitle(f"Project files - {self._project.name}")
        window.resize(1450, 900)
        central = QWidget(window)
        layout = QVBoxLayout(central)
        splitter = QSplitter(Qt.Orientation.Vertical, central)
        tree = QTreeWidget(splitter)
        tree.setHeaderLabels(["Path", "Kind", "Category", "Size", "Purpose"])
        tree.setColumnWidth(0, 420)
        tree.setColumnWidth(4, 420)
        self._clone_tree(self.project_tree, tree)
        self._detached_project_tree = tree
        details = QPlainTextEdit(splitter)
        details.setReadOnly(True)
        details.setMinimumHeight(150)
        details.setPlaceholderText("Select a file to inspect its classification and purpose.")
        tree.itemSelectionChanged.connect(
            lambda: self._show_entry_in_view(tree, details)
        )
        splitter.addWidget(tree)
        splitter.addWidget(details)
        splitter.setSizes([650, 220])
        layout.addWidget(splitter)
        window.setCentralWidget(central)
        self._remember_view("files", window)

    @staticmethod
    def _show_entry_in_view(tree: QTreeWidget, details: QPlainTextEdit) -> None:
        selected = tree.selectedItems()
        if not selected:
            return
        entry = selected[0].data(0, Qt.ItemDataRole.UserRole)
        if not isinstance(entry, ProjectFileEntry):
            return
        details.setPlainText(
            f"Path: {entry.relative_path}\nKind: {entry.kind}\n"
            f"Category: {entry.category}\nSize: {entry.size_bytes:,} bytes\n"
            f"Purpose: {entry.purpose}"
        )

    def _open_run_changes(self) -> None:
        if self._activate_existing_view("changes"):
            return
        window = QMainWindow(self)
        window.setWindowTitle(f"Run changes - {self._project.name}")
        window.resize(1450, 900)
        central = QWidget(window)
        layout = QVBoxLayout(central)
        summary = QLabel(self.changes_summary.text())
        summary.setWordWrap(True)
        unchanged = QCheckBox("Show unchanged files")
        unchanged.setChecked(self.show_unchanged.isChecked())
        tree = QTreeWidget()
        tree.setHeaderLabels(["Path", "Change", "Category", "Size", "Purpose"])
        tree.setColumnWidth(0, 420)
        tree.setColumnWidth(4, 420)
        layout.addWidget(summary)
        layout.addWidget(unchanged)
        layout.addWidget(tree, 1)
        self._detached_changes_tree = tree
        self._detached_changes_summary = summary
        self._detached_show_unchanged = unchanged
        unchanged.toggled.connect(self._render_run_changes)
        self._render_run_changes()
        window.setCentralWidget(central)
        self._remember_view("changes", window)

    def _open_run_log(self) -> None:
        if self._activate_existing_view("log"):
            return
        window = QMainWindow(self)
        window.setWindowTitle(f"Run log - {self._project.name}")
        window.resize(1200, 850)
        log = QPlainTextEdit()
        log.setReadOnly(True)
        log.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        log.setPlainText(self.log_view.toPlainText())
        window.setCentralWidget(log)
        self._detached_log = log
        self._remember_view("log", window)

    def _refresh_project_inventory(self) -> None:
        self.project_tree.clear()
        entries = scan_project(self._project.root)
        items: dict[str, QTreeWidgetItem] = {}
        ordered = sorted(
            entries.values(),
            key=lambda entry: (len(Path(entry.relative_path).parts), entry.relative_path.casefold()),
        )
        for entry in ordered:
            relative = Path(entry.relative_path)
            size = "" if entry.kind == "directory" else f"{entry.size_bytes:,} B"
            item = QTreeWidgetItem([
                relative.name, entry.kind, entry.category, size, entry.purpose
            ])
            item.setData(0, Qt.ItemDataRole.UserRole, entry)
            item.setToolTip(0, entry.relative_path)
            parent_path = relative.parent.as_posix()
            parent_item = items.get(parent_path)
            if parent_item is None or parent_path == ".":
                self.project_tree.addTopLevelItem(item)
            else:
                parent_item.addChild(item)
            items[entry.relative_path] = item
        self.project_tree.sortItems(0, Qt.SortOrder.AscendingOrder)
        self.project_tree.expandToDepth(1)
        if self._detached_project_tree is not None:
            self._detached_project_tree.clear()
            self._clone_tree(self.project_tree, self._detached_project_tree)

    @staticmethod
    def _add_inventory_item(tree: QTreeWidget, entry: ProjectFileEntry, status: str) -> None:
        size = "" if entry.kind == "directory" else f"{entry.size_bytes:,} B"
        item = QTreeWidgetItem([
            entry.relative_path, entry.kind, entry.category, size, entry.purpose
        ])
        item.setData(0, Qt.ItemDataRole.UserRole, entry)
        item.setToolTip(0, entry.relative_path)
        tree.addTopLevelItem(item)

    def _show_selected_file(self) -> None:
        selected = self.project_tree.selectedItems()
        if not selected:
            return
        entry = selected[0].data(0, Qt.ItemDataRole.UserRole)
        if not isinstance(entry, ProjectFileEntry):
            return
        self.file_details.setPlainText(
            f"Path: {entry.relative_path}\n"
            f"Kind: {entry.kind}\n"
            f"Category: {entry.category}\n"
            f"Size: {entry.size_bytes:,} bytes\n"
            f"Purpose: {entry.purpose}"
        )

    def _refresh_run_changes(self, execution_directory: Path, run_id: str) -> None:
        self.changes_tree.clear()
        current = scan_project(execution_directory)
        changes = list(compare_inventory(self._inventory_baseline, current))
        run_root = execution_directory / ".oneflow-runs" / run_id
        run_outputs = scan_workflow_run_outputs(execution_directory, run_root)
        counts = {"unchanged": 0, "added": 0, "modified": 0, "deleted": 0}
        changes = [
            change for change in changes
            if change.entry.kind != "directory"
        ]
        for change in changes:
            if change.status in counts:
                counts[change.status] += 1
        changes.extend(
            FileChange(entry, "generated in run")
            for entry in run_outputs
        )
        self._run_changes = changes
        self._run_output_count = len(run_outputs)
        self._run_change_counts = counts
        self.changes_summary.setText(
            f"Added: {counts['added']}    Modified: {counts['modified']}    "
            f"Deleted: {counts['deleted']}    Unchanged: {counts['unchanged']}    "
            f"Workflow run outputs: {len(run_outputs)}"
        )
        self._render_run_changes()
        self.inventory_tabs.setCurrentIndex(1)

    def _render_run_changes(self) -> None:
        self.changes_tree.clear()
        for change in self._run_changes:
            if change.status == "unchanged" and not self.show_unchanged.isChecked():
                continue
            entry = change.entry
            size = "" if entry.kind == "directory" else f"{entry.size_bytes:,} B"
            item = QTreeWidgetItem([
                entry.relative_path, change.status, entry.category, size, entry.purpose
            ])
            item.setData(0, Qt.ItemDataRole.UserRole, entry)
            item.setToolTip(0, entry.relative_path)
            self.changes_tree.addTopLevelItem(item)
        self.changes_tree.sortItems(0, Qt.SortOrder.AscendingOrder)
        if self._detached_changes_tree is not None:
            self._detached_changes_tree.clear()
            for change in self._run_changes:
                show_unchanged = (
                    self._detached_show_unchanged.isChecked()
                    if self._detached_show_unchanged is not None
                    else self.show_unchanged.isChecked()
                )
                if change.status == "unchanged" and not show_unchanged:
                    continue
                entry = change.entry
                size = "" if entry.kind == "directory" else f"{entry.size_bytes:,} B"
                item = QTreeWidgetItem([
                    entry.relative_path, change.status, entry.category, size, entry.purpose
                ])
                item.setData(0, Qt.ItemDataRole.UserRole, entry)
                item.setToolTip(0, entry.relative_path)
                self._detached_changes_tree.addTopLevelItem(item)
            self._detached_changes_tree.sortItems(0, Qt.SortOrder.AscendingOrder)
        if self._detached_changes_summary is not None:
            self._detached_changes_summary.setText(self.changes_summary.text())

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
        self._inventory_baseline = scan_project(self._project.root)
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
        self.project_status.setText(f"Running: {self._project.name}  |  Run ID: {run_id}")
        self.inventory_tabs.setCurrentWidget(self.log_view)
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
        if self._detached_log is not None:
            self._detached_log.appendPlainText(line)
        if self._log_path is not None:
            try:
                with self._log_path.open("a", encoding="utf-8") as stream:
                    stream.write(line + "\n")
            except OSError as error:
                self.statusBar().showMessage(f"Could not write run log: {error}")

    def _run_completed(self, outcome: ProjectRunOutcome) -> None:
        self._append_log(f"{outcome.message}\nExecution directory: {outcome.project_directory}")
        self._refresh_project_inventory()
        self._refresh_run_changes(outcome.project_directory, outcome.run_id)
        self.inventory_tabs.setCurrentWidget(self.changes_tree)
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
        self.project_status.setText(
            f"Last run: {outcome.status}. Execution directory: {outcome.project_directory}"
        )
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
