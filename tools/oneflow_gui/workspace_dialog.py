"""Create and edit OneFLOW workspace definitions."""

from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
)

from .user_data_paths import UserDataPaths
from .workspace_workflow_adapter import store_workspace_case_path
from .workspace_repository import (
    WorkspaceArtifactInput,
    WorkspaceArtifactOutput,
    WorkspaceCase,
    WorkspaceDefinition,
    load_workspace,
    save_workspace,
)


class WorkspaceDialog(QDialog):
    """Edit case membership and execution policy for a .ofw workspace."""

    def __init__(self, parent=None, workspace_path: Path | None = None) -> None:
        super().__init__(parent)
        self._workspace_path = workspace_path
        self.run_requested = False
        self._run_after_save = False
        self.setWindowTitle("OneFLOW Multi-Case Workspace (.ofw)")
        self.resize(820, 540)

        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("This .ofw file groups existing case directories and records their execution policy. It is separate from a project folder created or opened through New Project / Open Project."))
        layout.addWidget(QLabel("In dependency mode, declare artifact inputs and outputs as JSON on each case row."))

        name_row = QHBoxLayout()
        name_row.addWidget(QLabel("Workspace name"))
        self.name_edit = QLineEdit()
        self.name_edit.setPlaceholderText("My OneFLOW workspace")
        name_row.addWidget(self.name_edit, 1)
        layout.addLayout(name_row)

        mode_row = QHBoxLayout()
        mode_row.addWidget(QLabel("Execution mode"))
        self.mode_combo = QComboBox()
        self.mode_combo.addItem("Single case", "single_case")
        self.mode_combo.addItem("Sequential cases", "sequential_cases")
        self.mode_combo.addItem("Independent cases (one process)", "independent_cases")
        self.mode_combo.addItem("Dependency workflow", "dependency_workflow")
        mode_row.addWidget(self.mode_combo, 1)
        layout.addLayout(mode_row)

        self.case_table = QTableWidget(0, 6)
        self.case_table.setHorizontalHeaderLabels(["Case ID", "Case directory (editable)", "Depends on (comma-separated IDs)", "Inputs (JSON)", "Outputs (JSON)", "Browse"])
        self.case_table.horizontalHeader().setStretchLastSection(True)
        layout.addWidget(self.case_table, 1)

        case_buttons = QHBoxLayout()
        add_case = QPushButton("Add case directory...")
        add_case.clicked.connect(self._add_case)
        remove_case = QPushButton("Remove selected")
        remove_case.clicked.connect(self._remove_case)
        case_buttons.addWidget(add_case)
        case_buttons.addWidget(remove_case)
        case_buttons.addStretch(1)
        layout.addLayout(case_buttons)

        footer = QHBoxLayout()
        footer.addStretch(1)
        save_button = QPushButton("Save workspace...")
        save_button.clicked.connect(self._save)
        run_button = QPushButton("Save and Run...")
        run_button.clicked.connect(self._save_and_run)
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        footer.addWidget(save_button)
        footer.addWidget(run_button)
        footer.addWidget(close_button)
        layout.addLayout(footer)

        if workspace_path is not None:
            self._load(workspace_path)

    def _load(self, path: Path) -> None:
        try:
            workspace = load_workspace(path)
        except (OSError, ValueError) as error:
            QMessageBox.critical(self, "Cannot open workspace", str(error))
            return
        self._workspace_path = path.resolve()
        self.name_edit.setText(workspace.name)
        index = self.mode_combo.findData(workspace.execution_mode)
        if index >= 0:
            self.mode_combo.setCurrentIndex(index)
        for case in workspace.cases:
            self._append_case(
                case.case_id,
                case.path,
                ", ".join(case.depends_on),
                json.dumps([item.to_dict() for item in case.inputs]),
                json.dumps([item.to_dict() for item in case.outputs]),
            )

    def _append_case(
        self, case_id: str, path: str, dependencies: str = "",
        inputs: str = "[]", outputs: str = "[]",
    ) -> None:
        row = self.case_table.rowCount()
        self.case_table.insertRow(row)
        self.case_table.setItem(row, 0, QTableWidgetItem(case_id))
        self.case_table.setItem(row, 1, QTableWidgetItem(path))
        self.case_table.setItem(row, 2, QTableWidgetItem(dependencies))
        self.case_table.setItem(row, 3, QTableWidgetItem(inputs))
        self.case_table.setItem(row, 4, QTableWidgetItem(outputs))
        self._set_browse_button(row)

    def _set_browse_button(self, row: int) -> None:
        """Attach a directory picker to the row without making paths read-only."""
        button = QPushButton("Browse...")
        button.clicked.connect(lambda _checked=False, index=row: self._browse_case_directory(index))
        self.case_table.setCellWidget(row, 5, button)

    def _browse_case_directory(self, row: int) -> None:
        """Replace one case path using a native directory picker."""
        path_item = self.case_table.item(row, 1)
        current_path = path_item.text().strip() if path_item is not None else ""
        start = Path(current_path).expanduser() if current_path else UserDataPaths.discover().projects
        if not start.is_dir():
            start = Path.home()
        if start.is_file():
            start = start.parent
        if not start.is_dir():
            start = Path.home()
        directory = QFileDialog.getExistingDirectory(
            self, "Select case directory", str(start)
        )
        if directory:
            if path_item is None:
                path_item = QTableWidgetItem()
                self.case_table.setItem(row, 1, path_item)
            path_item.setText(directory)

    def _add_case(self) -> None:
        projects_directory = UserDataPaths.discover().projects
        start_directory = projects_directory if projects_directory.is_dir() else Path.home()
        directory = QFileDialog.getExistingDirectory(
            self, "Add case directory", str(start_directory)
        )
        if not directory:
            return
        root = Path(directory).resolve()
        case_id = "".join(
            character.lower() if character.isascii() and character.isalnum() else "-"
            for character in root.name
        ).strip("-") or "case"
        existing = {
            self.case_table.item(row, 0).text()
            for row in range(self.case_table.rowCount())
            if self.case_table.item(row, 0) is not None
        }
        base_id = case_id or "case"
        suffix = 2
        while case_id in existing:
            case_id = f"{base_id}-{suffix}"
            suffix += 1
        self._append_case(case_id, str(root))

    def _remove_case(self) -> None:
        rows = sorted({index.row() for index in self.case_table.selectedIndexes()}, reverse=True)
        for row in rows:
            self.case_table.removeRow(row)
        # Row indexes change after removal, so reconnect the remaining buttons.
        for row in range(self.case_table.rowCount()):
            self._set_browse_button(row)

    @property
    def workspace_path(self) -> Path | None:
        return self._workspace_path

    def _save_and_run(self) -> None:
        self._run_after_save = True
        self._save()

    def _save(self) -> None:
        run_after_save = self._run_after_save
        self._run_after_save = False
        name = self.name_edit.text().strip()
        if not name:
            QMessageBox.warning(self, "Workspace name required", "Enter a workspace name.")
            return
        cases = []
        for row in range(self.case_table.rowCount()):
            case_id_item = self.case_table.item(row, 0)
            path_item = self.case_table.item(row, 1)
            dependency_item = self.case_table.item(row, 2)
            if case_id_item is None or path_item is None:
                continue
            case_id = case_id_item.text().strip()
            case_path = path_item.text().strip()
            if not case_id or not case_path:
                QMessageBox.warning(self, "Incomplete case", f"Complete the fields in row {row + 1}.")
                return
            dependency_text = dependency_item.text().strip() if dependency_item else ""
            dependencies = tuple(
                dependency.strip() for dependency in dependency_text.split(",")
                if dependency.strip()
            )
            try:
                raw_inputs = json.loads(self.case_table.item(row, 3).text() or "[]")
                raw_outputs = json.loads(self.case_table.item(row, 4).text() or "[]")
                if not isinstance(raw_inputs, list) or not isinstance(raw_outputs, list):
                    raise ValueError("Inputs and outputs must each be JSON arrays.")
                if any(not isinstance(item, dict) for item in raw_inputs + raw_outputs):
                    raise ValueError("Each input and output binding must be a JSON object.")
            except (json.JSONDecodeError, ValueError) as error:
                QMessageBox.warning(
                    self, "Invalid artifact bindings",
                    f"Check the inputs/outputs JSON in row {row + 1}: {error}",
                )
                return
            cases.append(
                WorkspaceCase(
                    case_id, case_path, dependencies,
                    tuple(
                        WorkspaceArtifactInput(
                            item.get("artifact"), item.get("setting"),
                            item.get("parameter_file"),
                        )
                        for item in raw_inputs
                    ),
                    tuple(
                        WorkspaceArtifactOutput(
                            item.get("id"), item.get("path"),
                            item.get("parameter_file"), item.get("setting"),
                        )
                        for item in raw_outputs
                    ),
                )
            )

        if not cases:
            QMessageBox.warning(self, "No cases", "Add at least one case directory.")
            return

        if self._workspace_path is None:
            suggested = (name.lower().replace(" ", "-") or "workspace") + ".ofw"
            workspace_directory = UserDataPaths.discover().workspaces
            workspace_directory.mkdir(parents=True, exist_ok=True)
            filename, _ = QFileDialog.getSaveFileName(
                self, "Save OneFLOW workspace", str(workspace_directory / suggested),
                "OneFLOW workspace (*.ofw)",
            )
            if not filename:
                return
            self._workspace_path = Path(filename)
            if self._workspace_path.suffix.lower() != ".ofw":
                self._workspace_path = self._workspace_path.with_suffix(".ofw")

        cases = [
            WorkspaceCase(
                case.case_id,
                store_workspace_case_path(case.path, self._workspace_path),
                case.depends_on,
                case.inputs,
                case.outputs,
            )
            for case in cases
        ]

        workspace_id = "".join(
            character.lower() if character.isascii() and character.isalnum() else "-"
            for character in self._workspace_path.stem
        ).strip("-._") or "workspace"
        try:
            workspace = WorkspaceDefinition(
                workspace_id=workspace_id,
                name=name,
                execution_mode=self.mode_combo.currentData(),
                cases=tuple(cases),
            )
            save_workspace(self._workspace_path, workspace)
        except (OSError, ValueError) as error:
            QMessageBox.critical(self, "Cannot save workspace", str(error))
            return
        QMessageBox.information(
            self, "Workspace saved", f"Saved workspace to:\n{self._workspace_path}"
        )
        self.run_requested = run_after_save
        self.accept()
