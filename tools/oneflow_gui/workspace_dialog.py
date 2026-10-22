"""Create and edit OneFLOW workspace definitions."""

from __future__ import annotations

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

from .workspace_repository import (
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
        self.setWindowTitle("OneFLOW Workspace")
        self.resize(820, 540)

        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("A workspace records case paths and how its cases should execute."))

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
        self.mode_combo.addItem("Independent cases (parallel-ready)", "independent_cases")
        self.mode_combo.addItem("Dependency workflow", "dependency_workflow")
        mode_row.addWidget(self.mode_combo, 1)
        layout.addLayout(mode_row)

        self.case_table = QTableWidget(0, 3)
        self.case_table.setHorizontalHeaderLabels(["Case ID", "Case directory", "Depends on (comma-separated IDs)"])
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
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.reject)
        footer.addWidget(save_button)
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
            self._append_case(case.case_id, case.path, ", ".join(case.depends_on))

    def _append_case(self, case_id: str, path: str, dependencies: str = "") -> None:
        row = self.case_table.rowCount()
        self.case_table.insertRow(row)
        self.case_table.setItem(row, 0, QTableWidgetItem(case_id))
        self.case_table.setItem(row, 1, QTableWidgetItem(path))
        self.case_table.setItem(row, 2, QTableWidgetItem(dependencies))

    def _add_case(self) -> None:
        directory = QFileDialog.getExistingDirectory(
            self, "Add case directory", str(Path.home())
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

    def _save(self) -> None:
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
            if self._workspace_path is not None:
                try:
                    case_path = str(Path(case_path).expanduser().resolve().relative_to(self._workspace_path.parent.resolve()))
                except ValueError:
                    case_path = str(Path(case_path).expanduser().resolve())
            cases.append(WorkspaceCase(case_id, case_path, dependencies))

        if not cases:
            QMessageBox.warning(self, "No cases", "Add at least one case directory.")
            return

        if self._workspace_path is None:
            suggested = (name.lower().replace(" ", "-") or "workspace") + ".ofw"
            filename, _ = QFileDialog.getSaveFileName(
                self, "Save OneFLOW workspace", str(Path.home() / suggested),
                "OneFLOW workspace (*.ofw)",
            )
            if not filename:
                return
            self._workspace_path = Path(filename)
            if self._workspace_path.suffix.lower() != ".ofw":
                self._workspace_path = self._workspace_path.with_suffix(".ofw")

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
        self.accept()
