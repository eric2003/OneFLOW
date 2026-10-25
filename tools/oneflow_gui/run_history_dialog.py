"""Browse and manage persisted workflow run history without deleting artifacts."""

from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
)

from .workflow_repository import WorkflowRepository, WorkflowRunRecord


class RunHistoryDialog(QDialog):
    """Inspect run details/logs, edit notes, and remove selected history records."""

    def __init__(self, repository_root: Path, parent=None) -> None:
        super().__init__(parent)
        self._repository = WorkflowRepository(repository_root)
        self._records: list[WorkflowRunRecord] = []
        self.setWindowTitle("OneFLOW Run History")
        self.resize(1000, 720)

        root = QVBoxLayout(self)
        root.addWidget(QLabel("Filter history records; cleanup affects records only, not run artifacts."))
        filters = QFormLayout()
        self.workflow_filter = QLineEdit()
        self.workflow_filter.setPlaceholderText("Any workflow ID")
        self.status_filter = QComboBox()
        self.status_filter.addItems(["all", "running", "succeeded", "failed", "cancelled", "pending"])
        self.before_filter = QLineEdit()
        self.before_filter.setPlaceholderText("Optional ISO timestamp, e.g. 2026-10-15T00:00:00Z")
        filters.addRow("Workflow ID", self.workflow_filter)
        filters.addRow("Status", self.status_filter)
        filters.addRow("Started before", self.before_filter)
        root.addLayout(filters)

        filter_actions = QHBoxLayout()
        refresh = QPushButton("Refresh")
        refresh.clicked.connect(self._reload)
        self.cleanup_button = QPushButton("Delete matching history")
        self.cleanup_button.clicked.connect(self._cleanup)
        filter_actions.addWidget(refresh)
        filter_actions.addWidget(self.cleanup_button)
        filter_actions.addStretch(1)
        root.addLayout(filter_actions)

        self.table = QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels(["Started", "Status", "Task Type", "Workflow", "Run ID", "Note"])
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setColumnWidth(0, 190)
        self.table.setColumnWidth(1, 100)
        self.table.setColumnWidth(2, 150)
        self.table.setColumnWidth(3, 170)
        self.table.setColumnWidth(4, 210)
        self.table.itemSelectionChanged.connect(self._show_selected)
        root.addWidget(self.table, 2)

        root.addWidget(QLabel("Run details / log"))
        self.details = QPlainTextEdit()
        self.details.setReadOnly(True)
        root.addWidget(self.details, 3)

        note_actions = QHBoxLayout()
        note_actions.addWidget(QLabel("Note"))
        self.note_edit = QLineEdit()
        note_actions.addWidget(self.note_edit, 1)
        self.save_note_button = QPushButton("Save note")
        self.save_note_button.clicked.connect(self._save_note)
        self.delete_button = QPushButton("Delete selected record")
        self.delete_button.clicked.connect(self._delete_selected)
        note_actions.addWidget(self.save_note_button)
        note_actions.addWidget(self.delete_button)
        root.addLayout(note_actions)
        self._reload()

    def _selected_record(self) -> WorkflowRunRecord | None:
        row = self.table.currentRow()
        if row < 0 or row >= len(self._records):
            return None
        return self._records[row]

    def _reload(self) -> None:
        workflow_id = self.workflow_filter.text().strip() or None
        status = self.status_filter.currentText()
        started_before = self.before_filter.text().strip() or None
        try:
            self._records = self._repository.list_runs(
                workflow_id,
                status=None if status == "all" else status,
                started_before=started_before,
            )
        except (OSError, ValueError) as error:
            QMessageBox.critical(self, "Cannot load run history", str(error))
            return
        self.table.setRowCount(0)
        for record in self._records:
            row = self.table.rowCount()
            self.table.insertRow(row)
            for column, value in enumerate((
                record.started_at, record.status, record.task_type,
                record.workflow_id, record.run_id, record.note
            )):
                self.table.setItem(row, column, QTableWidgetItem(value))
        if self._records:
            self.table.selectRow(0)
        else:
            self.details.clear()
            self.note_edit.clear()

    def _show_selected(self) -> None:
        record = self._selected_record()
        if record is None:
            return
        details = ["Run record:", json.dumps(record.to_dict(), indent=2, ensure_ascii=True)]
        if record.log_path:
            log_path = Path(record.log_path)
            details.append("")
            details.append(f"Log path: {log_path}")
            try:
                details.append(log_path.read_text(encoding="utf-8", errors="replace"))
            except OSError as error:
                details.append(f"Log unavailable: {error}")
        self.details.setPlainText("\n".join(details))
        self.note_edit.setText(record.note)

    def _save_note(self) -> None:
        record = self._selected_record()
        if record is None:
            QMessageBox.information(self, "No run selected", "Select a history record first.")
            return
        try:
            self._repository.update_run_note(record.run_id, self.note_edit.text())
        except (OSError, ValueError) as error:
            QMessageBox.critical(self, "Cannot save note", str(error))
            return
        self._reload()

    def _delete_selected(self) -> None:
        record = self._selected_record()
        if record is None:
            return
        answer = QMessageBox.question(
            self,
            "Delete history record",
            f"Delete history record '{record.run_id}'? The log and run artifacts will be kept.",
        )
        if answer != QMessageBox.StandardButton.Yes:
            return
        try:
            self._repository.delete_run(record.run_id)
        except (OSError, ValueError) as error:
            QMessageBox.critical(self, "Cannot delete history", str(error))
            return
        self._reload()

    def _cleanup(self) -> None:
        workflow_id = self.workflow_filter.text().strip() or None
        status = self.status_filter.currentText()
        statuses = None if status == "all" else {status}
        started_before = self.before_filter.text().strip() or None
        description = []
        if workflow_id:
            description.append(f"workflow={workflow_id}")
        if statuses:
            description.append(f"status={status}")
        if started_before:
            description.append(f"started before {started_before}")
        scope = ", ".join(description) if description else "all history records"
        answer = QMessageBox.question(
            self,
            "Delete matching history",
            f"Delete {scope}? Run logs and artifacts will not be deleted.",
        )
        if answer != QMessageBox.StandardButton.Yes:
            return
        try:
            count = self._repository.clear_runs(
                workflow_id=workflow_id,
                statuses=statuses,
                started_before=started_before,
            )
        except (OSError, ValueError) as error:
            QMessageBox.critical(self, "Cannot clean history", str(error))
            return
        QMessageBox.information(self, "History cleaned", f"Deleted {count} history record(s).")
        self._reload()
