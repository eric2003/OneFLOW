"""Dialog for reviewing and controlling registered background tasks."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QTimer
from PySide6.QtWidgets import (
    QDialog,
    QHBoxLayout,
    QMessageBox,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
)

from .background_task_manager import request_background_stop
from .background_task_registry import ACTIVE_STATUSES, BackgroundTask, BackgroundTaskRegistry


class BackgroundTaskDialog(QDialog):
    """List persisted tasks and provide explicit stop and unregister actions."""

    def __init__(self, root: Path, parent=None) -> None:
        super().__init__(parent)
        self._registry = BackgroundTaskRegistry(root)
        self._tasks: list[BackgroundTask] = []
        self.setWindowTitle("Background Tasks")
        self.resize(1100, 520)

        layout = QVBoxLayout(self)
        self.table = QTableWidget(0, 6, self)
        self.table.setHorizontalHeaderLabels(
            ["Task ID", "Status", "Project", "Started", "Host PID", "Message"]
        )
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QTableWidget.SelectionMode.SingleSelection)
        self.table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.horizontalHeader().setSectionResizeMode(2, self.table.horizontalHeader().ResizeMode.Stretch)
        layout.addWidget(self.table, 1)

        actions = QHBoxLayout()
        self.refresh_button = QPushButton("Refresh")
        self.refresh_button.clicked.connect(self.refresh)
        actions.addWidget(self.refresh_button)
        self.stop_button = QPushButton("Stop Task")
        self.stop_button.clicked.connect(self._stop_selected)
        actions.addWidget(self.stop_button)
        self.unregister_button = QPushButton("Unregister Finished Task")
        self.unregister_button.clicked.connect(self._unregister_selected)
        actions.addWidget(self.unregister_button)
        actions.addStretch(1)
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.accept)
        actions.addWidget(close_button)
        layout.addLayout(actions)

        self._timer = QTimer(self)
        self._timer.setInterval(1000)
        self._timer.timeout.connect(self.refresh)
        self._timer.start()
        self.refresh()

    def refresh(self) -> None:
        selected_id = self._selected_task_id()
        self._tasks = self._registry.list_tasks()
        self.table.setRowCount(len(self._tasks))
        selected_row = -1
        for row, task in enumerate(self._tasks):
            values = (
                task.task_id,
                task.status,
                task.project_directory,
                task.started_at,
                str(task.host_pid or ""),
                task.message,
            )
            for column, value in enumerate(values):
                self.table.setItem(row, column, QTableWidgetItem(value))
            if task.task_id == selected_id:
                selected_row = row
        if selected_row >= 0:
            self.table.selectRow(selected_row)
        self._update_actions()

    def _selected_task_id(self) -> str | None:
        row = self.table.currentRow()
        if row < 0 or row >= len(self._tasks):
            return None
        return self._tasks[row].task_id

    def _selected_task(self) -> BackgroundTask | None:
        task_id = self._selected_task_id()
        return next((task for task in self._tasks if task.task_id == task_id), None)

    def _update_actions(self) -> None:
        task = self._selected_task()
        self.stop_button.setEnabled(task is not None and task.status in ACTIVE_STATUSES)
        self.unregister_button.setEnabled(task is not None and task.status not in ACTIVE_STATUSES)

    def _stop_selected(self) -> None:
        task = self._selected_task()
        if task is None or task.status not in ACTIVE_STATUSES:
            return
        answer = QMessageBox.question(
            self,
            "Stop background task",
            f"Request stop for task {task.task_id}? If it is a workflow, subsequent cases should not start.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if answer != QMessageBox.StandardButton.Yes:
            return
        try:
            request_background_stop(self._registry.root, task.task_id)
        except (OSError, ValueError) as error:
            QMessageBox.critical(self, "Cannot stop task", str(error))
        self.refresh()

    def _unregister_selected(self) -> None:
        task = self._selected_task()
        if task is None:
            return
        if task.status in ACTIVE_STATUSES:
            QMessageBox.warning(
                self,
                "Task is active",
                "Stop the task and wait for a terminal state before unregistering it.",
            )
            return
        answer = QMessageBox.question(
            self,
            "Unregister finished task",
            f"Remove the history record for {task.task_id}? The log file will be kept.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if answer != QMessageBox.StandardButton.Yes:
            return
        self._registry.remove(task.task_id)
        self.refresh()
