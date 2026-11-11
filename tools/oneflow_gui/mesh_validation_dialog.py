"""Dedicated Qt workbench for end-to-end Mesh-to-CFD validation."""

from __future__ import annotations

import subprocess
import sys
import threading
from pathlib import Path

from PySide6.QtCore import QThread, Signal, Qt
from PySide6.QtWidgets import (
    QComboBox, QDialog, QFormLayout, QHBoxLayout, QHeaderView, QLabel,
    QMessageBox, QPlainTextEdit, QPushButton, QSpinBox, QTreeWidget,
    QTreeWidgetItem, QVBoxLayout,
)

from .executable_discovery import executable_candidates
from .mesh_validation import cfd_template_choices, mesh_case_groups
from .process_control import process_group_options, terminate_process_tree
from .user_data_paths import UserDataPaths


class MeshValidationWorker(QThread):
    """Run all selected validation Cases, continuing after individual failures."""

    output_line = Signal(str)
    completed = Signal(bool, str)

    def __init__(self, commands: list[tuple[str, tuple[str, ...], Path]]) -> None:
        super().__init__()
        self._commands = commands
        self._cancel = threading.Event()
        self._process: subprocess.Popen[str] | None = None
        self._lock = threading.Lock()

    def request_stop(self) -> None:
        self._cancel.set()
        with self._lock:
            process = self._process
        if process is not None:
            terminate_process_tree(process)

    def run(self) -> None:
        passed: list[str] = []
        failed: list[str] = []
        for label, arguments, working_directory in self._commands:
            if self._cancel.is_set():
                break
            self.output_line.emit(f"\n=== {label} ===")
            self.output_line.emit("$ " + subprocess.list2cmdline(list(arguments)))
            try:
                with subprocess.Popen(
                    list(arguments), cwd=working_directory,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    text=True, encoding="utf-8", errors="replace", bufsize=1,
                    **process_group_options(),
                ) as process:
                    with self._lock:
                        self._process = process
                    try:
                        if process.stdout is not None:
                            for line in process.stdout:
                                self.output_line.emit(line.rstrip("\r\n"))
                                if self._cancel.is_set():
                                    terminate_process_tree(process)
                                    break
                        return_code = process.wait()
                    finally:
                        with self._lock:
                            if self._process is process:
                                self._process = None
            except OSError as error:
                self.output_line.emit(f"ERROR: {error}")
                failed.append(label)
                continue
            if self._cancel.is_set():
                break
            if return_code == 0:
                passed.append(label)
                self.output_line.emit(f"PASSED: {label}")
            else:
                failed.append(label)
                self.output_line.emit(f"FAILED: {label} (exit {return_code})")
        cancelled = self._cancel.is_set()
        summary = [
            "Mesh-to-CFD Validation summary:",
            f"  Selected: {len(self._commands)}",
            f"  Passed: {len(passed)}",
            f"  Failed: {len(failed)}",
            f"  Not run: {len(self._commands) - len(passed) - len(failed)}",
        ]
        if failed:
            summary.append("  Failed Cases: " + ", ".join(failed))
        if cancelled:
            summary.append("  Run cancelled.")
        self.completed.emit(not failed and not cancelled, "\n".join(summary))


class MeshValidationDialog(QDialog):
    """Select multiple mesh producer Cases and run each complete validation chain."""

    def __init__(self, repository_root: Path, *, executable: str = "", parent=None) -> None:
        super().__init__(parent)
        self._repository_root = Path(repository_root).resolve()
        self._worker: MeshValidationWorker | None = None
        self.setWindowTitle("Mesh-to-CFD Validation")
        self.resize(920, 700)

        layout = QVBoxLayout(self)
        layout.addWidget(QLabel(
            "Select one or more native C++ or Python STL-to-Gmsh Cases. Each selected "
            "Case runs its applicable "
            "STL → Gmsh → CGNS → OneFLOW validation chain and a short solver smoke test. For 3D meshes, the shared consumer configuration is dimension-adjusted; this is not a numerical-accuracy claim. Re-running a Case "
            "clears only that Case's previous validation directory."
        ))

        self.case_tree = QTreeWidget()
        self.case_tree.setColumnCount(2)
        self.case_tree.setHeaderLabels(["Run", "STL / Mesh Case"])
        self.case_tree.setRootIsDecorated(True)
        self.case_tree.setAlternatingRowColors(True)
        self.case_tree.header().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        self._case_items: dict[str, QTreeWidgetItem] = {}
        try:
            dimension_groups = mesh_case_groups(self._repository_root)
            for dimension_label, choices in dimension_groups:
                dimension_item = QTreeWidgetItem(self.case_tree)
                dimension_item.setText(1, dimension_label)
                dimension_item.setFlags(
                    dimension_item.flags() | Qt.ItemFlag.ItemIsUserCheckable
                )
                category_item = QTreeWidgetItem(dimension_item)
                category_item.setText(1, "STL to Gmsh")
                category_item.setFlags(
                    category_item.flags() | Qt.ItemFlag.ItemIsUserCheckable
                )
                for case_id, title in choices:
                    item = QTreeWidgetItem(category_item)
                    item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                    item.setText(1, title)
                    item.setData(0, Qt.ItemDataRole.UserRole, case_id)
                    item.setCheckState(0, Qt.CheckState.Checked)
                    self._case_items[case_id] = item
                category_item.setCheckState(0, Qt.CheckState.Checked)
                dimension_item.setCheckState(0, Qt.CheckState.Checked)
            self.case_tree.itemChanged.connect(self._update_group_state)
            self.case_tree.expandToDepth(1)
        except (OSError, ValueError, KeyError) as error:
            self.case_tree.setEnabled(False)
            layout.addWidget(QLabel(f"Mesh Case catalog unavailable: {error}"))
        layout.addWidget(self.case_tree, 1)

        selection_actions = QHBoxLayout()
        select_all = QPushButton("Select all")
        select_all.clicked.connect(lambda: self._set_all(True))
        select_none = QPushButton("Select none")
        select_none.clicked.connect(lambda: self._set_all(False))
        selection_actions.addWidget(select_all)
        selection_actions.addWidget(select_none)
        selection_actions.addStretch(1)
        layout.addLayout(selection_actions)

        form = QFormLayout()
        self.executable_combo = QComboBox()
        self.executable_combo.setEditable(True)
        self.executable_combo.setMinimumContentsLength(48)
        self.executable_combo.addItems(
            [str(path) for path in executable_candidates(self._repository_root)]
        )
        if executable:
            self.executable_combo.setEditText(executable)
        elif self.executable_combo.count():
            self.executable_combo.setCurrentIndex(0)
        self.template_combo = QComboBox()
        for template_id, title in cfd_template_choices():
            self.template_combo.addItem(title, template_id)
        self.steps_spin = QSpinBox()
        self.steps_spin.setRange(1, 100000)
        self.steps_spin.setValue(2)
        self.steps_spin.setToolTip("Short CFD solver smoke-test iteration count.")
        form.addRow("OneFLOW executable", self.executable_combo)
        form.addRow("CFD consumer template (3D dimension is overridden from mesh metadata)", self.template_combo)
        form.addRow("CFD smoke-test steps", self.steps_spin)
        layout.addLayout(form)

        workspace = UserDataPaths.discover().mesh_validation
        self.workspace_label = QLabel("Workspace: " + str(workspace))
        self.workspace_label.setWordWrap(True)
        layout.addWidget(self.workspace_label)

        self.log_view = QPlainTextEdit()
        self.log_view.setReadOnly(True)
        self.log_view.setMaximumBlockCount(20000)
        self.log_view.setPlaceholderText("Validation output will appear here.")
        layout.addWidget(self.log_view, 1)

        buttons = QHBoxLayout()
        buttons.addStretch(1)
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self._start)
        self.stop_button = QPushButton("Stop")
        self.stop_button.setEnabled(False)
        self.stop_button.clicked.connect(self._stop)
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.close)
        buttons.addWidget(self.run_button)
        buttons.addWidget(self.stop_button)
        buttons.addWidget(close_button)
        layout.addLayout(buttons)

    def _set_subtree_state(self, item: QTreeWidgetItem, state: Qt.CheckState) -> None:
        item.setCheckState(0, state)
        for index in range(item.childCount()):
            self._set_subtree_state(item.child(index), state)

    def _refresh_ancestor_states(self, item: QTreeWidgetItem | None) -> None:
        while item is not None:
            states = [
                item.child(index).checkState(0)
                for index in range(item.childCount())
            ]
            if states and all(state == Qt.CheckState.Checked for state in states):
                item.setCheckState(0, Qt.CheckState.Checked)
            elif states and all(state == Qt.CheckState.Unchecked for state in states):
                item.setCheckState(0, Qt.CheckState.Unchecked)
            else:
                item.setCheckState(0, Qt.CheckState.PartiallyChecked)
            item = item.parent()

    def _set_all(self, checked: bool) -> None:
        state = Qt.CheckState.Checked if checked else Qt.CheckState.Unchecked
        self.case_tree.blockSignals(True)
        for index in range(self.case_tree.topLevelItemCount()):
            self._set_subtree_state(self.case_tree.topLevelItem(index), state)
        self.case_tree.blockSignals(False)

    def _update_group_state(self, item: QTreeWidgetItem, column: int) -> None:
        if column != 0:
            return
        self.case_tree.blockSignals(True)
        if item.childCount():
            state = item.checkState(0)
            if state == Qt.CheckState.PartiallyChecked:
                state = Qt.CheckState.Unchecked
            for index in range(item.childCount()):
                self._set_subtree_state(item.child(index), state)
            item.setCheckState(0, state)
        self._refresh_ancestor_states(item.parent())
        self.case_tree.blockSignals(False)

    def _start(self) -> None:
        if self._worker is not None and self._worker.isRunning():
            return
        selected = [
            (case_id, item.text(1))
            for case_id, item in self._case_items.items()
            if item.checkState(0) == Qt.CheckState.Checked
        ]
        if not selected:
            QMessageBox.information(self, "No mesh Cases selected", "Select at least one mesh Case.")
            return
        executable_text = self.executable_combo.currentText().strip()
        executable = Path(executable_text).expanduser().resolve() if executable_text else None
        if executable is None or not executable.is_file():
            QMessageBox.warning(self, "OneFLOW executable required", "Select a valid OneFLOW executable.")
            return
        commands = []
        for case_id, title in selected:
            arguments = (
                sys.executable, "-m", "tools.oneflow_gui.mesh_validation",
                "--oneflow-executable", str(executable),
                "--mesh-case-id", case_id,
                "--cfd-template-id", str(self.template_combo.currentData()),
                "--max-steps", str(self.steps_spin.value()),
            )
            commands.append((f"Mesh-to-CFD: {title}", arguments, self._repository_root))
        self.log_view.clear()
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self._worker = MeshValidationWorker(commands)
        self._worker.output_line.connect(self.log_view.appendPlainText)
        self._worker.completed.connect(self._completed)
        self._worker.finished.connect(self._worker_finished)
        self._worker.start()

    def _stop(self) -> None:
        if self._worker is not None:
            self._worker.request_stop()
            self.stop_button.setEnabled(False)

    def _completed(self, passed: bool, message: str) -> None:
        self.log_view.appendPlainText("\n" + message)
        self.log_view.appendPlainText("OVERALL: " + ("PASSED" if passed else "FAILED / CANCELLED"))

    def _worker_finished(self) -> None:
        self._worker = None
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)

    def closeEvent(self, event) -> None:
        if self._worker is not None and self._worker.isRunning():
            self._stop()
            event.ignore()
            return
        event.accept()
