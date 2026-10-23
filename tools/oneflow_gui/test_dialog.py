"""Qt dialog that launches the repository's existing test commands."""

from __future__ import annotations

import subprocess
import threading
from pathlib import Path

from PySide6.QtCore import QThread, Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QVBoxLayout,
)

from .user_data_paths import UserDataPaths
from .test_runner import (
    TestCommand,
    available_suites,
    build_cmake_command,
    build_default_test_commands,
    build_regression_command,
    stage_regression_tests,
)


class TestCommandWorker(QThread):
    """Execute existing test commands sequentially and stream their output."""

    output_line = Signal(str)
    completed = Signal(bool, str)

    def __init__(self, commands: list[TestCommand]) -> None:
        super().__init__()
        self._commands = commands
        self._cancel = threading.Event()
        self._process: subprocess.Popen[str] | None = None
        self._process_lock = threading.Lock()

    def request_stop(self) -> None:
        self._cancel.set()
        with self._process_lock:
            process = self._process
        if process is not None and process.poll() is None:
            process.terminate()

    def run(self) -> None:
        failed_labels: list[str] = []
        for command in self._commands:
            if self._cancel.is_set():
                self.completed.emit(False, "Test run cancelled.")
                return
            self.output_line.emit(f"\n=== {command.label} ===")
            self.output_line.emit("$ " + subprocess.list2cmdline(list(command.arguments)))
            try:
                if command.label.startswith("CMake configure"):
                    command.working_directory.mkdir(parents=True, exist_ok=True)
                with subprocess.Popen(
                    list(command.arguments),
                    cwd=command.working_directory,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    bufsize=1,
                ) as process:
                    with self._process_lock:
                        self._process = process
                    if process.stdout is not None:
                        for line in process.stdout:
                            self.output_line.emit(line.rstrip("\r\n"))
                            if self._cancel.is_set() and process.poll() is None:
                                process.terminate()
                    return_code = process.wait()
                    with self._process_lock:
                        self._process = None
            except OSError as error:
                self.output_line.emit(f"ERROR: {error}")
                failed_labels.append(command.label)
                break
            if self._cancel.is_set():
                self.completed.emit(False, "Test run cancelled.")
                return
            if return_code != 0:
                failed_labels.append(command.label)
                self.output_line.emit(f"FAILED: exit code {return_code}")
                break
            self.output_line.emit("PASSED")
        if failed_labels:
            self.completed.emit(False, "Failed: " + ", ".join(failed_labels))
        else:
            self.completed.emit(True, f"Passed all {len(self._commands)} command(s).")


class TestRunnerDialog(QDialog):
    """Offer unit tests by default and an optional existing CFD regression suite."""

    def __init__(
        self,
        repository_root: Path,
        parent=None,
        *,
        build_config: str = "Debug",
        executable: str | None = None,
    ) -> None:
        super().__init__(parent)
        if build_config not in {"Debug", "Release"}:
            raise ValueError(f"Unsupported build configuration: {build_config!r}")
        self._repository_root = Path(repository_root).resolve()
        self._build_config = build_config
        self._worker: TestCommandWorker | None = None
        self.setWindowTitle("OneFLOW Tests")
        self.resize(900, 650)

        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Run the existing unit-test commands; CFD regression is optional."))

        self.unit_tests_check = QCheckBox("Unit tests (CTest and Python tests)")
        self.unit_tests_check.setChecked(True)
        layout.addWidget(self.unit_tests_check)

        self.regression_check = QCheckBox("Optional CFD / grid regression")
        self.regression_check.setChecked(False)
        layout.addWidget(self.regression_check)

        form = QFormLayout()
        self.executable_combo = QComboBox()
        self.executable_combo.setEditable(True)
        self.executable_combo.setMinimumContentsLength(42)
        if executable:
            self.executable_combo.setEditText(executable)
        self.suite_combo = QComboBox()
        self.suite_combo.addItems(available_suites(self._repository_root))
        default_index = self.suite_combo.findText("ci-serial.txt")
        if default_index >= 0:
            self.suite_combo.setCurrentIndex(default_index)
        form.addRow("OneFLOW executable", self.executable_combo)
        form.addRow("Existing suite", self.suite_combo)
        layout.addLayout(form)

        self.log_view = QPlainTextEdit()
        self.log_view.setReadOnly(True)
        self.log_view.setMaximumBlockCount(20000)
        self.log_view.setPlaceholderText("Test output will appear here.")
        layout.addWidget(self.log_view, 1)

        buttons = QHBoxLayout()
        buttons.addStretch(1)
        self.build_button = QPushButton("Build")
        self.build_button.setToolTip("Build the selected CMake configuration.")
        self.build_button.clicked.connect(lambda: self._start_build(rebuild=False))
        self.rebuild_button = QPushButton("Rebuild")
        self.rebuild_button.setToolTip("Clean and rebuild the selected CMake configuration.")
        self.rebuild_button.clicked.connect(lambda: self._start_build(rebuild=True))
        self.run_button = QPushButton("Run tests")
        self.run_button.clicked.connect(self._start)
        self.stop_button = QPushButton("Stop")
        self.stop_button.setEnabled(False)
        self.stop_button.clicked.connect(self._stop)
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.close)
        buttons.addWidget(self.build_button)
        buttons.addWidget(self.rebuild_button)
        buttons.addWidget(self.run_button)
        buttons.addWidget(self.stop_button)
        buttons.addWidget(close_button)
        layout.addLayout(buttons)

    def _start(self) -> None:
        if self._worker is not None and self._worker.isRunning():
            return
        commands: list[TestCommand] = []
        if self.unit_tests_check.isChecked():
            commands.extend(
                build_default_test_commands(
                    self._repository_root,
                    build_config=self._build_config,
                )
            )
        if self.regression_check.isChecked():
            executable = self.executable_combo.currentText().strip()
            suite_name = self.suite_combo.currentText().strip()
            try:
                if not suite_name:
                    raise ValueError("No regression suite is available.")
                staged_test_directory = stage_regression_tests(
                    self._repository_root,
                    UserDataPaths.discover().cfd_regression_workspace,
                )
                executable_path = Path(executable).expanduser().resolve()
                expected_build_executable = (
                    self._repository_root
                    / "build"
                    / "bin"
                    / self._build_config
                    / "OneFLOW.exe"
                ).resolve()
                executable_will_be_built = (
                    executable_path == expected_build_executable
                    and any(command.label.startswith("CMake build") for command in commands)
                )
                commands.append(
                    build_regression_command(
                        self._repository_root,
                        executable,
                        suite_name,
                        test_directory=staged_test_directory,
                        validate_executable=not executable_will_be_built,
                    )
                )
            except (OSError, ValueError) as error:
                QMessageBox.warning(self, "Cannot start regression", str(error))
                return
        if not commands:
            QMessageBox.information(self, "No tests selected", "Select unit tests or CFD regression.")
            return

        self._launch_commands(commands)

    def _start_build(self, *, rebuild: bool) -> None:
        if self._worker is not None and self._worker.isRunning():
            return
        command = build_cmake_command(
            self._repository_root,
            self._build_config,
            rebuild=rebuild,
        )
        self._launch_commands([command])

    def _launch_commands(self, commands: list[TestCommand]) -> None:
        self.log_view.clear()
        self.run_button.setEnabled(False)
        self.build_button.setEnabled(False)
        self.rebuild_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self._worker = TestCommandWorker(commands)
        self._worker.output_line.connect(self.log_view.appendPlainText)
        self._worker.completed.connect(self._completed)
        self._worker.finished.connect(self._worker_finished)
        self._worker.start()

    def _stop(self) -> None:
        if self._worker is not None:
            self._worker.request_stop()
            self.stop_button.setEnabled(False)

    def _completed(self, passed: bool, message: str) -> None:
        self.log_view.appendPlainText(f"\n{message}")
        self.log_view.appendPlainText("FINAL STATUS: PASS" if passed else "FINAL STATUS: FAIL")

    def _worker_finished(self) -> None:
        self.run_button.setEnabled(True)
        self.build_button.setEnabled(True)
        self.rebuild_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self._worker = None

    def closeEvent(self, event) -> None:
        if self._worker is not None and self._worker.isRunning():
            self._worker.request_stop()
            event.ignore()
            return
        event.accept()
