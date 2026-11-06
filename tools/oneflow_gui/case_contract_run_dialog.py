"""Dialog that launches residual evidence for one Case contract."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QObject, QThread, Signal
from PySide6.QtWidgets import (
    QDialog,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QTextEdit,
    QVBoxLayout,
)

from .case_contract_run import (
    CaseContractRunRequest,
    records_path_for_repository,
    run_case_contract_and_capture,
)
from .executable_discovery import executable_candidates


class _CaseContractRunWorker(QObject):
    log_line = Signal(str)
    finished_ok = Signal(object)
    failed = Signal(str)

    def __init__(self, request: CaseContractRunRequest, records_path: Path) -> None:
        super().__init__()
        self._request = request
        self._records_path = records_path

    def run(self) -> None:
        try:
            records = run_case_contract_and_capture(
                self._request,
                records_path=self._records_path,
                log_callback=lambda line: self.log_line.emit(line),
            )
            self.finished_ok.emit(records)
        except (OSError, ValueError, FileNotFoundError) as error:
            self.failed.emit(str(error))


class CaseContractRunDialog(QDialog):
    """Pick an executable and run residual evidence linked to one Case contract."""

    def __init__(
        self,
        repository_root: Path,
        *,
        contract_id: str,
        case_path: str,
        suite: str,
        residual_baseline: str,
        residual_absolute_tolerance: float,
        comparison_outputs: tuple[str, ...] = (),
        case_id: str = "",
        parent=None,
        contract_ids: tuple[str, ...] | None = None,
    ) -> None:
        super().__init__(parent)
        self._repository_root = repository_root.resolve()
        # Explicit contract link only; never inferred from paths.
        explicit = contract_ids if contract_ids is not None else (contract_id,)
        if any(not isinstance(item, str) or not item.strip() for item in explicit):
            raise ValueError("contract_ids must contain non-empty strings")
        self._contract_ids = tuple(explicit)
        self._contract_id = self._contract_ids[0]
        self._case_path = case_path
        self._suite = suite
        self._residual_baseline = residual_baseline
        self._residual_absolute_tolerance = residual_absolute_tolerance
        self._comparison_outputs = tuple(comparison_outputs)
        self._case_id = case_id
        self._thread: QThread | None = None
        self._worker: _CaseContractRunWorker | None = None

        self.setWindowTitle(f"Case contract runner — {self._contract_id}")
        self.resize(720, 480)
        layout = QVBoxLayout(self)
        layout.addWidget(
            QLabel(
                "Runs the residual runner for this Case contract and records "
                "execution evidence with an explicit contract link. "
                "Case names and paths alone never create a contract link."
            )
        )

        form = QFormLayout()
        self.contract_label = QLabel(self._contract_id)
        self.case_label = QLabel(case_path)
        self.suite_label = QLabel(suite or "(none)")
        self.baseline_label = QLabel(residual_baseline or "(none)")
        form.addRow("Contract", self.contract_label)
        form.addRow("Case path", self.case_label)
        form.addRow("Suite", self.suite_label)
        form.addRow("Residual baseline", self.baseline_label)

        self.executable_edit = QLineEdit()
        candidates = [str(path) for path in executable_candidates(self._repository_root)]
        if candidates:
            self.executable_edit.setText(candidates[0])
        browse = QPushButton("Browse...")
        browse.clicked.connect(self._browse_executable)
        exe_row = QHBoxLayout()
        exe_row.addWidget(self.executable_edit, 1)
        exe_row.addWidget(browse)
        form.addRow("OneFLOW executable", exe_row)

        self.mpi_edit = QLineEdit()
        self.mpi_edit.setPlaceholderText("Optional, e.g. mpirun -np 1 (empty = serial)")
        form.addRow("MPI launcher", self.mpi_edit)
        layout.addLayout(form)

        self.log_view = QTextEdit()
        self.log_view.setReadOnly(True)
        layout.addWidget(self.log_view, 1)

        actions = QHBoxLayout()
        actions.addStretch(1)
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self._start_run)
        self.close_button = QPushButton("Close")
        self.close_button.clicked.connect(self.accept)
        actions.addWidget(self.run_button)
        actions.addWidget(self.close_button)
        layout.addLayout(actions)

    def _browse_executable(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select OneFLOW executable",
            str(self._repository_root),
        )
        if path:
            self.executable_edit.setText(path)

    def _start_run(self) -> None:
        executable = self.executable_edit.text().strip()
        if not executable:
            QMessageBox.warning(self, "Case contract runner", "Select an executable.")
            return
        exe_path = Path(executable).expanduser()
        if not exe_path.is_file():
            QMessageBox.warning(
                self,
                "Case contract runner",
                f"Executable not found:\n{exe_path}",
            )
            return

        request = CaseContractRunRequest(
            contract_id=self._contract_id,
            case_id=self._case_id,
            case_path=self._case_path,
            suite=self._suite,
            residual_baseline=self._residual_baseline,
            residual_absolute_tolerance=self._residual_absolute_tolerance,
            comparison_outputs=self._comparison_outputs,
            repository_root=self._repository_root,
            oneflow_executable=exe_path,
            mpi_cmd=self.mpi_edit.text().strip(),
        )
        records_path = records_path_for_repository(self._repository_root)
        self.run_button.setEnabled(False)
        self.log_view.append("Starting residual runner...")
        self._thread = QThread(self)
        self._worker = _CaseContractRunWorker(request, records_path)
        self._worker.moveToThread(self._thread)
        self._thread.started.connect(self._worker.run)
        self._worker.log_line.connect(self.log_view.append)
        self._worker.finished_ok.connect(self._on_finished)
        self._worker.failed.connect(self._on_failed)
        self._worker.finished_ok.connect(self._thread.quit)
        self._worker.failed.connect(self._thread.quit)
        self._thread.finished.connect(self._cleanup_worker)
        self._thread.start()

    def _on_finished(self, records) -> None:
        self.run_button.setEnabled(True)
        if not records:
            self.log_view.append("No execution record was written.")
            return
        record = records[0]
        self.log_view.append(
            f"Recorded {record.test_id} outcome={record.outcome} "
            f"contracts={list(record.contract_ids)}"
        )

    def _on_failed(self, message: str) -> None:
        self.run_button.setEnabled(True)
        self.log_view.append("Failed: " + message)
        QMessageBox.warning(self, "Case contract runner", message)

    def _cleanup_worker(self) -> None:
        if self._worker is not None:
            self._worker.deleteLater()
            self._worker = None
        if self._thread is not None:
            self._thread.deleteLater()
            self._thread = None
