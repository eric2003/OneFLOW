"""Unified entry point and read-only browser for the Contract Workbench."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (
    QAbstractItemView,
    QDialog,
    QDialogButtonBox,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QPushButton,
    QPlainTextEdit,
    QSplitter,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from .user_data_paths import UserDataPaths


def _ensure_tools_on_path(repository_root: Path) -> None:
    tools_root = (repository_root / "tools").resolve()
    tools_text = str(tools_root)
    if tools_text not in sys.path:
        sys.path.insert(0, tools_text)


class RepositoryFileSelectionDialog(QDialog):
    """Let users choose one, several, or all repository files to open."""

    def __init__(self, paths: tuple[str, ...], parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Select files to open")
        self.resize(620, 440)
        self._selected_paths: tuple[str, ...] = ()

        layout = QVBoxLayout(self)
        description = QLabel(
            "Select one or more files. Use Select all to open every listed file."
        )
        description.setWordWrap(True)
        layout.addWidget(description)

        self.file_list = QListWidget()
        self.file_list.setSelectionMode(
            QAbstractItemView.SelectionMode.ExtendedSelection
        )
        for path in paths:
            item = QListWidgetItem(path)
            item.setData(Qt.ItemDataRole.UserRole, path)
            self.file_list.addItem(item)
        layout.addWidget(self.file_list, 1)

        actions = QHBoxLayout()
        select_all_button = QPushButton("Select all")
        select_all_button.clicked.connect(self.file_list.selectAll)
        clear_button = QPushButton("Clear selection")
        clear_button.clicked.connect(self.file_list.clearSelection)
        actions.addWidget(select_all_button)
        actions.addWidget(clear_button)
        actions.addStretch(1)
        layout.addLayout(actions)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Open
            | QDialogButtonBox.StandardButton.Cancel
        )
        self.open_button = buttons.button(QDialogButtonBox.StandardButton.Open)
        self.open_button.setEnabled(False)
        self.file_list.itemSelectionChanged.connect(
            lambda: self.open_button.setEnabled(bool(self.file_list.selectedItems()))
        )
        buttons.accepted.connect(self._accept_selection)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    @property
    def selected_paths(self) -> tuple[str, ...]:
        return self._selected_paths

    def _accept_selection(self) -> None:
        self._selected_paths = tuple(
            item.data(Qt.ItemDataRole.UserRole)
            for item in self.file_list.selectedItems()
        )
        if self._selected_paths:
            self.accept()


class ContractDomainDialog(QDialog):
    """Browse contracts in one domain using the shared catalog model."""

    def __init__(
        self,
        repository_root: Path,
        domain: str,
        parent=None,
        *,
        execution_records_path: Path | None = None,
        current_commit: str | None = None,
    ) -> None:
        super().__init__(parent)
        self._repository_root = repository_root.resolve()
        self._domain = domain
        self._execution_records_path = (
            Path(execution_records_path)
            if execution_records_path is not None
            else UserDataPaths.discover().verification_execution_records
        )
        self._current_commit = (
            current_commit
            if current_commit is not None
            else self._read_repository_commit()
        )
        _ensure_tools_on_path(self._repository_root)

        from verification.cmake_test_catalog import discover_cmake_test_targets
        from verification.contract_catalog import (
            load_contract_catalog,
            resolve_contract_evidence,
        )

        self._resolve_contract_evidence = resolve_contract_evidence
        self._discover_cmake_test_targets = discover_cmake_test_targets
        summaries, diagnostics = load_contract_catalog(self._repository_root)
        self._summaries = tuple(
            item for item in summaries if item.domain == domain
        )
        self._diagnostics = tuple(diagnostics)

        self.setWindowTitle(f"Contract Workbench — {domain.title()}")
        self.resize(1120, 740)
        root = QVBoxLayout(self)
        intro = QLabel(
            f"{domain.title()} contracts use the shared Contract Workbench navigation. "
            "Evidence status describes what is declared or recorded; it is not inferred "
            "from the existence of a linked test."
        )
        intro.setWordWrap(True)
        root.addWidget(intro)

        split = QSplitter()
        self.contract_list = QListWidget()
        self.contract_list.setMinimumWidth(280)
        for summary in self._summaries:
            item = QListWidgetItem(f"{summary.title}  [{summary.status}]")
            item.setData(Qt.ItemDataRole.UserRole, summary.id)
            self.contract_list.addItem(item)
        split.addWidget(self.contract_list)

        detail = QWidget()
        detail_layout = QVBoxLayout(detail)
        self.heading = QLabel("Select a contract")
        self.heading.setObjectName("contractTitle")
        detail_layout.addWidget(self.heading)
        self.tabs = QTabWidget()
        self.overview_view = QPlainTextEdit()
        self.overview_view.setReadOnly(True)
        self.principle_view = QPlainTextEdit()
        self.principle_view.setReadOnly(True)
        self.architecture_view = QPlainTextEdit()
        self.architecture_view.setReadOnly(True)
        self.implementation_view = QPlainTextEdit()
        self.implementation_view.setReadOnly(True)
        self.implementation_page = QWidget()
        implementation_layout = QVBoxLayout(self.implementation_page)
        implementation_layout.addWidget(self.implementation_view, 1)
        self.open_source_button = QPushButton("Open selected source files...")
        self.open_source_button.clicked.connect(self._open_sources)
        implementation_layout.addWidget(
            self.open_source_button, 0, Qt.AlignmentFlag.AlignRight
        )
        self.documentation_view = QPlainTextEdit()
        self.documentation_view.setReadOnly(True)
        self.evidence_view = QPlainTextEdit()
        self.evidence_view.setReadOnly(True)
        self.case_view = QPlainTextEdit()
        self.case_view.setReadOnly(True)

        self.tabs.addTab(self.overview_view, "Overview")
        self.tabs.addTab(self.principle_view, "Principle")
        self.tabs.addTab(self.architecture_view, "Architecture")
        self.tabs.addTab(self.implementation_page, "Implementation / Code")
        self.tabs.addTab(self.documentation_view, "Documentation / References")
        self.tabs.addTab(self.evidence_view, "Verification")
        if domain == "case":
            self.tabs.addTab(self.case_view, "Case Setup")
        detail_layout.addWidget(self.tabs, 1)

        open_layout = QHBoxLayout()
        self.open_docs_button = QPushButton("Open documentation...")
        self.open_docs_button.clicked.connect(self._open_documentation)
        self.open_tests_button = QPushButton("Open declared test sources...")
        self.open_tests_button.setEnabled(False)
        self.open_tests_button.clicked.connect(self._open_verification_sources)
        open_layout.addWidget(self.open_docs_button)
        open_layout.addWidget(self.open_tests_button)
        self.open_case_directory_button = None
        self.open_case_references_button = None
        if domain == "case":
            self.open_case_directory_button = QPushButton("Open case directory...")
            self.open_case_directory_button.clicked.connect(self._open_case_directory)
            self.open_case_references_button = QPushButton("Open suite and references...")
            self.open_case_references_button.clicked.connect(self._open_case_references)
            open_layout.addWidget(self.open_case_directory_button)
            open_layout.addWidget(self.open_case_references_button)
        open_layout.addStretch(1)
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.accept)
        open_layout.addWidget(close_button)
        detail_layout.addLayout(open_layout)
        split.addWidget(detail)
        split.setStretchFactor(0, 1)
        split.setStretchFactor(1, 4)
        root.addWidget(split, 1)

        self.contract_list.currentRowChanged.connect(self._select_contract)
        if self._summaries:
            self.contract_list.setCurrentRow(0)
        else:
            self.heading.setText(f"No {domain} contracts are registered yet.")
            self._set_diagnostics_only()

        if self._diagnostics:
            self.status_label = QLabel(
                f"Catalog diagnostics: {len(self._diagnostics)} issue(s). "
                "Details are available in the Verification tab."
            )
            self.status_label.setWordWrap(True)
            root.addWidget(self.status_label)

    def _read_repository_commit(self) -> str | None:
        """Read the checked-out revision; return unknown outside a Git worktree."""
        try:
            result = subprocess.run(
                ["git", "rev-parse", "HEAD"],
                cwd=self._repository_root,
                check=True,
                capture_output=True,
                text=True,
                timeout=2,
            )
        except (OSError, subprocess.SubprocessError):
            return None
        commit = result.stdout.strip()
        return commit or None

    def _execution_history_lines(self, contract_id: str) -> list[str]:
        """Show only explicit contract links and keep unlinked runs separate."""
        from verification.execution_records import ExecutionRecordStore

        lines = ["", "Execution history", "-----------------"]
        try:
            records = ExecutionRecordStore(self._execution_records_path).list_records()
        except (OSError, ValueError) as exc:
            lines.append(f"Execution history unavailable: {exc}")
            return lines

        linked = [item for item in records if contract_id in item.contract_ids]
        unlinked = [item for item in records if not item.contract_ids]
        lines.extend(self._format_execution_records("Linked executions", linked))
        lines.extend(
            self._format_execution_records(
                "Recent unlinked executions (not attributed to this contract)",
                unlinked,
            )
        )
        return lines

    def _format_execution_records(self, heading: str, records) -> list[str]:
        lines = ["", heading, "-" * len(heading)]
        if not records:
            lines.append("None recorded.")
            return lines
        for record in reversed(records[-10:]):
            freshness = record.freshness(self._current_commit)
            lines.append(
                f"[{record.outcome.upper()} / {freshness.upper()}] "
                f"{record.test_name}"
            )
            lines.append(
                f"  Finished: {record.finished_at}; framework: {record.framework}; "
                f"commit: {record.commit or 'unknown'}"
            )
            lines.append(f"  Command: {' '.join(record.command)}")
            for artifact in record.artifacts:
                lines.append(
                    f"  Artifact ({artifact.role}): {artifact.path}"
                    + (f" [sha256: {artifact.sha256}]" if artifact.sha256 else "")
                )
        return lines

    def _set_diagnostics_only(self) -> None:
        self.evidence_view.setPlainText("\n".join(self._diagnostics))

    def _current_summary(self):
        row = self.contract_list.currentRow()
        if row < 0 or row >= len(self._summaries):
            return None
        return self._summaries[row]

    def _select_contract(self, row: int) -> None:
        if row < 0 or row >= len(self._summaries):
            return
        summary = self._summaries[row]
        self.heading.setText(f"{summary.title}  ·  {summary.id}")
        self.overview_view.setPlainText(
            f"Domain: {summary.domain}\n"
            f"Kind: {summary.kind}\n"
            f"Status: {summary.status}\n\n"
            f"{summary.overview}\n\n"
            "Claims\n------\n"
            + "\n".join(f"- {claim}" for claim in summary.claims)
            + "\n\nLimitations\n-----------\n"
            + "\n".join(f"- {item}" for item in summary.limitations)
        )
        self.principle_view.setPlainText(summary.principle)
        self.architecture_view.setPlainText(summary.architecture)
        self.implementation_view.setPlainText(
            "\n".join(summary.implementation_paths) or "No implementation paths declared."
        )
        self.open_source_button.setEnabled(bool(summary.implementation_paths))
        self.documentation_view.setPlainText(
            "References\n----------\n"
            + "\n".join(summary.references)
            + "\n\nDocumentation\n-------------\n"
            + "\n".join(summary.documentation_paths)
        )
        resolutions = self._resolve_contract_evidence(
            self._repository_root, summary
        )
        evidence_lines = []
        for resolution in resolutions:
            item = resolution.evidence
            evidence_lines.append(
                f"[{resolution.status.upper()}] {item.kind}: {item.reference}"
                f"\n  Declaration: {item.status}"
                f"\n  Scope: {item.scope}"
            )
            if resolution.resolved_path:
                evidence_lines.append(f"  Resolved file: {resolution.resolved_path}")
            if resolution.diagnostic:
                evidence_lines.append(f"  Diagnostic: {resolution.diagnostic}")
            if item.kind in ("linked_test", "google_test") and resolution.resolved_path:
                try:
                    source_relative = Path(resolution.resolved_path).resolve().relative_to(
                        self._repository_root
                    ).as_posix()
                except ValueError:
                    source_relative = ""
                if source_relative:
                    targets = self._discover_cmake_test_targets(
                        self._repository_root, source_relative
                    )
                    if not targets:
                        evidence_lines.append("  Build declaration: target_not_found")
                    for target in targets:
                        evidence_lines.append(
                            f"  Build target: {target.target}"
                            f"\n  CMake: {target.cmake_file}"
                            f"\n  Build declaration: target_resolved"
                            f"\n  CTest discovery declaration: "
                            f"{'registered' if target.discovery_registered else 'not_registered'}"
                        )
                        if target.discovery_mode:
                            evidence_lines.append(
                                f"  CTest discovery mode: {target.discovery_mode}"
                            )
        if self._diagnostics:
            evidence_lines.extend(
                ("", "Catalog diagnostics", "-------------------", *self._diagnostics)
            )
        evidence_lines.extend(self._execution_history_lines(summary.id))
        self.evidence_view.setPlainText("\n".join(evidence_lines) or "No evidence declared.")
        self.open_tests_button.setEnabled(
            any(
                resolution.evidence.kind
                in ("python_test", "google_test", "linked_test")
                and resolution.status == "resolved"
                and bool(resolution.resolved_path)
                for resolution in resolutions
            )
        )
        if summary.case_data is not None:
            case = summary.case_data
            self.case_view.setPlainText(
                f"Suite: {case.suite}\n"
                f"Case path: {case.case_path}\n"
                f"Solver: {case.solver}\n"
                f"Dimension: {case.dimension}D\n"
                f"Flow model: {case.flow_model}\n"
                f"Flux scheme: {case.flux_scheme}\n"
                f"Maximum steps: {case.max_steps}\n"
                f"Residual absolute tolerance: {case.residual_absolute_tolerance:g}\n\n"
                "Inputs\n------\n"
                + "\n".join(case.input_paths)
                + "\n\nCompared outputs\n-----------------\n"
                + "\n".join(case.comparison_outputs)
                + "\n\nResidual outputs\n-----------------\n"
                + "\n".join(case.residual_outputs)
                + "\n\nResidual baseline\n-----------------\n"
                + case.residual_baseline
                + "\n\nAnalytical references (declared; not automatically compared)\n"
                "------------------------------------------------------------\n"
                + "\n".join(
                    self._format_analytical_reference(path)
                    for path in case.analytical_references
                )
            )

    @staticmethod
    def _format_analytical_reference(relative_path: str) -> str:
        """Explain the physical meaning of a declared analytical reference."""
        name = Path(relative_path).name.lower()
        descriptions = {
            "cf.blasius": "Blasius analytical reference - skin-friction coefficient (cf)",
            "u.blasius": "Blasius analytical reference - streamwise velocity profile (u)",
            "v.blasius": "Blasius analytical reference - wall-normal velocity profile (v)",
        }
        description = descriptions.get(name)
        if description is None:
            description = f"Analytical reference: {Path(relative_path).name}"
        return f"{description}: {relative_path}"

    def _open_paths(self, paths: tuple[str, ...], label: str) -> None:
        """Resolve safe repository files, then open only the user's selection."""
        existing_paths = []
        for relative in dict.fromkeys(paths):
            candidate = (self._repository_root / relative).resolve()
            try:
                candidate.relative_to(self._repository_root)
            except ValueError:
                continue
            if candidate.is_file():
                existing_paths.append(
                    candidate.relative_to(self._repository_root).as_posix()
                )
        if not existing_paths:
            QMessageBox.information(self, label, "No existing repository files to open.")
            return

        selection = RepositoryFileSelectionDialog(tuple(existing_paths), self)
        selection.setWindowTitle(f"{label}: select files")
        if selection.exec() != QDialog.DialogCode.Accepted:
            return
        for relative in selection.selected_paths:
            candidate = (self._repository_root / relative).resolve()
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(candidate)))

    def _open_sources(self) -> None:
        summary = self._current_summary()
        if summary is not None:
            self._open_paths(summary.implementation_paths, "Implementation")

    def _open_documentation(self) -> None:
        summary = self._current_summary()
        if summary is not None:
            self._open_paths(
                summary.documentation_paths + summary.references,
                "Documentation",
            )

    def _open_verification_sources(self) -> None:
        """Open only repository files resolved from declared test evidence."""
        summary = self._current_summary()
        if summary is None:
            return
        resolutions = self._resolve_contract_evidence(
            self._repository_root, summary
        )
        relative_paths = []
        for resolution in resolutions:
            if resolution.evidence.kind not in (
                "python_test",
                "google_test",
                "linked_test",
            ):
                continue
            if resolution.status != "resolved" or not resolution.resolved_path:
                continue
            try:
                relative = Path(resolution.resolved_path).resolve().relative_to(
                    self._repository_root
                )
            except ValueError:
                continue
            relative_paths.append(relative.as_posix())
        self._open_paths(
            tuple(dict.fromkeys(relative_paths)),
            "Verification test sources",
        )

    def _open_case_directory(self) -> None:
        summary = self._current_summary()
        if summary is None or summary.case_data is None:
            return
        candidate = (self._repository_root / summary.case_data.case_path).resolve()
        try:
            candidate.relative_to(self._repository_root)
        except ValueError:
            QMessageBox.warning(
                self,
                "Case directory",
                "The declared case directory is outside the repository.",
            )
            return
        if not candidate.is_dir():
            QMessageBox.information(
                self,
                "Case directory",
                f"Case directory does not exist: {summary.case_data.case_path}",
            )
            return
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(candidate)))

    def _open_case_references(self) -> None:
        summary = self._current_summary()
        if summary is None or summary.case_data is None:
            return
        self._open_paths(
            (summary.case_data.suite, summary.case_data.residual_baseline)
            + summary.case_data.analytical_references,
            "Case references",
        )


class ContractWorkbenchDialog(QDialog):
    """Domain launcher for the shared Contract Workbench experience."""

    def __init__(self, repository_root: Path, parent=None) -> None:
        super().__init__(parent)
        self._repository_root = repository_root.resolve()
        self.setWindowTitle("Contract Workbench")
        self.resize(760, 620)

        root = QVBoxLayout(self)
        title = QLabel("Contract Workbench")
        title.setObjectName("pageTitle")
        root.addWidget(title)
        intro = QLabel(
            "Explore OneFLOW knowledge and engineering contracts through a common "
            "navigation model. Each domain keeps its own typed content and evidence."
        )
        intro.setWordWrap(True)
        root.addWidget(intro)

        self.domain_buttons: dict[str, QPushButton] = {}
        entries = (
            (
                "formula",
                "Formula & Knowledge",
                "Numerical methods, equations, theory curves, principles and references.",
                self._open_formula,
            ),
            (
                "geometry",
                "Geometry",
                "Topology invariants, geometry scenarios, fixtures and linked C++ tests.",
                self._open_geometry,
            ),
            (
                "case",
                "Case",
                "CFD case setup, suite membership, input files, output manifests and regression evidence.",
                self._open_case,
            ),
            (
                "runtime",
                "Runtime",
                "Process lifecycle, resource ownership, state transitions and verification evidence.",
                self._open_runtime,
            ),
        )
        for key, label, description, callback in entries:
            row = QHBoxLayout()
            button = QPushButton(label)
            button.setMinimumWidth(210)
            button.setMinimumHeight(54)
            button.clicked.connect(callback)
            self.domain_buttons[key] = button
            text_column = QVBoxLayout()
            description_label = QLabel(description)
            description_label.setWordWrap(True)
            text_column.addWidget(description_label)
            row.addWidget(button)
            row.addLayout(text_column, 1)
            root.addLayout(row)

        future = QLabel(
            "Workflow and Performance domains are planned as additional typed views; "
            "they will be added when their evidence models are grounded in real repository use."
        )
        future.setWordWrap(True)
        future.setObjectName("subtitle")
        root.addWidget(future)
        root.addStretch(1)

        close_button = QPushButton("Close")
        close_button.clicked.connect(self.accept)
        root.addWidget(close_button)
        self.setStyleSheet(
            """
            QDialog, QWidget { background: #f5f7fa; color: #202938; }
            QLabel#pageTitle { font-size: 24px; font-weight: 650; color: #172b4d; }
            QLabel#subtitle { color: #65738a; }
            QPushButton {
                background: white; border: 1px solid #cbd4e1; border-radius: 7px;
                padding: 10px; text-align: left; font-weight: 600;
            }
            QPushButton:hover { background: #eaf0f8; border-color: #8fa9cf; }
            """
        )

    def _open_formula(self) -> None:
        from .verification_curve_dialog import VerificationCurveDialog

        VerificationCurveDialog(self._repository_root, self).exec()

    def _open_geometry(self) -> None:
        from .geometry_contract_dialog import GeometryContractDialog

        GeometryContractDialog(self._repository_root, self).exec()

    def _open_case(self) -> None:
        ContractDomainDialog(self._repository_root, "case", self).exec()

    def _open_runtime(self) -> None:
        ContractDomainDialog(self._repository_root, "runtime", self).exec()
