"""Unified entry point and read-only browser for the Contract Workbench."""

from __future__ import annotations

import sys
from pathlib import Path

from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (
    QDialog,
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


def _ensure_tools_on_path(repository_root: Path) -> None:
    tools_root = (repository_root / "tools").resolve()
    tools_text = str(tools_root)
    if tools_text not in sys.path:
        sys.path.insert(0, tools_text)


class ContractDomainDialog(QDialog):
    """Browse contracts in one domain using the shared catalog model."""

    def __init__(self, repository_root: Path, domain: str, parent=None) -> None:
        super().__init__(parent)
        self._repository_root = repository_root.resolve()
        self._domain = domain
        _ensure_tools_on_path(self._repository_root)

        from verification.contract_catalog import load_contract_catalog

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
        self.documentation_view = QPlainTextEdit()
        self.documentation_view.setReadOnly(True)
        self.evidence_view = QPlainTextEdit()
        self.evidence_view.setReadOnly(True)
        self.case_view = QPlainTextEdit()
        self.case_view.setReadOnly(True)

        self.tabs.addTab(self.overview_view, "Overview")
        self.tabs.addTab(self.principle_view, "Principle")
        self.tabs.addTab(self.architecture_view, "Architecture")
        self.tabs.addTab(self.implementation_view, "Implementation / Code")
        self.tabs.addTab(self.documentation_view, "Documentation / References")
        self.tabs.addTab(self.evidence_view, "Verification")
        if domain == "case":
            self.tabs.addTab(self.case_view, "Case Setup")
        detail_layout.addWidget(self.tabs, 1)

        open_layout = QHBoxLayout()
        self.open_source_button = QPushButton("Open source files...")
        self.open_source_button.clicked.connect(self._open_sources)
        self.open_docs_button = QPushButton("Open documentation...")
        self.open_docs_button.clicked.connect(self._open_documentation)
        open_layout.addWidget(self.open_source_button)
        open_layout.addWidget(self.open_docs_button)
        self.open_case_directory_button = None
        self.open_case_references_button = None
        if domain == "case":
            self.open_case_directory_button = QPushButton("Open case directory...")
            self.open_case_directory_button.clicked.connect(self._open_case_directory)
            self.open_case_references_button = QPushButton("Open analytical references...")
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
        self.documentation_view.setPlainText(
            "References\n----------\n"
            + "\n".join(summary.references)
            + "\n\nDocumentation\n-------------\n"
            + "\n".join(summary.documentation_paths)
        )
        evidence_lines = [
            f"[{item.status}] {item.kind}: {item.reference}\n  Scope: {item.scope}"
            for item in summary.evidence
        ]
        if self._diagnostics:
            evidence_lines.extend(
                ("", "Catalog diagnostics", "-------------------", *self._diagnostics)
            )
        self.evidence_view.setPlainText("\n".join(evidence_lines) or "No evidence declared.")
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
                + "\n".join(case.analytical_references)
            )

    def _open_paths(self, paths: tuple[str, ...], label: str) -> None:
        opened = 0
        for relative in paths:
            candidate = (self._repository_root / relative).resolve()
            try:
                candidate.relative_to(self._repository_root)
            except ValueError:
                continue
            if candidate.is_file():
                QDesktopServices.openUrl(QUrl.fromLocalFile(str(candidate)))
                opened += 1
        if not opened:
            QMessageBox.information(self, label, "No existing repository files to open.")

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
            (summary.case_data.residual_baseline,)
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
