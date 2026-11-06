"""Unified entry point and read-only browser for the Contract Workbench."""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import QDesktopServices, QPixmap
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QDialog,
    QDialogButtonBox,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QPlainTextEdit,
    QScrollArea,
    QSplitter,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from .application_settings import ApplicationSettings
from .source_editor import open_source_files
from .user_data_paths import UserDataPaths

from .parameter_consumer_catalog import (
    ParameterConsumerCatalog,
    validate_parameter_consumer_sources,
)


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
        self._summary_by_id = {item.id: item for item in self._summaries}
        self._diagnostics = tuple(diagnostics)
        self._parameter_consumer_catalog = (
            ParameterConsumerCatalog() if domain == "parameter" else None
        )

        self.setWindowTitle(f"Contract Workbench - {domain.title()}")
        self.resize(1120, 740)
        root = QVBoxLayout(self)
        intro = QLabel(
            f"{domain.title()} contracts use the shared Contract Workbench navigation. "
            "Evidence status describes what is declared or recorded; it is not inferred "
            "from the existence of a linked test."
        )
        intro.setWordWrap(True)
        root.addWidget(intro)

        if domain == "parameter":
            self.search_edit = QLineEdit()
            self.search_edit.setObjectName("parameterContractSearch")
            self.search_edit.setPlaceholderText(
                "Search parameters by name, description, group, source, or reference..."
            )
            self.search_edit.textChanged.connect(self._filter_contracts)
            root.addWidget(self.search_edit)
        else:
            self.search_edit = None

        split = QSplitter()
        self.contract_list = QListWidget()
        self.contract_list.setObjectName("contractList")
        self.contract_list.setMinimumWidth(300)
        self._populate_contract_list()
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
        self.implementation_page = QWidget()
        implementation_layout = QVBoxLayout(self.implementation_page)
        implementation_intro = QLabel(
            "Select a source file to inspect a bounded, read-only code excerpt. "
            "Use the source editor action below to open the complete file."
        )
        implementation_intro.setWordWrap(True)
        implementation_layout.addWidget(implementation_intro)
        implementation_split = QSplitter()
        self.implementation_file_list = QListWidget()
        self.implementation_file_list.setObjectName("implementationFileList")
        self.implementation_file_list.setMinimumWidth(220)
        self.implementation_file_list.currentRowChanged.connect(
            self._preview_implementation_source
        )
        implementation_split.addWidget(self.implementation_file_list)
        implementation_preview = QWidget()
        implementation_preview_layout = QVBoxLayout(implementation_preview)
        self.implementation_source_label = QLabel("Select a source file")
        self.implementation_source_label.setObjectName("implementationSourceLabel")
        implementation_preview_layout.addWidget(self.implementation_source_label)
        self.implementation_view = QPlainTextEdit()
        self.implementation_view.setObjectName("implementationCodePreview")
        self.implementation_view.setReadOnly(True)
        self.implementation_view.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        self.implementation_view.setStyleSheet(
            "font-family: Consolas, 'Courier New', monospace;"
        )
        implementation_preview_layout.addWidget(self.implementation_view, 1)
        implementation_navigation = QHBoxLayout()
        self.previous_implementation_chunk_button = QPushButton("Previous 80 lines")
        self.previous_implementation_chunk_button.clicked.connect(
            lambda: self._change_implementation_chunk(-80)
        )
        self.next_implementation_chunk_button = QPushButton("Next 80 lines")
        self.next_implementation_chunk_button.clicked.connect(
            lambda: self._change_implementation_chunk(80)
        )
        implementation_navigation.addWidget(self.previous_implementation_chunk_button)
        implementation_navigation.addWidget(self.next_implementation_chunk_button)
        implementation_navigation.addStretch(1)
        implementation_preview_layout.addLayout(implementation_navigation)
        implementation_split.addWidget(implementation_preview)
        implementation_split.setStretchFactor(0, 1)
        implementation_split.setStretchFactor(1, 4)
        implementation_layout.addWidget(implementation_split, 1)
        self._implementation_paths: tuple[str, ...] = ()
        self._implementation_start_line = 1
        self._implementation_chunk_size = 80
        self.open_source_button = QPushButton("Open source files...")
        self.open_source_button.clicked.connect(self._open_sources)
        self.documentation_view = QPlainTextEdit()
        self.documentation_view.setReadOnly(True)
        self.evidence_view = QPlainTextEdit()
        self.evidence_view.setReadOnly(True)
        self.examples_view = QPlainTextEdit()
        self.examples_view.setReadOnly(True)
        self.runtime_consumers_view = QPlainTextEdit()
        self.runtime_consumers_view.setReadOnly(True)
        self.case_view = QPlainTextEdit()
        self.case_view.setReadOnly(True)
        self.workflow_view = QPlainTextEdit()
        self.workflow_view.setReadOnly(True)
        self.parameter_formula_page = None
        self.parameter_formula_text = None
        self.parameter_formula_latex = None
        self.parameter_formula_preview = None
        self.parameter_formula_scroll = None
        self.render_parameter_formula_button = None
        self.copy_parameter_formula_button = None
        self._parameter_formula_pixmap = QPixmap()
        self._parameter_formula_zoom = 1.0
        if domain == "parameter":
            self.parameter_formula_page = QWidget()
            formula_layout = QVBoxLayout(self.parameter_formula_page)
            self.parameter_formula_text = QPlainTextEdit()
            self.parameter_formula_text.setReadOnly(True)
            self.parameter_formula_latex = QPlainTextEdit()
            self.parameter_formula_latex.setReadOnly(True)
            self.parameter_formula_latex.setStyleSheet("font-family: Consolas, monospace;")
            formula_layout.addWidget(QLabel("Definitions and symbol meanings"))
            formula_layout.addWidget(self.parameter_formula_text, 2)
            formula_layout.addWidget(QLabel("LaTeX source"))
            formula_layout.addWidget(self.parameter_formula_latex, 1)
            formula_actions = QHBoxLayout()
            self.render_parameter_formula_button = QPushButton("Render formula")
            self.render_parameter_formula_button.clicked.connect(self._render_parameter_formula)
            self.zoom_out_parameter_formula_button = QPushButton("Zoom out")
            self.zoom_out_parameter_formula_button.clicked.connect(
                self._zoom_parameter_formula_out
            )
            self.zoom_in_parameter_formula_button = QPushButton("Zoom in")
            self.zoom_in_parameter_formula_button.clicked.connect(
                self._zoom_parameter_formula_in
            )
            self.copy_parameter_formula_button = QPushButton("Copy LaTeX")
            self.copy_parameter_formula_button.clicked.connect(self._copy_parameter_formula)
            formula_actions.addWidget(self.render_parameter_formula_button)
            formula_actions.addWidget(self.zoom_out_parameter_formula_button)
            formula_actions.addWidget(self.zoom_in_parameter_formula_button)
            formula_actions.addWidget(self.copy_parameter_formula_button)
            formula_actions.addStretch(1)
            formula_layout.addLayout(formula_actions)
            formula_preview_label = QLabel("Formula preview")
            formula_preview_label.setObjectName("parameterFormulaPreviewLabel")
            formula_layout.addWidget(formula_preview_label)
            self.parameter_formula_preview = QLabel("Select a parameter and render its formula.")
            self.parameter_formula_preview.setAlignment(Qt.AlignmentFlag.AlignCenter)
            self.parameter_formula_preview.setMinimumSize(520, 260)
            self.parameter_formula_scroll = QScrollArea()
            self.parameter_formula_scroll.setWidgetResizable(True)
            self.parameter_formula_scroll.setMinimumHeight(300)
            self.parameter_formula_scroll.setWidget(self.parameter_formula_preview)
            formula_layout.addWidget(self.parameter_formula_scroll, 5)

        self.tabs.addTab(self.overview_view, "Overview")
        if domain == "parameter":
            self.tabs.addTab(self.parameter_formula_page, "Formula")
            self.tabs.addTab(self.implementation_page, "Implementation / Code")
            self.tabs.addTab(self.runtime_consumers_view, "Runtime Consumers")
            self.tabs.addTab(self.examples_view, "Examples")
        else:
            self.tabs.addTab(self.principle_view, "Principle")
            self.tabs.addTab(self.architecture_view, "Architecture")
            self.tabs.addTab(self.implementation_page, "Implementation / Code")
        self.tabs.addTab(self.documentation_view, "Documentation / References")
        self.tabs.addTab(self.evidence_view, "Verification")
        if domain == "case":
            self.tabs.addTab(self.case_view, "Case Setup")
        if domain == "workflow":
            self.tabs.addTab(self.workflow_view, "Workflow Graph")
        detail_layout.addWidget(self.tabs, 1)

        self.open_actions_layout = QHBoxLayout()
        self.open_source_button.setObjectName("openSourceFilesButton")
        self.open_docs_button = QPushButton("Open documentation...")
        self.open_docs_button.clicked.connect(self._open_documentation)
        self.open_tests_button = QPushButton("Open declared test sources...")
        self.open_tests_button.setEnabled(False)
        self.open_tests_button.clicked.connect(self._open_verification_sources)
        self.open_actions_layout.addWidget(self.open_source_button)
        self.open_actions_layout.addWidget(self.open_docs_button)
        self.open_actions_layout.addWidget(self.open_tests_button)
        self.open_related_formula_button = None
        if domain == "parameter":
            self.open_related_formula_button = QPushButton(
                "Open related Formula & Knowledge..."
            )
            self.open_related_formula_button.clicked.connect(
                self._open_related_formula
            )
            self.open_actions_layout.addWidget(self.open_related_formula_button)
        self.open_case_directory_button = None
        self.open_case_references_button = None
        self.open_case_runner_button = None
        self.open_workflow_definition_button = None
        self.open_workflow_runner_button = None
        if domain == "workflow":
            self.open_workflow_definition_button = QPushButton(
                "Open workflow definition..."
            )
            self.open_workflow_definition_button.clicked.connect(
                self._open_workflow_definition
            )
            self.open_actions_layout.addWidget(
                self.open_workflow_definition_button
            )
            self.open_workflow_runner_button = QPushButton(
                "Run in workflow runner..."
            )
            self.open_workflow_runner_button.clicked.connect(self._run_workflow)
            self.open_actions_layout.addWidget(self.open_workflow_runner_button)
        if domain == "case":
            self.open_case_directory_button = QPushButton("Open case directory...")
            self.open_case_directory_button.clicked.connect(self._open_case_directory)
            self.open_case_references_button = QPushButton("Open suite and references...")
            self.open_case_references_button.clicked.connect(self._open_case_references)
            self.open_case_runner_button = QPushButton("Run in case runner...")
            self.open_case_runner_button.clicked.connect(self._run_case)
            self.open_actions_layout.addWidget(self.open_case_directory_button)
            self.open_actions_layout.addWidget(self.open_case_references_button)
            self.open_actions_layout.addWidget(self.open_case_runner_button)
        self.open_actions_layout.addStretch(1)
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.accept)
        self.open_actions_layout.addWidget(close_button)
        detail_layout.addLayout(self.open_actions_layout)
        split.addWidget(detail)
        split.setStretchFactor(0, 1)
        split.setStretchFactor(1, 4)
        root.addWidget(split, 1)

        self.contract_list.currentRowChanged.connect(self._select_contract)
        if self._summaries:
            self._select_contract(self.contract_list.currentRow())
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
        from verification.execution_records import group_execution_records

        lines = ["", heading, "-" * len(heading)]
        if not records:
            lines.append("None recorded.")
            return lines
        invocation_groups = group_execution_records(records)
        for invocation_id, invocation_records in reversed(invocation_groups[-10:]):
            if invocation_id:
                outcomes = {item.outcome for item in invocation_records}
                if "failed" in outcomes:
                    outcome = "failed"
                elif "blocked" in outcomes:
                    outcome = "blocked"
                elif "cancelled" in outcomes:
                    outcome = "cancelled"
                elif outcomes == {"passed"}:
                    outcome = "passed"
                else:
                    outcome = "mixed"
                lines.append(
                    f"Invocation group: {invocation_id} "
                    f"({len(invocation_records)} steps; {outcome.upper()})"
                )
            for record in invocation_records:
                freshness = record.freshness(self._current_commit)
                lines.append(
                    f"[{record.outcome.upper()} / {freshness.upper()}] "
                    f"{record.test_name}"
                )
                lines.append(
                    f"  Finished: {record.finished_at}; framework: {record.framework}; "
                    f"commit: {record.commit or 'unknown'}"
                )
                if record.invocation_id:
                    lines.append(f"  Invocation: {record.invocation_id}")
                lines.append(f"  Command: {' '.join(record.command)}")
                for artifact in record.artifacts:
                    lines.append(
                        f"  Artifact ({artifact.role}): {artifact.path}"
                        + (f" [sha256: {artifact.sha256}]" if artifact.sha256 else "")
                    )
        return lines

    def _set_diagnostics_only(self) -> None:
        self.evidence_view.setPlainText("\n".join(self._diagnostics))

    def _populate_contract_list(self, query: str = "", selected_id: str = "") -> None:
        self.contract_list.clear()
        normalized_query = query.casefold().strip()
        if self._domain == "parameter":
            groups = {}
            for summary in self._summaries:
                searchable = " ".join((
                    summary.id, summary.title, summary.overview,
                    summary.principle, summary.architecture, summary.group_title,
                    " ".join(summary.references),
                    " ".join(summary.implementation_paths),
                    " ".join(summary.documentation_paths),
                    *(
                        " ".join((
                            formula.id, formula.title, formula.role,
                            formula.plain, formula.latex, formula.assumptions,
                            " ".join(
                                symbol.symbol + " " + symbol.meaning + " "
                                + " ".join(symbol.notations) + " "
                                + symbol.oneflow_name + " " + symbol.notes
                                for symbol in formula.symbols
                            ),
                        ))
                        for formula in (
                            summary.parameter_knowledge.formulas
                            if summary.parameter_knowledge else ()
                        )
                    ),
                    *(
                        " ".join((
                            symbol.symbol, symbol.meaning,
                            " ".join(symbol.notations), symbol.oneflow_name, symbol.notes,
                        ))
                        for symbol in (
                            summary.parameter_knowledge.symbols
                            if summary.parameter_knowledge else ()
                        )
                    ),
                    *(
                        (" ".join(summary.parameter_knowledge.references)
                         + " " + " ".join(summary.parameter_knowledge.caveats))
                        if summary.parameter_knowledge else "",
                    ),
                )).casefold()
                if normalized_query:
                    if re.fullmatch(r"[a-z_][a-z0-9_]*", normalized_query):
                        token_pattern = (
                            r"(?<![a-z0-9_])" + re.escape(normalized_query)
                            + r"(?![a-z0-9_])"
                        )
                        if re.search(token_pattern, searchable) is None:
                            continue
                    elif normalized_query not in searchable:
                        continue
                group_key = summary.group_id or "ungrouped"
                groups.setdefault(group_key, []).append(summary)

            ordered_groups = sorted(
                groups.items(),
                key=lambda pair: (
                    min(item.group_order for item in pair[1]),
                    pair[1][0].group_title.casefold() or "ungrouped",
                    pair[0],
                ),
            )
            for group_key, summaries in ordered_groups:
                group_title = summaries[0].group_title or "Other / uncategorized"
                header = QListWidgetItem(f"{group_title}  ({len(summaries)})")
                header.setData(Qt.ItemDataRole.UserRole, None)
                header.setFlags(Qt.ItemFlag.ItemIsEnabled)
                self.contract_list.addItem(header)
                for summary in sorted(summaries, key=lambda item: item.title.casefold()):
                    item = QListWidgetItem(f"    {summary.title}")
                    item.setData(Qt.ItemDataRole.UserRole, summary.id)
                    item.setToolTip(
                        f"{summary.id}\n{summary.overview}\n"
                        + "\n".join(summary.implementation_paths)
                    )
                    self.contract_list.addItem(item)
        else:
            for summary in self._summaries:
                item = QListWidgetItem(f"{summary.title}  [{summary.status}]")
                item.setData(Qt.ItemDataRole.UserRole, summary.id)
                self.contract_list.addItem(item)

        selected_row = -1
        for row in range(self.contract_list.count()):
            if self.contract_list.item(row).data(Qt.ItemDataRole.UserRole) == selected_id:
                selected_row = row
                break
        if selected_row < 0:
            candidate_ids = set()
            if self._domain == "parameter":
                for summary in self._summaries:
                    parameter_key = summary.parameter_key or summary.source_id
                    has_consumers = (
                        self._parameter_consumer_catalog is not None
                        and bool(self._parameter_consumer_catalog.find(parameter_key))
                    )
                    if summary.implementation_paths or has_consumers:
                        candidate_ids.add(summary.id)
            selected_row = next(
                (
                    row for row in range(self.contract_list.count())
                    if self.contract_list.item(row).data(Qt.ItemDataRole.UserRole)
                    in candidate_ids
                ),
                -1,
            )
            if selected_row < 0:
                selected_row = next(
                    (row for row in range(self.contract_list.count())
                     if self.contract_list.item(row).data(Qt.ItemDataRole.UserRole)),
                    -1,
                )
        self.contract_list.setCurrentRow(selected_row)

    def _filter_contracts(self, query: str) -> None:
        current = self.contract_list.currentItem()
        selected_id = current.data(Qt.ItemDataRole.UserRole) if current else ""
        self._populate_contract_list(query, selected_id or "")

    def _current_summary(self):
        item = self.contract_list.currentItem()
        if item is None:
            return None
        return self._summary_by_id.get(item.data(Qt.ItemDataRole.UserRole))

    def _populate_parameter_consumers(self, parameter_name: str, consumers) -> None:
        """Show audited source reads without implying complete consumer discovery."""
        if not consumers:
            self.runtime_consumers_view.setPlainText(
                f"No audited runtime consumer is registered for {parameter_name}.\n\n"
                "This means the parameter is outside the current curated audit slice; "
                "it does not prove that the parameter is unused. Trace task, solver, "
                "model, and lazy-read conditions before drawing that conclusion."
            )
            return

        catalog = ParameterConsumerCatalog(consumers=tuple(consumers))
        issues = validate_parameter_consumer_sources(self._repository_root, catalog)
        lines = [
            f"Parameter: {parameter_name}",
            f"Audited read sites: {len(consumers)}",
            "",
            "Evidence boundary",
            "-----------------",
            "The source check confirms that recorded read expressions, consumer symbols, "
            "and named downstream symbols still exist. It does not prove that a run "
            "entered the path or that the parameter changed a numerical result.",
            "",
        ]
        for index, consumer in enumerate(consumers, start=1):
            consumer_issues = [
                issue for issue in issues
                if issue.startswith(f"{consumer.parameter_name}:")
            ]
            lines.extend((
                f"{index}. {consumer.symbol} [{consumer.read_kind.value}]",
                f"   Direct read: {consumer.read_expression}",
                f"   Source: {consumer.source_path}",
                f"   Applicability: {consumer.applicability}",
                f"   Behavior: {consumer.behavior}",
                f"   Evidence: {consumer.evidence_note}",
                "   Source audit: "
                + ("FAILED" if consumer_issues else "references resolved"),
            ))
            for path in consumer.behavior_source_paths:
                lines.append(f"   Downstream source: {path}")
            if consumer.behavior_symbols:
                lines.append(
                    "   Downstream symbols: " + ", ".join(consumer.behavior_symbols)
                )
            if consumer.knowledge_ids:
                lines.append(
                    "   Formula & Knowledge: " + ", ".join(consumer.knowledge_ids)
                )
            lines.extend(f"   Audit issue: {issue}" for issue in consumer_issues)
            lines.append("")
        self.runtime_consumers_view.setPlainText("\n".join(lines).rstrip())


    def _select_contract(self, row: int) -> None:
        if row < 0:
            return
        item = self.contract_list.item(row)
        if item is None:
            return
        summary = self._summary_by_id.get(item.data(Qt.ItemDataRole.UserRole))
        if summary is None:
            return
        self.heading.setText(f"{summary.title} | {summary.id}")
        self.overview_view.setPlainText(
            f"Domain: {summary.domain}\n"
            f"Kind: {summary.kind}\n"
            f"Status: {summary.status}\n\n"
            f"{summary.overview}\n\n"
            + (
                "Principle\n---------\n" + summary.principle + "\n\n"
                if summary.domain == "parameter" and summary.principle.strip() else ""
            )
            + (
                "Architecture\n------------\n" + summary.architecture + "\n\n"
                if summary.domain == "parameter" and summary.architecture.strip() else ""
            )
            + "Claims\n------\n"
            + "\n".join(f"- {claim}" for claim in summary.claims)
            + "\n\nLimitations\n-----------\n"
            + "\n".join(f"- {item}" for item in summary.limitations)
        )
        self.principle_view.setPlainText(summary.principle)
        self.architecture_view.setPlainText(summary.architecture)
        parameter_key = summary.parameter_key or summary.source_id
        consumer_records = (
            self._parameter_consumer_catalog.find(parameter_key)
            if self._parameter_consumer_catalog is not None else ()
        )
        if summary.domain == "parameter":
            self._populate_parameter_consumers(parameter_key, consumer_records)
        consumer_paths = tuple(
            path
            for consumer in consumer_records
            for path in (consumer.source_path, *consumer.behavior_source_paths)
        )
        self._implementation_paths = tuple(dict.fromkeys(
            (*summary.implementation_paths, *consumer_paths)
        ))
        self._implementation_start_line = 1
        self.implementation_file_list.clear()
        for relative_path in self._implementation_paths:
            file_item = QListWidgetItem(relative_path)
            file_item.setData(Qt.ItemDataRole.UserRole, relative_path)
            self.implementation_file_list.addItem(file_item)
        first_implementation_row = 0 if self._implementation_paths else -1
        self.implementation_file_list.setCurrentRow(first_implementation_row)
        if first_implementation_row >= 0:
            # Do not rely only on currentRowChanged: the first row may already
            # be current when items are populated, leaving the preview untouched.
            self._preview_implementation_source(first_implementation_row)
        if not self._implementation_paths:
            self.implementation_source_label.setText("No implementation files declared")
            self.implementation_view.setPlainText(
                "No implementation paths are declared for this contract."
            )
            self.previous_implementation_chunk_button.setEnabled(False)
            self.next_implementation_chunk_button.setEnabled(False)
        if summary.domain == "parameter":
            if summary.examples:
                example_blocks = []
                for index, example in enumerate(summary.examples, start=1):
                    example_blocks.append(
                        f"{index}. {example.title}\n"
                        f"{'=' * (len(str(index)) + len(example.title) + 3)}\n"
                        f"Input\n-----\n{example.input}\n\n"
                        f"Calculation\n-----------\n{example.calculation}\n\n"
                        f"Result\n------\n{example.result}\n\n"
                        f"Limitations\n-----------\n{example.limitations}"
                    )
                self.examples_view.setPlainText("\n\n".join(example_blocks))
            else:
                self.examples_view.setPlainText(
                    "No reviewed worked example is registered for this parameter.\n"
                    "Add a source-backed example with input, calculation, result, and limitations."
                )
        if self.parameter_formula_text is not None:
            knowledge = summary.parameter_knowledge
            if knowledge is None:
                self.parameter_formula_text.setPlainText(
                    "No parameter-specific mathematical knowledge is registered."
                )
                self.parameter_formula_latex.setPlainText("")
                self.render_parameter_formula_button.setEnabled(False)
                self.copy_parameter_formula_button.setEnabled(False)
                self.parameter_formula_preview.setText("No formula preview available.")
                self.parameter_formula_preview.setPixmap(QPixmap())
            else:
                formula_lines = []
                latex_lines = []
                for formula in knowledge.formulas:
                    formula_lines.extend((
                        f"{formula.title} [{formula.id}]",
                        f"Role: {formula.role}",
                        f"Plain: {formula.plain}",
                        f"Assumptions: {formula.assumptions}",
                        "Terms:",
                    ))
                    for symbol in formula.symbols:
                        formula_lines.append(
                            f"  {symbol.symbol}: {symbol.meaning}; "
                            f"other notation: {', '.join(symbol.notations) or 'not specified'}; "
                            f"OneFLOW: {symbol.oneflow_name or 'no direct mapping'}"
                            + (f"; {symbol.notes}" if symbol.notes else "")
                        )
                    formula_lines.append("")
                    latex_lines.append(formula.latex)
                formula_lines.extend(("Symbol mapping", "--------------"))
                for symbol in knowledge.symbols:
                    formula_lines.append(
                        f"{symbol.symbol}: {symbol.meaning}; "
                        f"other notation: {', '.join(symbol.notations) or 'not specified'}; "
                        f"OneFLOW: {symbol.oneflow_name or 'no direct mapping'}"
                        + (f"; {symbol.notes}" if symbol.notes else "")
                    )
                formula_lines.extend(("", "References", "----------"))
                formula_lines.extend(knowledge.references)
                formula_lines.extend(("", "Caveats", "-------"))
                formula_lines.extend(knowledge.caveats)
                self.parameter_formula_text.setPlainText("\n".join(formula_lines))
                self.parameter_formula_latex.setPlainText("\n".join(latex_lines))
                has_formula = bool(latex_lines)
                self.render_parameter_formula_button.setEnabled(has_formula)
                self.copy_parameter_formula_button.setEnabled(has_formula)
                self.parameter_formula_preview.setText(
                    "Click Render formula to generate a preview."
                )
                self.parameter_formula_preview.setPixmap(QPixmap())
                self._parameter_formula_pixmap = QPixmap()
        self.open_source_button.setEnabled(bool(self._implementation_paths))
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
                            f"\n  CMake source mapping: target_resolved"
                            f"\n  CTest discovery declaration: "
                            f"{'declared' if target.discovery_declared else 'not_declared'}"
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
        if self.open_workflow_runner_button is not None:
            definition_path = (
                self._repository_root / summary.workflow_data.definition_path
                if summary.workflow_data is not None else None
            )
            self.open_workflow_runner_button.setEnabled(
                definition_path is not None and definition_path.is_file()
            )
        if self.open_related_formula_button is not None:
            self.open_related_formula_button.setEnabled(bool(summary.references))
        if summary.workflow_data is not None:
            workflow = summary.workflow_data
            lines = [
                f"Workflow ID: {workflow.workflow_id}",
                f"Execution mode: {workflow.execution_mode}",
                f"Definition: {workflow.definition_path}",
                "",
                "Steps",
                "-----",
            ]
            for step in workflow.steps:
                lines.append(f"{step.case_id} ({step.directory})")
                lines.append(
                    "  Depends on: "
                    + (", ".join(step.depends_on) if step.depends_on else "None")
                )
                for output in step.outputs:
                    lines.append(
                        f"  Produces: {output.artifact_id} -> {output.path}"
                    )
                    if output.parameter_file:
                        lines.append(
                            f"    Producer setting: {output.setting} in {output.parameter_file}"
                        )
                for item in step.inputs:
                    lines.append(
                        f"  Consumes: {item.artifact_id} -> {item.setting} "
                        f"in {item.parameter_file}"
                    )
            self.workflow_view.setPlainText("\\n".join(lines))

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
                f"Output absolute tolerance: {case.output_absolute_tolerance:g}\n"
                f"Residual absolute tolerance: {case.residual_absolute_tolerance:g}\n"
                f"Residual relative tolerance: {case.residual_relative_tolerance:g}\n"
                f"Residual NaN policy: {case.residual_nan_policy}\n\n"
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

    def _copy_parameter_formula(self) -> None:
        if self.parameter_formula_latex is not None:
            QApplication.clipboard().setText(
                self.parameter_formula_latex.toPlainText()
            )

    def _render_parameter_formula(self) -> None:
        summary = self._current_summary()
        if summary is None or summary.parameter_knowledge is None:
            return
        latex = "\n".join(
            formula.latex for formula in summary.parameter_knowledge.formulas
        )
        if not latex:
            return
        try:
            from verification.latex_preview import render_latex_preview_png
            preview_root = UserDataPaths.discover().contract_workbench_preview
            preview_root.mkdir(parents=True, exist_ok=True)
            safe_id = "".join(
                character if character.isalnum() else "_"
                for character in summary.source_id
            )
            output_path = preview_root / f"parameter_{safe_id}_formula.png"
            render_latex_preview_png(latex, output_path, dpi=300, fontsize=26)
            pixmap = QPixmap(str(output_path))
            if pixmap.isNull():
                raise RuntimeError("formula renderer did not produce a readable image")
            self._parameter_formula_pixmap = pixmap
            self._update_parameter_formula_preview()
            self.parameter_formula_preview.setToolTip(
                "Rendered from the exact LaTeX source shown above."
            )
        except Exception as error:
            self._parameter_formula_pixmap = QPixmap()
            self.parameter_formula_preview.setText(
                f"Formula preview unavailable ({error}). The LaTeX source remains available."
            )

    def _zoom_parameter_formula_in(self) -> None:
        self._change_parameter_formula_zoom(1.25)

    def _zoom_parameter_formula_out(self) -> None:
        self._change_parameter_formula_zoom(0.8)

    def _change_parameter_formula_zoom(self, factor: float) -> None:
        self._parameter_formula_zoom = min(
            3.0, max(0.5, self._parameter_formula_zoom * factor)
        )
        self._update_parameter_formula_preview()

    def _update_parameter_formula_preview(self) -> None:
        if self._parameter_formula_pixmap.isNull():
            return
        viewport = self.parameter_formula_scroll.viewport().size()
        target_width = max(
            viewport.width(),
            int(self._parameter_formula_pixmap.width() * self._parameter_formula_zoom),
        )
        target_height = max(
            viewport.height(),
            int(self._parameter_formula_pixmap.height() * self._parameter_formula_zoom),
        )
        self.parameter_formula_preview.resize(target_width, target_height)
        self.parameter_formula_preview.setPixmap(
            self._parameter_formula_pixmap.scaled(
                target_width,
                target_height,
                Qt.AspectRatioMode.KeepAspectRatio,
                Qt.TransformationMode.SmoothTransformation,
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
        selected_files = tuple((self._repository_root / relative).resolve() for relative in selection.selected_paths)
        settings = ApplicationSettings.load()
        try:
            open_source_files(settings.source_editor, selected_files, working_directory=self._repository_root)
        except (OSError, ValueError) as error:
            QMessageBox.warning(self, f"{label}: unable to launch editor",
                f"Editor: {settings.source_editor}\n{error}\n\nChoose another editor in OneFLOW Configuration.")

    def _preview_implementation_source(self, row: int) -> None:
        """Display one bounded source excerpt without opening an external editor."""
        item = self.implementation_file_list.item(row) if row >= 0 else None
        if item is None:
            return
        relative_path = item.data(Qt.ItemDataRole.UserRole)
        candidate = (self._repository_root / relative_path).resolve()
        try:
            candidate.relative_to(self._repository_root)
        except ValueError:
            self.implementation_source_label.setText("Rejected unsafe source path")
            self.implementation_view.setPlainText("The source path escapes the repository root.")
            return
        if not candidate.is_file():
            self.implementation_source_label.setText(f"{relative_path} — file not found")
            self.implementation_view.setPlainText("The declared source file does not exist.")
            return
        try:
            source_lines = candidate.read_text(encoding="utf-8").splitlines()
        except (OSError, UnicodeError) as error:
            self.implementation_source_label.setText(f"{relative_path} — unable to read")
            self.implementation_view.setPlainText(str(error))
            return
        self._implementation_start_line = max(
            1, min(self._implementation_start_line, max(1, len(source_lines)))
        )
        self._render_implementation_chunk(relative_path, source_lines)

    def _render_implementation_chunk(
        self, relative_path: str, source_lines: list[str]
    ) -> None:
        total_lines = len(source_lines)
        start = self._implementation_start_line
        end = min(total_lines, start + self._implementation_chunk_size - 1)
        numbered_lines = [
            f"{line_number:>5}  {source_lines[line_number - 1]}"
            for line_number in range(start, end + 1)
        ]
        self.implementation_source_label.setText(
            f"{relative_path}  |  Lines {start}-{end} of {total_lines}"
        )
        self.implementation_view.setPlainText(
            "\n".join(numbered_lines) if numbered_lines else "(empty source file)"
        )
        self.previous_implementation_chunk_button.setEnabled(start > 1)
        self.next_implementation_chunk_button.setEnabled(end < total_lines)

    def _change_implementation_chunk(self, offset: int) -> None:
        item = self.implementation_file_list.currentItem()
        if item is None:
            return
        self._implementation_start_line = max(
            1, self._implementation_start_line + offset
        )
        self._preview_implementation_source(self.implementation_file_list.currentRow())

    def _open_sources(self) -> None:
        summary = self._current_summary()
        if summary is not None:
            self._open_paths(self._implementation_paths, "Implementation")

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

    def _open_workflow_definition(self) -> None:
        summary = self._current_summary()
        if summary is not None and summary.workflow_data is not None:
            self._open_paths(
                (summary.workflow_data.definition_path,),
                "Workflow definition",
            )

    def _run_workflow(self) -> None:
        """Open the existing workflow runner with an explicit contract link."""
        summary = self._current_summary()
        if summary is None or summary.workflow_data is None:
            return
        definition_path = (
            self._repository_root / summary.workflow_data.definition_path
        ).resolve()
        try:
            definition_path.relative_to(self._repository_root)
            definition = json.loads(definition_path.read_text(encoding="utf-8"))
            from .workflow_dialog import WorkflowDialog
            from .workflow_repository import WorkflowDefinition

            workflow = WorkflowDefinition.from_dict(definition)
            dialog = WorkflowDialog(
                self._repository_root,
                parent=self,
                initial_workflow=workflow,
                contract_ids=(summary.id,),
            )
        except (OSError, ValueError, TypeError, json.JSONDecodeError) as error:
            QMessageBox.warning(self, "Cannot open workflow runner", str(error))
            return
        dialog.exec()

    def _open_related_formula(self) -> None:
        summary = self._current_summary()
        if summary is None or not summary.references:
            return
        from .verification_curve_dialog import VerificationCurveDialog

        dialog = VerificationCurveDialog(self._repository_root, self)
        for entry_id in summary.references:
            if dialog.select_entry(entry_id):
                dialog.exec()
                return
        QMessageBox.information(
            self,
            "Formula & Knowledge",
            "No related formula entry is currently available in the knowledge catalog.",
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

    def _run_case(self) -> None:
        """Open the residual Case runner with an explicit contract link."""
        summary = self._current_summary()
        if summary is None or summary.case_data is None:
            return
        case = summary.case_data
        try:
            from .case_contract_run_dialog import CaseContractRunDialog

            dialog = CaseContractRunDialog(
                self._repository_root,
                contract_id=summary.id,
                case_path=case.case_path,
                suite=case.suite,
                residual_baseline=case.residual_baseline,
                residual_absolute_tolerance=case.residual_absolute_tolerance,
                comparison_outputs=tuple(case.comparison_outputs),
                case_id=Path(case.case_path).name,
                parent=self,
                contract_ids=(summary.id,),
            )
        except (OSError, ValueError, TypeError) as error:
            QMessageBox.warning(self, "Cannot open case runner", str(error))
            return
        dialog.exec()


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
                "parameter",
                "Parameter",
                "Parameter definitions, runtime consumers, constraints, formula links and verification evidence.",
                self._open_parameter,
            ),
            (
                "case",
                "Case",
                "CFD case setup, suite membership, input files, output manifests and regression evidence.",
                self._open_case,
            ),
            (
                "workflow",
                "Workflow",
                "Dependency graphs, producer-consumer artifact bindings, execution ordering and failure blocking.",
                self._open_workflow,
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
            "Performance remains a planned domain; it will be added when a reproducible benchmark contract and evidence model are available."
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

    def _open_parameter(self) -> None:
        ContractDomainDialog(self._repository_root, "parameter", self).exec()

    def _open_case(self) -> None:
        ContractDomainDialog(self._repository_root, "case", self).exec()

    def _open_runtime(self) -> None:
        ContractDomainDialog(self._repository_root, "runtime", self).exec()

    def _open_workflow(self) -> None:
        ContractDomainDialog(self._repository_root, "workflow", self).exec()
