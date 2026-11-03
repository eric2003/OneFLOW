"""Tests for the unified Contract Workbench domain navigation."""

from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication

from tools.verification.execution_records import (
    ExecutionArtifact,
    ExecutionRecord,
    ExecutionRecordStore,
)

from tools.oneflow_gui.contract_workbench_dialog import (
    ContractDomainDialog,
    ContractWorkbenchDialog,
)


class ContractWorkbenchDialogTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.application = QApplication.instance() or QApplication([])
        cls.repository_root = Path(__file__).resolve().parents[2]

    def test_domain_hub_exposes_current_contract_domains(self) -> None:
        dialog = ContractWorkbenchDialog(self.repository_root)
        try:
            self.assertEqual(
                set(dialog.domain_buttons),
                {"formula", "geometry", "parameter", "case", "runtime", "workflow"},
            )
            self.assertEqual(
                dialog.domain_buttons["formula"].text(),
                "Formula & Knowledge",
            )
            self.assertEqual(dialog.domain_buttons["parameter"].text(), "Parameter")
            self.assertEqual(dialog.domain_buttons["case"].text(), "Case")
            self.assertEqual(dialog.domain_buttons["runtime"].text(), "Runtime")
            self.assertEqual(dialog.domain_buttons["workflow"].text(), "Workflow")
        finally:
            dialog.close()

    def test_parameter_domain_lists_source_backed_parameter_contracts(self) -> None:
        dialog = ContractDomainDialog(self.repository_root, "parameter")
        try:
            contract_ids = {
                dialog.contract_list.item(index).data(Qt.ItemDataRole.UserRole)
                for index in range(dialog.contract_list.count())
            }
            self.assertTrue(
                {
                    "parameter/time_integral",
                    "parameter/rk_stage",
                    "parameter/mach_ref",
                }.issubset(contract_ids)
            )
            dialog.contract_list.setCurrentRow(
                next(
                    index
                    for index in range(dialog.contract_list.count())
                    if dialog.contract_list.item(index).data(256)
                    == "parameter/rk_stage"
                )
            )
            self.assertIn("parameter/rk_stage", dialog.heading.text())
            self.assertIn(
                "time_integral_rk_stage_weight",
                dialog.documentation_view.toPlainText(),
            )
            self.assertIn("Ctrl.cpp", dialog.implementation_view.toPlainText())
            self.assertIn("declared_not_executed", dialog.evidence_view.toPlainText())
            self.assertTrue(dialog.open_related_formula_button.isEnabled())
            with patch(
                "tools.oneflow_gui.verification_curve_dialog.VerificationCurveDialog"
            ) as formula_dialog_type:
                formula_dialog_type.return_value.select_entry.return_value = True
                dialog._open_related_formula()
            formula_dialog_type.return_value.select_entry.assert_called_once_with(
                "time_integral_rk_stage_weight"
            )
            formula_dialog_type.return_value.exec.assert_called_once()
        finally:
            dialog.close()

    def test_mach_parameter_exposes_formula_symbols_aliases_and_latex(self) -> None:
        dialog = ContractDomainDialog(self.repository_root, "parameter")
        try:
            mach_row = next(
                index for index in range(dialog.contract_list.count())
                if dialog.contract_list.item(index).data(Qt.ItemDataRole.UserRole)
                == "parameter/mach_ref"
            )
            dialog.contract_list.setCurrentRow(mach_row)
            self.assertEqual(
                dialog.tabs.tabText(dialog.tabs.indexOf(dialog.parameter_formula_page)),
                "Formula",
            )
            tab_names = [
                dialog.tabs.tabText(index) for index in range(dialog.tabs.count())
            ]
            self.assertEqual(tab_names[:2], ["Overview", "Formula"])
            self.assertNotIn("Principle", tab_names)
            self.assertNotIn("Architecture", tab_names)
            self.assertNotIn("Implementation / Code", tab_names)
            self.assertIn("NsCom.cpp", dialog.overview_view.toPlainText())
            self.assertNotIn("[implemented]", dialog.contract_list.item(mach_row).text())
            self.assertTrue(dialog.zoom_in_parameter_formula_button.isEnabled())
            self.assertTrue(dialog.zoom_out_parameter_formula_button.isEnabled())
            self.assertGreaterEqual(
                dialog.parameter_formula_scroll.minimumHeight(), 300
            )
            formula_layout = dialog.parameter_formula_page.layout()
            preview_label = dialog.parameter_formula_page.findChild(
                type(dialog.parameter_formula_preview),
                "parameterFormulaPreviewLabel",
            )
            self.assertIsNotNone(preview_label)
            self.assertLess(
                formula_layout.indexOf(dialog.parameter_formula_latex),
                formula_layout.indexOf(preview_label),
            )
            self.assertLess(
                formula_layout.indexOf(preview_label),
                formula_layout.indexOf(dialog.parameter_formula_scroll),
            )
            formula_text = dialog.parameter_formula_text.toPlainText()
            latex_text = dialog.parameter_formula_latex.toPlainText()
            self.assertIn("Local Mach number", formula_text)
            self.assertIn("Free-stream / reference Mach number", formula_text)
            self.assertIn("gama_ref", formula_text)
            self.assertIn("M_{\\infty}", latex_text)
            self.assertIn("const_cp", formula_text)
            self.assertIn(r"\mathrm{const\_cp}", latex_text)
            self.assertTrue(dialog.render_parameter_formula_button.isEnabled())
            self.assertTrue(dialog.copy_parameter_formula_button.isEnabled())

            dialog.search_edit.setText("Ma_ref")
            visible_ids = {
                dialog.contract_list.item(index).data(Qt.ItemDataRole.UserRole)
                for index in range(dialog.contract_list.count())
                if dialog.contract_list.item(index).data(Qt.ItemDataRole.UserRole)
            }
            self.assertEqual(visible_ids, {"parameter/mach_ref"})
        finally:
            dialog.close()

    def test_parameter_contracts_are_grouped_and_searchable(self) -> None:
        dialog = ContractDomainDialog(self.repository_root, "parameter")
        try:
            self.assertIsNotNone(dialog.search_edit)
            headers = [
                dialog.contract_list.item(index).text()
                for index in range(dialog.contract_list.count())
                if dialog.contract_list.item(index).data(Qt.ItemDataRole.UserRole) is None
            ]
            self.assertIn("Run Control  (3)", headers)
            self.assertIn("Time Integration  (2)", headers)
            self.assertIn("Compressible Flow  (1)", headers)

            dialog.search_edit.setText("mach_ref")
            visible_ids = {
                dialog.contract_list.item(index).data(Qt.ItemDataRole.UserRole)
                for index in range(dialog.contract_list.count())
                if dialog.contract_list.item(index).data(Qt.ItemDataRole.UserRole)
            }
            self.assertEqual(visible_ids, {"parameter/mach_ref"})
            dialog.search_edit.setText("Maximum iteration steps")
            visible_ids = {
                dialog.contract_list.item(index).data(Qt.ItemDataRole.UserRole)
                for index in range(dialog.contract_list.count())
                if dialog.contract_list.item(index).data(Qt.ItemDataRole.UserRole)
            }
            self.assertEqual(visible_ids, {"parameter/max-steps"})
            dialog.search_edit.setText("Maximum simulation time")
            visible_ids = {
                dialog.contract_list.item(index).data(Qt.ItemDataRole.UserRole)
                for index in range(dialog.contract_list.count())
                if dialog.contract_list.item(index).data(Qt.ItemDataRole.UserRole)
            }
            self.assertEqual(visible_ids, {"parameter/max-time"})
            dialog.search_edit.setText("maxTime")
            visible_ids = {
                dialog.contract_list.item(index).data(Qt.ItemDataRole.UserRole)
                for index in range(dialog.contract_list.count())
                if dialog.contract_list.item(index).data(Qt.ItemDataRole.UserRole)
            }
            self.assertEqual(
                visible_ids, {"parameter/max-steps", "parameter/max-time"}
            )
            dialog.search_edit.setText("NsCom.cpp")
            visible_ids = {
                dialog.contract_list.item(index).data(Qt.ItemDataRole.UserRole)
                for index in range(dialog.contract_list.count())
                if dialog.contract_list.item(index).data(Qt.ItemDataRole.UserRole)
            }
            self.assertEqual(
                visible_ids, {"parameter/mach_ref", "parameter/cfl-ramp"}
            )
            dialog.search_edit.clear()
            visible_ids = {
                dialog.contract_list.item(index).data(Qt.ItemDataRole.UserRole)
                for index in range(dialog.contract_list.count())
                if dialog.contract_list.item(index).data(Qt.ItemDataRole.UserRole)
            }
            self.assertEqual(visible_ids, {
                "parameter/time_integral", "parameter/rk_stage", "parameter/mach_ref",
                "parameter/max-steps", "parameter/max-time", "parameter/cfl-ramp"
            })
        finally:
            dialog.close()

    def test_case_domain_displays_typed_case_setup_tab(self) -> None:
        dialog = ContractDomainDialog(self.repository_root, "case")
        try:
            self.assertGreaterEqual(dialog.contract_list.count(), 1)
            self.assertEqual(dialog.tabs.tabText(dialog.tabs.count() - 1), "Case Setup")
            self.assertIn("SLAU2", dialog.case_view.toPlainText())
            self.assertIn("Output absolute tolerance: 1e-08", dialog.case_view.toPlainText())
            self.assertIn("Residual relative tolerance: 0", dialog.case_view.toPlainText())
            self.assertIn("Residual NaN policy: reject", dialog.case_view.toPlainText())
            self.assertIn("Blasius", dialog.case_view.toPlainText())
            case_text = dialog.case_view.toPlainText()
            self.assertIn(
                "Blasius analytical reference - skin-friction coefficient (cf)",
                case_text,
            )
            self.assertIn(
                "Blasius analytical reference - streamwise velocity profile (u)",
                case_text,
            )
            self.assertIn(
                "Blasius analytical reference - wall-normal velocity profile (v)",
                case_text,
            )
            self.assertIsNotNone(dialog.open_case_directory_button)
            self.assertIsNotNone(dialog.open_case_references_button)
            desktop_services = Mock()
            with patch(
                "tools.oneflow_gui.contract_workbench_dialog.QDesktopServices",
                desktop_services,
            ):
                dialog._open_case_directory()
            expected_directory = (
                self.repository_root / "test/cfd/plateuns2dslau2"
            ).resolve()
            desktop_services.openUrl.assert_called_once()
            self.assertEqual(
                Path(
                    desktop_services.openUrl.call_args.args[0].toLocalFile()
                ).resolve(),
                expected_directory,
            )
            selected_references = (
                "test/suites/cpu-serial.txt",
                "test/baselines/residual-baseline.json",
                "test/cfd/plateuns2dslau2/results/cf.blasius",
                "test/cfd/plateuns2dslau2/results/u.blasius",
                "test/cfd/plateuns2dslau2/results/v.blasius",
            )
            with patch(
                "tools.oneflow_gui.contract_workbench_dialog.RepositoryFileSelectionDialog"
            ) as selection_type:
                selection = selection_type.return_value
                selection.exec.return_value = 1
                selection.selected_paths = selected_references
                with patch(
                    "tools.oneflow_gui.contract_workbench_dialog.ApplicationSettings.load",
                    return_value=Mock(source_editor="vscode"),
                ), patch(
                    "tools.oneflow_gui.contract_workbench_dialog.open_source_files"
                ) as open_files:
                    dialog._open_case_references()
            open_files.assert_called_once_with(
                "vscode",
                tuple((self.repository_root / item).resolve() for item in selected_references),
                working_directory=self.repository_root,
            )
        finally:
            dialog.close()

    def test_workflow_domain_displays_typed_graph_and_definition_navigation(self) -> None:
        dialog = ContractDomainDialog(self.repository_root, "workflow")
        try:
            contract_index = next(
                index
                for index, summary in enumerate(dialog._summaries)
                if summary.id == "workflow/rae2822-grid-to-solver"
            )
            dialog.contract_list.setCurrentRow(contract_index)
            self.assertEqual(
                dialog.tabs.tabText(dialog.tabs.count() - 1), "Workflow Graph"
            )
            graph_text = dialog.workflow_view.toPlainText()
            self.assertIn("rae2822-grid-to-solver", graph_text)
            self.assertIn(
                "Produces: rae2822_grid -> results/rae2822_2d.ofl", graph_text
            )
            self.assertIn(
                "Consumes: rae2822_grid -> gridFileName in script/cfd.txt",
                graph_text,
            )
            self.assertIn("Depends on: rae2822grid", graph_text)
            self.assertIsNotNone(dialog.open_workflow_definition_button)
            self.assertTrue(dialog.open_workflow_definition_button.isEnabled())
            self.assertIsNotNone(dialog.open_workflow_runner_button)
            self.assertTrue(dialog.open_workflow_runner_button.isEnabled())
            self.assertTrue(dialog.open_tests_button.isEnabled())

            definition_path = "tools/oneflow_gui/examples/rae2822_multicase.json"
            with patch(
                "tools.oneflow_gui.contract_workbench_dialog.RepositoryFileSelectionDialog"
            ) as selection_type:
                selection = selection_type.return_value
                selection.exec.return_value = 1
                selection.selected_paths = (definition_path,)
                with patch(
                    "tools.oneflow_gui.contract_workbench_dialog.ApplicationSettings.load",
                    return_value=Mock(source_editor="vscode"),
                ), patch(
                    "tools.oneflow_gui.contract_workbench_dialog.open_source_files"
                ) as open_files:
                    dialog._open_workflow_definition()
            open_files.assert_called_once_with(
                "vscode",
                ((self.repository_root / definition_path).resolve(),),
                working_directory=self.repository_root,
            )
        finally:
            dialog.close()

    def test_workflow_runner_receives_explicit_contract_link(self) -> None:
        dialog = ContractDomainDialog(self.repository_root, "workflow")
        try:
            contract_index = next(
                index for index, summary in enumerate(dialog._summaries)
                if summary.id == "workflow/rae2822-grid-to-solver"
            )
            dialog.contract_list.setCurrentRow(contract_index)
            with patch("tools.oneflow_gui.workflow_dialog.WorkflowDialog") as runner_type:
                runner_type.return_value.exec.return_value = 0
                dialog._run_workflow()
            runner_type.assert_called_once()
            args, kwargs = runner_type.call_args
            self.assertEqual(args[0], self.repository_root)
            self.assertEqual(kwargs["contract_ids"], ("workflow/rae2822-grid-to-solver",))
            self.assertEqual(kwargs["initial_workflow"].workflow_id, "rae2822-grid-to-solver")
            self.assertEqual(kwargs["initial_workflow"].steps[0]["type"], "case_graph")
        finally:
            dialog.close()

    def test_runtime_verification_shows_resolution_separately_from_declaration(self) -> None:
        dialog = ContractDomainDialog(self.repository_root, "runtime")
        try:
            contract_index = next(
                index
                for index, summary in enumerate(dialog._summaries)
                if summary.id == "runtime/process-lifecycle"
            )
            dialog.contract_list.setCurrentRow(contract_index)
            evidence_text = dialog.evidence_view.toPlainText()
            self.assertIn("[RESOLVED] google_test:", evidence_text)
            self.assertIn("Declaration: linked_not_executed", evidence_text)
            self.assertIn("Resolved file:", evidence_text)
            self.assertIn(
                "Build target: simu_imp_case_lifecycle_integration_test",
                evidence_text,
            )
            self.assertIn("Build declaration: target_resolved", evidence_text)
            self.assertIn("CTest discovery declaration: registered", evidence_text)
            self.assertIn("CTest discovery mode: PRE_TEST", evidence_text)
        finally:
            dialog.close()

    def test_gui_process_contract_displays_python_test_declarations(self) -> None:
        dialog = ContractDomainDialog(self.repository_root, "runtime")
        try:
            contract_index = next(
                index
                for index, summary in enumerate(dialog._summaries)
                if summary.id == "runtime/gui-process-execution"
            )
            dialog.contract_list.setCurrentRow(contract_index)
            evidence_text = dialog.evidence_view.toPlainText()
            self.assertIn("[RESOLVED] python_test:", evidence_text)
            self.assertIn(
                "ProcessControlTest.test_terminate_and_reap_waits_for_direct_child",
                evidence_text,
            )
            self.assertIn("Declaration: linked_not_executed", evidence_text)
            self.assertIn("Resolved file:", evidence_text)
        finally:
            dialog.close()

    def test_gui_process_contract_exposes_implementation_sources(self) -> None:
        dialog = ContractDomainDialog(self.repository_root, "runtime")
        try:
            contract_index = next(
                index
                for index, summary in enumerate(dialog._summaries)
                if summary.id == "runtime/gui-process-execution"
            )
            dialog.contract_list.setCurrentRow(contract_index)
            implementation_text = dialog.implementation_view.toPlainText()
            self.assertIn("tools/oneflow_gui/process_control.py", implementation_text)
            self.assertIn("tools/oneflow_gui/test_dialog.py", implementation_text)
            self.assertFalse(dialog.open_source_button.isHidden())
            self.assertEqual(
                dialog.open_actions_layout.indexOf(dialog.open_source_button), 0
            )
            self.assertLess(
                dialog.open_actions_layout.indexOf(dialog.open_source_button),
                dialog.open_actions_layout.indexOf(dialog.open_docs_button),
            )
            self.assertLess(
                dialog.open_actions_layout.indexOf(dialog.open_docs_button),
                dialog.open_actions_layout.indexOf(dialog.open_tests_button),
            )
            self.assertEqual(dialog.open_source_button.text(), "Open source files...")
            self.assertTrue(dialog.open_source_button.isEnabled())
        finally:
            dialog.close()

    def test_file_selection_supports_single_multiple_and_select_all(self) -> None:
        from tools.oneflow_gui.contract_workbench_dialog import (
            RepositoryFileSelectionDialog,
        )

        dialog = RepositoryFileSelectionDialog(
            ("src/one.cpp", "src/two.cpp", "src/three.cpp")
        )
        try:
            dialog.file_list.item(1).setSelected(True)
            self.assertEqual(
                [item.text() for item in dialog.file_list.selectedItems()],
                ["src/two.cpp"],
            )
            dialog.file_list.item(2).setSelected(True)
            self.assertEqual(
                {item.text() for item in dialog.file_list.selectedItems()},
                {"src/two.cpp", "src/three.cpp"},
            )
            dialog.file_list.clearSelection()
            select_all = next(
                button
                for button in dialog.findChildren(type(dialog.open_button))
                if button.text() == "Select all"
            )
            select_all.click()
            selected = {
                item.text() for item in dialog.file_list.selectedItems()
            }
            self.assertEqual(selected, {"src/one.cpp", "src/two.cpp", "src/three.cpp"})
            dialog._accept_selection()
            self.assertEqual(
                set(dialog.selected_paths),
                {"src/one.cpp", "src/two.cpp", "src/three.cpp"},
            )
        finally:
            dialog.close()

    def test_open_sources_opens_only_selected_files(self) -> None:
        dialog = ContractDomainDialog(self.repository_root, "runtime")
        try:
            contract_index = next(
                index
                for index, summary in enumerate(dialog._summaries)
                if summary.id == "runtime/gui-process-execution"
            )
            dialog.contract_list.setCurrentRow(contract_index)
            with patch(
                "tools.oneflow_gui.contract_workbench_dialog.RepositoryFileSelectionDialog"
            ) as selection_type:
                selection = selection_type.return_value
                selection.exec.return_value = 1
                selection.selected_paths = ("tools/oneflow_gui/process_control.py",)
                with patch("tools.oneflow_gui.contract_workbench_dialog.ApplicationSettings.load",
                           return_value=Mock(source_editor="vscode")), patch(
                    "tools.oneflow_gui.contract_workbench_dialog.open_source_files"
                ) as open_files:
                    dialog._open_sources()
            open_files.assert_called_once_with("vscode",
                ((self.repository_root / "tools/oneflow_gui/process_control.py").resolve(),),
                working_directory=self.repository_root)
        finally:
            dialog.close()

    def test_open_paths_reports_editor_launch_failures(self) -> None:
        dialog = ContractDomainDialog(self.repository_root, "runtime")
        try:
            with patch(
                "tools.oneflow_gui.contract_workbench_dialog.RepositoryFileSelectionDialog"
            ) as selection_type:
                selection = selection_type.return_value
                selection.exec.return_value = 1
                selection.selected_paths = ("tools/oneflow_gui/process_control.py",)
                with patch("tools.oneflow_gui.contract_workbench_dialog.ApplicationSettings.load",
                           return_value=Mock(source_editor="vscode")), patch(
                    "tools.oneflow_gui.contract_workbench_dialog.open_source_files",
                    side_effect=OSError("VS Code executable was not found."),
                ), patch(
                    "tools.oneflow_gui.contract_workbench_dialog.QMessageBox.warning"
                ) as warning:
                    dialog._open_paths(
                        ("tools/oneflow_gui/process_control.py",), "Implementation"
                    )
            warning.assert_called_once()
            self.assertIn("VS Code executable was not found", warning.call_args.args[2])
        finally:
            dialog.close()

    def test_open_declared_python_test_sources_from_verification(self) -> None:
        dialog = ContractDomainDialog(self.repository_root, "runtime")
        try:
            contract_index = next(
                index
                for index, summary in enumerate(dialog._summaries)
                if summary.id == "runtime/gui-process-execution"
            )
            dialog.contract_list.setCurrentRow(contract_index)
            self.assertTrue(dialog.open_tests_button.isEnabled())
            expected_paths = (
                "tools/oneflow_gui/test_process_control.py",
                "tools/oneflow_gui/test_runner_process.py",
                "tools/oneflow_gui/test_workflow_process.py",
            )
            with patch(
                "tools.oneflow_gui.contract_workbench_dialog.RepositoryFileSelectionDialog"
            ) as selection_type:
                selection = selection_type.return_value
                selection.exec.return_value = 1
                selection.selected_paths = expected_paths
                with patch(
                    "tools.oneflow_gui.contract_workbench_dialog.ApplicationSettings.load",
                    return_value=Mock(source_editor="vs2026"),
                ), patch(
                    "tools.oneflow_gui.contract_workbench_dialog.open_source_files"
                ) as open_files:
                    dialog.open_tests_button.click()
            open_files.assert_called_once_with(
                "vs2026",
                tuple((self.repository_root / item).resolve() for item in expected_paths),
                working_directory=self.repository_root,
            )
        finally:
            dialog.close()

    def test_execution_history_separates_explicit_links_from_unlinked_runs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            store_path = Path(temporary) / "execution-records.jsonl"
            store = ExecutionRecordStore(store_path)
            common = {
                "test_id": "unittest:sample",
                "framework": "unittest",
                "test_name": "SampleTest.test_contract",
                "source_path": "tools/sample_test.py",
                "command": ("python", "-m", "unittest", "tools.sample_test"),
                "outcome": "passed",
                "started_at": "2026-11-02T01:00:00+08:00",
                "finished_at": "2026-11-02T01:00:01+08:00",
                "exit_code": 0,
                "branch": "hexin",
            }
            store.append(
                ExecutionRecord(
                    **common,
                    commit="current-sha",
                    contract_ids=("runtime/gui-process-execution",),
                    artifacts=(ExecutionArtifact("log", ".oneflow/logs/run.log"),),
                )
            )
            store.append(
                ExecutionRecord(
                    **{
                        **common,
                        "test_id": "unittest:unlinked",
                        "test_name": "SampleTest.test_unlinked",
                    },
                    commit="old-sha",
                    contract_ids=(),
                )
            )
            dialog = ContractDomainDialog(
                self.repository_root,
                "runtime",
                execution_records_path=store_path,
                current_commit="current-sha",
            )
            try:
                index = next(
                    index
                    for index, summary in enumerate(dialog._summaries)
                    if summary.id == "runtime/gui-process-execution"
                )
                dialog.contract_list.setCurrentRow(index)
                text = dialog.evidence_view.toPlainText()
                linked_section = text.split("Linked executions", 1)[1].split(
                    "Recent unlinked executions", 1
                )[0]
                unlinked_section = text.split(
                    "Recent unlinked executions", 1
                )[1]
                self.assertIn("SampleTest.test_contract", linked_section)
                self.assertIn("CURRENT", linked_section)
                self.assertIn(".oneflow/logs/run.log", linked_section)
                self.assertNotIn("unlinked", linked_section)
                self.assertIn("SampleTest.test_unlinked", unlinked_section)
                self.assertIn("STALE", unlinked_section)
                self.assertIn("not attributed to this contract", text)
            finally:
                dialog.close()

    def test_execution_history_groups_workflow_steps_by_invocation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            store = ExecutionRecordStore(Path(temporary) / "execution-records.jsonl")
            common = {
                "framework": "workflow",
                "source_path": "tools/oneflow_gui/workflow_dialog.py",
                "command": ("python", "-m", "tools.oneflow_gui"),
                "outcome": "passed",
                "started_at": "2026-11-02T01:00:00+08:00",
                "finished_at": "2026-11-02T01:00:01+08:00",
                "exit_code": 0,
                "commit": "current-sha",
                "branch": "hexin",
                "contract_ids": ("workflow/rae2822-grid-to-solver",),
                "invocation_id": "workflow-run-42",
            }
            store.append(ExecutionRecord(
                **common, test_id="workflow:grid", test_name="Grid.run"
            ))
            store.append(ExecutionRecord(
                **common, test_id="workflow:solver", test_name="Solver.run"
            ))
            dialog = ContractDomainDialog(
                self.repository_root, "workflow",
                execution_records_path=store.path,
                current_commit="current-sha",
            )
            try:
                index = next(
                    index for index, summary in enumerate(dialog._summaries)
                    if summary.id == "workflow/rae2822-grid-to-solver"
                )
                dialog.contract_list.setCurrentRow(index)
                evidence = dialog.evidence_view.toPlainText()
                self.assertIn(
                    "Invocation group: workflow-run-42 (2 steps; PASSED)", evidence
                )
                self.assertLess(evidence.index("Grid.run"), evidence.index("Solver.run"))
            finally:
                dialog.close()

    def test_runtime_domain_uses_shared_contract_navigation(self) -> None:
        dialog = ContractDomainDialog(self.repository_root, "runtime")
        try:
            self.assertGreaterEqual(dialog.contract_list.count(), 1)
            tab_names = [
                dialog.tabs.tabText(index)
                for index in range(dialog.tabs.count())
            ]
            self.assertIn("Overview", tab_names)
            self.assertIn("Principle", tab_names)
            self.assertIn("Architecture", tab_names)
            self.assertIn("Implementation / Code", tab_names)
            self.assertIn("Verification", tab_names)
            self.assertNotIn("Case Setup", tab_names)
            self.assertIsNone(dialog.open_case_directory_button)
            self.assertIsNone(dialog.open_case_references_button)
        finally:
            dialog.close()


if __name__ == "__main__":
    unittest.main()
