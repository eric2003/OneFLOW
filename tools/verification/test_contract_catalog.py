from __future__ import annotations

from dataclasses import replace
import codecs
import json
import subprocess
import tempfile
import sys
import unittest
from pathlib import Path

from .contract_catalog import (
    ContractEvidence,
    _read_python_numeric_constant,
    ContractSummary,
    adapt_declared_contract,
    adapt_geometry_scenario,
    adapt_theory_entry,
    load_contract_catalog,
    resolve_contract_evidence,
    validate_contract_summaries,
)
from .curve_catalog import list_theory_entries


class ContractCatalogTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.repository_root = Path(__file__).resolve().parents[2]

    def test_adapts_numerical_entry_without_claiming_execution(self) -> None:
        entry = next(
            item for item in list_theory_entries() if item.id == "sst_f2_eddy"
        )

        summary = adapt_theory_entry(entry)

        self.assertEqual(summary.id, "numerical/sst_f2_eddy")
        self.assertEqual(summary.source_id, "sst_f2_eddy")
        self.assertEqual(summary.domain, "numerical")
        self.assertEqual(summary.kind, "numerical_formula")
        self.assertTrue(summary.implementation_paths)
        self.assertTrue(summary.documentation_paths)
        self.assertTrue(summary.claims)
        self.assertTrue(summary.limitations)
        self.assertTrue(summary.evidence)
        self.assertTrue(
            all(item.status == "declared_not_executed" for item in summary.evidence)
        )

    def test_adapts_geometry_fixture_and_linked_test_as_unexecuted(self) -> None:
        scenario_path = (
            self.repository_root
            / "tools"
            / "geometry_workbench"
            / "scenarios"
            / "registered_child_resolution.json"
        )
        scenario = json.loads(scenario_path.read_text(encoding="utf-8"))

        summary = adapt_geometry_scenario(scenario)

        self.assertEqual(summary.id, "geometry/registered_child_resolution")
        self.assertEqual(summary.source_id, "registered_child_resolution")
        self.assertEqual(summary.domain, "geometry")
        self.assertEqual(summary.kind, "scenario")
        self.assertEqual(summary.evidence[0].kind, "fixture")
        self.assertEqual(summary.evidence[0].status, "fixture_available")
        self.assertIn(
            "does not execute production C++",
            summary.evidence[0].scope,
        )
        linked_tests = [
            item for item in summary.evidence if item.kind == "linked_test"
        ]
        self.assertEqual(len(linked_tests), 1)
        self.assertEqual(linked_tests[0].status, "linked_not_executed")

    def test_adapts_declared_runtime_contract_with_typed_evidence(self) -> None:
        definition_path = (
            self.repository_root
            / "tools"
            / "verification"
            / "contracts"
            / "runtime"
            / "process_lifecycle.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))

        summary = adapt_declared_contract(definition)

        self.assertEqual(summary.id, "runtime/process-lifecycle")
        self.assertEqual(summary.domain, "runtime")
        self.assertEqual(summary.kind, "state_transition")
        self.assertEqual(summary.status, "implemented")
        self.assertTrue(summary.architecture)
        self.assertEqual(len(summary.evidence), 8)
        self.assertIn(
            "SimulationEntryIntegration.FullSimulationEntryReadsControlFilesAndFinalizesProcessRuntime",
            {item.reference for item in summary.evidence},
        )
        self.assertIn(
            "SimulationEntryIntegration.FailedFullSimulationEntryTearsDownCaseAndFinalizesProcessRuntime",
            {item.reference for item in summary.evidence},
        )
        self.assertTrue(
            all(item.status == "linked_not_executed" for item in summary.evidence)
        )
        resolutions = resolve_contract_evidence(self.repository_root, summary)
        self.assertEqual(
            [item.status for item in resolutions],
            ["resolved"] * len(summary.evidence),
        )
        self.assertTrue(
            all(
                item.evidence.status == "linked_not_executed"
                for item in resolutions
            )
        )
        issues = validate_contract_summaries(self.repository_root, [summary])
        self.assertEqual(issues, [], "\n".join(issues))

    def test_adapts_gui_process_execution_contract_and_resolves_python_tests(
        self,
    ) -> None:
        definition_path = (
            self.repository_root
            / "tools"
            / "verification"
            / "contracts"
            / "runtime"
            / "gui_process_execution.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))

        summary = adapt_declared_contract(definition)
        issues = validate_contract_summaries(self.repository_root, [summary])
        resolutions = resolve_contract_evidence(self.repository_root, summary)

        self.assertEqual(summary.id, "runtime/gui-process-execution")
        self.assertEqual(len(summary.evidence), 7)
        self.assertEqual(issues, [], "\n".join(issues))
        self.assertEqual(
            [item.status for item in resolutions],
            ["resolved"] * len(summary.evidence),
        )
        self.assertTrue(
            all(item.evidence.status == "linked_not_executed" for item in resolutions)
        )

    def test_adapts_laminar_flat_plate_case_contract(self) -> None:
        definition_path = (
            self.repository_root
            / "tools"
            / "verification"
            / "contracts"
            / "case"
            / "laminar_flat_plate_slau2.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))

        summary = adapt_declared_contract(definition)

        self.assertEqual(summary.id, "case/laminar-flat-plate-slau2")
        self.assertEqual(summary.domain, "case")
        self.assertEqual(summary.kind, "cfd_regression_case")
        self.assertIsNotNone(summary.case_data)
        case_data = summary.case_data
        self.assertEqual(case_data.solver, "NsSolver")
        self.assertEqual(case_data.dimension, 2)
        self.assertEqual(case_data.flux_scheme, "SLAU2")
        self.assertEqual(case_data.max_steps, 50)
        self.assertEqual(len(case_data.comparison_outputs), 5)
        self.assertEqual(case_data.residual_outputs, ("results/res.dat",))
        self.assertEqual(len(case_data.analytical_references), 3)
        self.assertEqual(case_data.residual_absolute_tolerance, 1e-15)
        self.assertEqual(case_data.output_absolute_tolerance, 1e-8)
        self.assertEqual(case_data.residual_relative_tolerance, 0.0)
        self.assertEqual(case_data.residual_nan_policy, "reject")

        issues = validate_contract_summaries(self.repository_root, [summary])
        self.assertEqual(issues, [], "\n".join(issues))

    def test_rejects_case_dimension_outside_supported_range(self) -> None:
        definition_path = (
            self.repository_root
            / "tools"
            / "verification"
            / "contracts"
            / "case"
            / "laminar_flat_plate_slau2.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))
        definition["case"]["dimension"] = 4

        with self.assertRaisesRegex(ValueError, "dimension must be 1, 2, or 3"):
            adapt_declared_contract(definition)

    def test_rejects_non_finite_case_residual_tolerance(self) -> None:
        definition_path = (
            self.repository_root
            / "tools"
            / "verification"
            / "contracts"
            / "case"
            / "laminar_flat_plate_slau2.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))
        definition["case"]["residual_absolute_tolerance"] = float("nan")

        with self.assertRaisesRegex(ValueError, "finite and non-negative"):
            adapt_declared_contract(definition)

    def test_rejects_case_paths_that_escape_repository_root(self) -> None:
        definition_path = (
            self.repository_root
            / "tools"
            / "verification"
            / "contracts"
            / "case"
            / "laminar_flat_plate_slau2.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))
        definition["case"]["input_paths"][0] = "../outside/control.txt"

        with self.assertRaisesRegex(ValueError, "must stay within its declared root"):
            adapt_declared_contract(definition)

    def test_rejects_duplicate_case_output_declarations(self) -> None:
        definition_path = (
            self.repository_root
            / "tools"
            / "verification"
            / "contracts"
            / "case"
            / "laminar_flat_plate_slau2.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))
        definition["case"]["comparison_outputs"].append(
            definition["case"]["comparison_outputs"][0]
        )

        with self.assertRaisesRegex(ValueError, "comparison_outputs must not contain duplicate"):
            adapt_declared_contract(definition)

    def test_adapts_rae2822_dependency_workflow_contract(self) -> None:
        definition_path = (
            self.repository_root
            / "tools"
            / "verification"
            / "contracts"
            / "workflow"
            / "rae2822_grid_to_solver.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))

        summary = adapt_declared_contract(definition)

        self.assertEqual(summary.id, "workflow/rae2822-grid-to-solver")
        self.assertEqual(summary.domain, "workflow")
        self.assertEqual(summary.kind, "dependency_workflow")
        self.assertIsNotNone(summary.workflow_data)
        workflow = summary.workflow_data
        self.assertEqual(workflow.workflow_id, "rae2822-grid-to-solver")
        self.assertEqual(workflow.execution_mode, "dependency_graph")
        self.assertEqual(
            [step.case_id for step in workflow.steps],
            ["rae2822grid", "rae2822_roe_sa"],
        )
        producer, consumer = workflow.steps
        self.assertEqual(producer.outputs[0].artifact_id, "rae2822_grid")
        self.assertEqual(producer.outputs[0].path, "results/rae2822_2d.ofl")
        self.assertEqual(consumer.depends_on, ("rae2822grid",))
        self.assertEqual(consumer.inputs[0].artifact_id, "rae2822_grid")
        self.assertEqual(consumer.inputs[0].setting, "gridFileName")

        issues = validate_contract_summaries(self.repository_root, [summary])
        self.assertEqual(issues, [], "\n".join(issues))

    def test_requires_typed_payload_for_workflow_contract(self) -> None:
        definition_path = (
            self.repository_root
            / "tools"
            / "verification"
            / "contracts"
            / "workflow"
            / "rae2822_grid_to_solver.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))
        definition.pop("workflow")

        with self.assertRaisesRegex(ValueError, "requires a workflow object"):
            adapt_declared_contract(definition)

    def test_rejects_workflow_input_without_direct_artifact_dependency(self) -> None:
        definition_path = (
            self.repository_root
            / "tools"
            / "verification"
            / "contracts"
            / "workflow"
            / "rae2822_grid_to_solver.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))
        definition["workflow"]["steps"][1]["depends_on"] = []

        with self.assertRaisesRegex(ValueError, "must depend directly"):
            adapt_declared_contract(definition)

    def test_detects_workflow_definition_drift(self) -> None:
        from .contract_catalog import WorkflowContractData

        definition_path = (
            self.repository_root
            / "tools"
            / "verification"
            / "contracts"
            / "workflow"
            / "rae2822_grid_to_solver.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))
        summary = adapt_declared_contract(definition)
        source_path = (
            self.repository_root
            / "tools"
            / "oneflow_gui"
            / "examples"
            / "rae2822_multicase.json"
        )
        source_definition = json.loads(source_path.read_text(encoding="utf-8"))
        source_definition["steps"][0]["cases"][1]["directory"] = "stale-case"

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "workflow.json").write_text(
                json.dumps(source_definition), encoding="utf-8"
            )
            workflow = WorkflowContractData(
                summary.workflow_data.workflow_id,
                "workflow.json",
                summary.workflow_data.execution_mode,
                summary.workflow_data.steps,
            )
            isolated = replace(
                summary,
                references=(),
                implementation_paths=(),
                documentation_paths=(),
                evidence=(),
                workflow_data=workflow,
            )
            issues = validate_contract_summaries(root, [isolated])

        self.assertTrue(
            any("workflow definition steps do not match" in issue for issue in issues),
            "\n".join(issues),
        )

    def test_detects_case_output_tolerance_drift(self) -> None:
        definition_path = (
            self.repository_root
            / "tools"
            / "verification"
            / "contracts"
            / "case"
            / "laminar_flat_plate_slau2.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))
        definition["case"]["output_absolute_tolerance"] = 1e-7
        summary = adapt_declared_contract(definition)

        issues = validate_contract_summaries(self.repository_root, [summary])

        self.assertTrue(
            any("output tolerance does not match" in issue for issue in issues),
            "\n".join(issues),
        )

    def test_detects_case_residual_relative_tolerance_drift(self) -> None:
        definition_path = (
            self.repository_root
            / "tools"
            / "verification"
            / "contracts"
            / "case"
            / "laminar_flat_plate_slau2.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))
        definition["case"]["residual_relative_tolerance"] = 1e-12
        summary = adapt_declared_contract(definition)

        issues = validate_contract_summaries(self.repository_root, [summary])

        self.assertTrue(
            any("residual relative tolerance does not match" in issue for issue in issues),
            "\n".join(issues),
        )

    def test_rejects_case_residual_nan_policy_other_than_reject(self) -> None:
        definition_path = (
            self.repository_root
            / "tools"
            / "verification"
            / "contracts"
            / "case"
            / "laminar_flat_plate_slau2.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))
        definition["case"]["residual_nan_policy"] = "allow"

        with self.assertRaisesRegex(ValueError, "residual_nan_policy must be reject"):
            adapt_declared_contract(definition)

    def test_detects_case_residual_tolerance_drift(self) -> None:
        definition_path = (
            self.repository_root
            / "tools"
            / "verification"
            / "contracts"
            / "case"
            / "laminar_flat_plate_slau2.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))
        definition["case"]["residual_absolute_tolerance"] = 1e-8
        summary = adapt_declared_contract(definition)

        issues = validate_contract_summaries(self.repository_root, [summary])

        self.assertTrue(
            any("residual tolerance does not match" in issue for issue in issues),
            "\n".join(issues),
        )

    def test_detects_residual_history_length_drift(self) -> None:
        definition_path = (
            self.repository_root
            / "tools"
            / "verification"
            / "contracts"
            / "case"
            / "laminar_flat_plate_slau2.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))
        definition["case"]["max_steps"] = 49
        summary = adapt_declared_contract(definition)

        issues = validate_contract_summaries(self.repository_root, [summary])

        self.assertTrue(
            any("residual row count does not match" in issue for issue in issues),
            "\n".join(issues),
        )

    def test_requires_typed_payload_for_case_contract(self) -> None:
        definition_path = (
            self.repository_root
            / "tools"
            / "verification"
            / "contracts"
            / "case"
            / "laminar_flat_plate_slau2.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))
        definition.pop("case")

        with self.assertRaisesRegex(ValueError, "requires a case object"):
            adapt_declared_contract(definition)

    def test_rejects_unknown_declared_evidence_status(self) -> None:
        definition_path = (
            self.repository_root
            / "tools"
            / "verification"
            / "contracts"
            / "runtime"
            / "process_lifecycle.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))
        definition["verification"][0]["status"] = "green"

        with self.assertRaisesRegex(ValueError, "unsupported status"):
            adapt_declared_contract(definition)

    def test_rejects_unknown_declared_evidence_kind(self) -> None:
        definition_path = (
            self.repository_root
            / "tools"
            / "verification"
            / "contracts"
            / "runtime"
            / "process_lifecycle.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))
        definition["verification"][0]["kind"] = "imaginary_result"

        with self.assertRaisesRegex(ValueError, "unsupported kind"):
            adapt_declared_contract(definition)

    def test_requires_contract_id_namespace_to_match_domain(self) -> None:
        definition_path = (
            self.repository_root
            / "tools"
            / "verification"
            / "contracts"
            / "runtime"
            / "process_lifecycle.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))
        definition["id"] = "numerical/process-lifecycle"

        with self.assertRaisesRegex(ValueError, "namespaced by its domain"):
            adapt_declared_contract(definition)

    def test_rejects_unsupported_declared_contract_schema(self) -> None:
        with self.assertRaisesRegex(ValueError, "unsupported contract schema_version"):
            adapt_declared_contract({"schema_version": 2})

    def test_rejects_geometry_scenario_without_stable_id(self) -> None:
        with self.assertRaisesRegex(ValueError, "non-empty id"):
            adapt_geometry_scenario({"title": "Missing identity"})

    def test_validates_adapted_repository_references(self) -> None:
        numerical = adapt_theory_entry(
            next(item for item in list_theory_entries() if item.id == "sst_f2_eddy")
        )
        scenario_path = (
            self.repository_root
            / "tools"
            / "geometry_workbench"
            / "scenarios"
            / "registered_child_resolution.json"
        )
        scenario = json.loads(scenario_path.read_text(encoding="utf-8"))
        geometry = adapt_geometry_scenario(scenario)

        issues = validate_contract_summaries(
            self.repository_root, [numerical, geometry]
        )
        self.assertEqual(issues, [], "\n".join(issues))

    def test_validates_repository_paths_for_all_artifact_evidence_kinds(self) -> None:
        definition_path = (
            self.repository_root
            / "tools"
            / "verification"
            / "contracts"
            / "runtime"
            / "process_lifecycle.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))
        summary = adapt_declared_contract(definition)
        summary = replace(
            summary,
            evidence=summary.evidence
            + (
                ContractEvidence(
                    kind="analytical_reference",
                    reference="missing/reference-profile.dat",
                    status="declared_not_executed",
                    scope="Reference data is declared but not executed.",
                ),
                ContractEvidence(
                    kind="residual_baseline",
                    reference="missing/residual-baseline.json",
                    status="declared_not_executed",
                    scope="Residual baseline is declared but not executed.",
                ),
            ),
        )

        issues = validate_contract_summaries(self.repository_root, [summary])

        self.assertTrue(
            any("missing/reference-profile.dat" in issue for issue in issues),
            "; ".join(issues),
        )
        self.assertTrue(
            any("missing/residual-baseline.json" in issue for issue in issues),
            "; ".join(issues),
        )

    def test_keeps_opaque_geometry_fixture_ids_out_of_path_validation(self) -> None:
        summary = adapt_geometry_scenario({"id": "scenario_without_path"})
        issues = validate_contract_summaries(self.repository_root, [summary])
        self.assertEqual(issues, [], "; ".join(issues))

    def test_reports_missing_declared_output_manifest(self) -> None:
        definition_path = (
            self.repository_root
            / "tools"
            / "verification"
            / "contracts"
            / "case"
            / "laminar_flat_plate_slau2.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))
        manifest_reference = (
            "test/cfd/plateuns2dslau2/autotest/missing-compare.txt"
        )
        output_manifest = next(
            item
            for item in definition["verification"]
            if item["kind"] == "output_manifest"
        )
        output_manifest["reference"] = manifest_reference
        summary = adapt_declared_contract(definition)

        issues = validate_contract_summaries(self.repository_root, [summary])

        self.assertTrue(
            any(
                f"missing repository file: {manifest_reference}" in issue
                for issue in issues
            ),
            "\n".join(issues),
        )

    def test_importing_scenarios_does_not_load_pipeline_runner(self) -> None:
        script = (
            "import sys; "
            "import tools.geometry_workbench.scenarios; "
            "assert 'tools.geometry_workbench.pipeline_stub' not in sys.modules"
        )
        subprocess.run(
            [sys.executable, "-c", script],
            cwd=self.repository_root,
            check=True,
            capture_output=True,
            text=True,
        )

    def test_loads_complete_combined_catalog_without_reference_errors(self) -> None:
        summaries, issues = load_contract_catalog(self.repository_root)

        self.assertTrue(any(item.domain == "numerical" for item in summaries))
        self.assertTrue(any(item.domain == "geometry" for item in summaries))
        runtime_contracts = [item for item in summaries if item.domain == "runtime"]
        self.assertGreaterEqual(len(runtime_contracts), 2)
        self.assertTrue(
            any(item.id == "runtime/gui-process-execution" for item in runtime_contracts)
        )
        case_contracts = [item for item in summaries if item.domain == "case"]
        self.assertEqual(len(case_contracts), 1)
        workflow_contracts = [
            item for item in summaries if item.domain == "workflow"
        ]
        self.assertEqual(len(workflow_contracts), 1)
        self.assertEqual(
            workflow_contracts[0].id, "workflow/rae2822-grid-to-solver"
        )
        self.assertEqual(
            case_contracts[0].id,
            "case/laminar-flat-plate-slau2",
        )
        self.assertIsNotNone(case_contracts[0].case_data)
        self.assertEqual(len({item.id for item in summaries}), len(summaries))
        self.assertEqual(issues, (), "\n".join(issues))

    def test_ignores_google_test_macro_text_in_comments_and_literals(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            source = root / "CommentedContractTest.cpp"
            source.write_text(
                "// TEST(FakeSuite, CommentedBehavior) {}\n"
                "/* TEST_F(FakeSuite, BlockCommentBehavior) {} */\n"
                'const char* text = "TEST(FakeSuite, StringBehavior) {}";\n'
                'const char* raw = R"tag(TEST(FakeSuite, RawStringBehavior) {})tag";\n',
                encoding="utf-8",
            )
            summary = ContractSummary(
                id="runtime/commented-tests",
                source_id="commented-tests",
                title="Commented tests",
                domain="runtime",
                kind="test_contract",
                overview="",
                principle="",
                architecture="",
                status="implemented",
                references=(),
                implementation_paths=(),
                documentation_paths=(),
                evidence=(
                    ContractEvidence(
                        kind="google_test",
                        reference="FakeSuite.CommentedBehavior",
                        status="linked_not_executed",
                        scope="CommentedContractTest.cpp",
                    ),
                    ContractEvidence(
                        kind="google_test",
                        reference="FakeSuite.BlockCommentBehavior",
                        status="linked_not_executed",
                        scope="CommentedContractTest.cpp",
                    ),
                    ContractEvidence(
                        kind="google_test",
                        reference="FakeSuite.StringBehavior",
                        status="linked_not_executed",
                        scope="CommentedContractTest.cpp",
                    ),
                    ContractEvidence(
                        kind="google_test",
                        reference="FakeSuite.RawStringBehavior",
                        status="linked_not_executed",
                        scope="CommentedContractTest.cpp",
                    ),
                ),
                claims=(),
                limitations=(),
            )

            issues = validate_contract_summaries(root, [summary])

        self.assertEqual(len(issues), 4, "\n".join(issues))
        self.assertTrue(
            all("linked GoogleTest declaration not found" in issue for issue in issues),
            "\n".join(issues),
        )

    def test_resolves_evidence_without_promoting_execution_status(self) -> None:
        from .contract_catalog import resolve_contract_evidence

        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            source = root / "EvidenceTest.cpp"
            source.write_text("TEST(EvidenceSuite, Present) {}\n", encoding="utf-8")
            summary = ContractSummary(
                id="runtime/evidence-resolution",
                source_id="evidence-resolution",
                title="Evidence resolution",
                domain="runtime",
                kind="test_contract",
                overview="",
                principle="",
                architecture="",
                status="implemented",
                references=(),
                implementation_paths=(),
                documentation_paths=(),
                evidence=(
                    ContractEvidence(
                        kind="google_test",
                        reference="EvidenceSuite.Present",
                        status="linked_not_executed",
                        scope="EvidenceTest.cpp; synthetic test source.",
                    ),
                    ContractEvidence(
                        kind="analytical_reference",
                        reference="missing.dat",
                        status="declared_not_executed",
                        scope="Reference profile.",
                    ),
                    ContractEvidence(
                        kind="pure_api",
                        reference="ComputeFlux",
                        status="declared_not_executed",
                        scope="API name only.",
                    ),
                ),
                claims=(),
                limitations=(),
            )

            resolutions = resolve_contract_evidence(root, summary)

        self.assertEqual(
            [item.status for item in resolutions], ["resolved", "stale", "declared"]
        )
        self.assertEqual(
            [item.evidence.status for item in resolutions],
            ["linked_not_executed", "declared_not_executed", "declared_not_executed"],
        )
        self.assertTrue(resolutions[0].resolved_path.endswith("EvidenceTest.cpp"))
        self.assertIn("missing repository file", resolutions[1].diagnostic)


    def test_reads_python_sources_with_utf8_bom(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            source = root / "BomContractTest.py"
            source.write_bytes(
                codecs.BOM_UTF8
                + (
                    "OUTPUT_ABSOLUTE_TOLERANCE = 1e-8\n"
                    "class BomContractTest(unittest.TestCase):\n"
                    "    def test_bom_source(self):\n"
                    "        pass\n"
                ).encode("utf-8")
            )

            tolerance = _read_python_numeric_constant(
                root, "BomContractTest.py", "OUTPUT_ABSOLUTE_TOLERANCE"
            )
            summary = ContractSummary(
                id="runtime/python-bom-source",
                source_id="python-bom-source",
                title="Python BOM source",
                domain="runtime",
                kind="test_contract",
                overview="",
                principle="",
                architecture="",
                status="implemented",
                references=(),
                implementation_paths=(),
                documentation_paths=(),
                evidence=(
                    ContractEvidence(
                        kind="python_test",
                        reference="BomContractTest.test_bom_source",
                        status="linked_not_executed",
                        scope="BomContractTest.py; UTF-8 BOM regression fixture.",
                    ),
                ),
                claims=(),
                limitations=(),
            )

            resolutions = resolve_contract_evidence(root, summary)
            issues = validate_contract_summaries(root, [summary])

        self.assertEqual(tolerance, 1e-8)
        self.assertEqual([item.status for item in resolutions], ["resolved"])
        self.assertEqual(issues, [], "\n".join(issues))

    def test_resolves_python_unittest_declarations_without_importing_them(self) -> None:
        from .contract_catalog import resolve_contract_evidence

        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            source = root / "ProcessControlTest.py"
            source.write_text(
                "class ProcessControlTest(unittest.TestCase):\n"
                "    def test_termination_contract(self):\n"
                "        pass\n"
                "    async def test_async_contract(self):\n"
                "        pass\n",
                encoding="utf-8",
            )
            summary = ContractSummary(
                id="runtime/python-test-resolution",
                source_id="python-test-resolution",
                title="Python test resolution",
                domain="runtime",
                kind="test_contract",
                overview="",
                principle="",
                architecture="",
                status="implemented",
                references=(),
                implementation_paths=(),
                documentation_paths=(),
                evidence=(
                    ContractEvidence(
                        kind="python_test",
                        reference="ProcessControlTest.test_termination_contract",
                        status="linked_not_executed",
                        scope="ProcessControlTest.py; test source.",
                    ),
                    ContractEvidence(
                        kind="python_test",
                        reference="ProcessControlTest.test_missing",
                        status="linked_not_executed",
                        scope="ProcessControlTest.py; test source.",
                    ),
                    ContractEvidence(
                        kind="python_test",
                        reference="ProcessControlTest.test_async_contract",
                        status="linked_not_executed",
                        scope="ProcessControlTest.py; test source.",
                    ),
                ),
                claims=(),
                limitations=(),
            )

            resolutions = resolve_contract_evidence(root, summary)

        self.assertEqual(
            [item.status for item in resolutions],
            ["resolved", "stale", "resolved"],
        )
        self.assertIn("Python unittest declaration not found", resolutions[1].diagnostic)

    def test_resolves_typed_google_test_declarations(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            source = root / "TypedContractTest.cpp"
            source.write_text(
                "TYPED_TEST(TypedSuite, TypedBehavior) {}\n"
                "TYPED_TEST_P(RegisteredTypedSuite, RegisteredBehavior) {}\n",
                encoding="utf-8",
            )
            summary = ContractSummary(
                id="runtime/typed-tests",
                source_id="typed-tests",
                title="Typed tests",
                domain="runtime",
                kind="test_contract",
                overview="",
                principle="",
                architecture="",
                status="implemented",
                references=(),
                implementation_paths=(),
                documentation_paths=(),
                evidence=(
                    ContractEvidence(
                        kind="google_test",
                        reference="TypedSuite.TypedBehavior",
                        status="linked_not_executed",
                        scope="TypedContractTest.cpp",
                    ),
                    ContractEvidence(
                        kind="google_test",
                        reference="RegisteredTypedSuite.RegisteredBehavior",
                        status="linked_not_executed",
                        scope="TypedContractTest.cpp",
                    ),
                ),
                claims=(),
                limitations=(),
            )

            issues = validate_contract_summaries(root, [summary])

        self.assertEqual(issues, [], "; ".join(issues))

    def test_reports_stale_linked_test_declaration(self) -> None:
        scenario_path = (
            self.repository_root
            / "tools"
            / "geometry_workbench"
            / "scenarios"
            / "registered_child_resolution.json"
        )
        scenario = json.loads(scenario_path.read_text(encoding="utf-8"))
        scenario["test"]["name"] = "MissingProductionTest"
        summary = adapt_geometry_scenario(scenario)
        summary = ContractSummary(
            id=summary.id,
            source_id=summary.source_id,
            title=summary.title,
            domain=summary.domain,
            kind=summary.kind,
            overview=summary.overview,
            principle=summary.principle,
            architecture=summary.architecture,
            status=summary.status,
            references=summary.references,
            implementation_paths=summary.implementation_paths,
            documentation_paths=summary.documentation_paths,
            evidence=(
                summary.evidence[0],
                ContractEvidence(
                    kind="linked_test",
                    reference="FaceSearchBasicTest.MissingProductionTest",
                    status="linked_not_executed",
                    scope="tests/geometry/FaceSearchBasicTest.cpp",
                ),
            ),
            claims=summary.claims,
            limitations=summary.limitations,
        )

        issues = validate_contract_summaries(self.repository_root, [summary])

        self.assertTrue(
            any("linked GoogleTest declaration not found" in issue for issue in issues),
            "\n".join(issues),
        )

    def test_reports_missing_and_escaping_paths(self) -> None:
        summary = ContractSummary(
            id="test/path_check",
            source_id="path_check",
            title="Path check",
            domain="test",
            kind="test",
            overview="",
            principle="",
            architecture="",
            status="implemented",
            references=(),
            implementation_paths=(
                "missing/file.cpp",
                "../outside.cpp",
                "C:\\outside\\source.cpp",
            ),
            documentation_paths=(),
            evidence=(),
            claims=(),
            limitations=(),
        )

        issues = validate_contract_summaries(self.repository_root, [summary])

        self.assertEqual(len(issues), 3)
        self.assertTrue(any("missing repository file" in issue for issue in issues))
        self.assertTrue(any("must stay relative to root" in issue for issue in issues))

    def test_reports_duplicate_contract_ids(self) -> None:
        summary = ContractSummary(
            id="test/duplicate",
            source_id="duplicate",
            title="Duplicate",
            domain="test",
            kind="test",
            overview="",
            principle="",
            architecture="",
            status="implemented",
            references=(),
            implementation_paths=(),
            documentation_paths=(),
            evidence=(),
            claims=(),
            limitations=(),
        )

        issues = validate_contract_summaries(
            self.repository_root, [summary, summary]
        )

        self.assertEqual(issues, ["duplicate contract id: test/duplicate"])


if __name__ == "__main__":
    unittest.main()
