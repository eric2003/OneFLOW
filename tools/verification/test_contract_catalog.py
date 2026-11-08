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
from .cmake_test_catalog import discover_cmake_test_targets_by_source


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

    def test_stl_native_io_contract_is_available_to_mesh_workbench(self) -> None:
        definitions, issues = load_contract_catalog(self.repository_root)
        self.assertEqual(issues, (), "\\n".join(issues))

        contract = next(
            item for item in definitions if item.id == "mesh/stl-native-io"
        )
        self.assertEqual(contract.domain, "mesh")
        self.assertEqual(contract.group_id, "mesh_formats")
        self.assertEqual(contract.status, "partially_verified")
        self.assertEqual(
            set(contract.implementation_paths),
            {
                "codes/geometry/include/StlMeshIO.h",
                "codes/geometry/src/StlMeshIO.cpp",
                "codes/geometry/src/GridConversion.cpp",
                "codes/geometry/src/GridTypes.cpp",
            },
        )
        self.assertIn("doc/source/stl_mesh_io.rst", contract.documentation_paths)
        self.assertEqual(len(contract.evidence), 12)

        resolutions = resolve_contract_evidence(self.repository_root, contract)
        self.assertEqual(
            [item.status for item in resolutions],
            ["resolved"] * len(contract.evidence),
        )
        self.assertTrue(
            all(
                item.evidence.status == "linked_not_executed"
                for item in resolutions
            )
        )

        geometry_cmake = (
            self.repository_root / "tests" / "geometry" / "CMakeLists.txt"
        ).read_text(encoding="utf-8")
        self.assertIn(
            "add_oneflow_gtest(StlMeshIOTest StlMeshIOTest.cpp)",
            geometry_cmake,
        )
        self.assertIn(
            "gtest_discover_tests(${TEST_NAME} DISCOVERY_MODE PRE_TEST)",
            geometry_cmake,
        )

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
        self.assertEqual(issues, [], "\\\\n".join(issues))

    def test_process_lifecycle_contract_tests_declare_ctest_discovery(self) -> None:
        definition_path = (
            self.repository_root / "tools" / "verification" / "contracts"
            / "runtime" / "process_lifecycle.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))
        summary = adapt_declared_contract(definition)
        targets_by_source = discover_cmake_test_targets_by_source(self.repository_root)

        google_tests = [item for item in summary.evidence if item.kind == "google_test"]
        self.assertEqual(len(google_tests), 8)
        for evidence in google_tests:
            with self.subTest(test=evidence.reference):
                source_path = evidence.scope.split(";", 1)[0].strip()
                targets = targets_by_source.get(source_path, ())
                self.assertTrue(targets, f"No CMake target maps {source_path}")
                self.assertTrue(
                    any(
                        item.discovery_declared and item.discovery_mode == "PRE_TEST"
                        for item in targets
                    ),
                    f"{evidence.reference} is not discoverable by CTest",
                )
                if source_path == "tests/main/simu_imp_case_lifecycle_integration_test.cpp":
                    self.assertTrue(
                        any(
                            item.target == "simu_imp_case_lifecycle_integration_test"
                            for item in targets
                        ),
                        f"{evidence.reference} is not mapped to the production integration target",
                    )

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

    def test_parameter_max_time_contract_records_time_stop_condition(self) -> None:
        definition_path = (
            self.repository_root / "tools" / "verification" / "contracts"
            / "parameter" / "max_time.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))
        summary = adapt_declared_contract(definition)

        self.assertEqual(summary.id, "parameter/max-time")
        self.assertEqual(summary.group_id, "run_control")
        self.assertEqual(summary.group_title, "Run Control")
        self.assertIn("currTime", summary.overview)
        self.assertIn("maxSteps", summary.overview)
        self.assertEqual(
            set(summary.implementation_paths),
            {
                "codes/global/src/Ctrl.cpp",
                "codes/global/src/Iteration.cpp",
                "codes/global/include/Ctrl.h",
                "codes/global/include/Iteration.h",
                "codes/multigrid/src/Multigrid.cpp",
            },
        )
        self.assertEqual(len(summary.evidence), 1)
        self.assertEqual(summary.evidence[0].status, "declared_not_executed")
        self.assertEqual(
            resolve_contract_evidence(self.repository_root, summary)[0].status,
            "declared",
        )
        issues = validate_contract_summaries(self.repository_root, [summary])
        self.assertEqual(issues, [], "\n".join(issues))

    def test_global_time_step_seeds_previous_step_interval(self) -> None:
        source_path = (
            self.repository_root
            / "codes"
            / "global"
            / "src"
            / "Ctrl.cpp"
        )
        source = source_path.read_text(encoding="utf-8")
        initialization = source.split(
            'pdt = GetDataValue< Real >( "global_dt" );', 1
        )[1]
        self.assertIn("pdt1 = pdt;", initialization.split("maxTime", 1)[0])

        definition_path = (
            self.repository_root
            / "tools"
            / "verification"
            / "contracts"
            / "parameter"
            / "global_time_step.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))
        self.assertTrue(
            any("pdt changes" in item for item in definition["claims"]["does_not_prove"])
        )
        self.assertIn(
            "codes/unsteady/src/TimeIntegration.cpp",
            definition["implementation"]["source_paths"],
        )

    def test_parameter_cfl_ramp_contract_records_source_and_test_links(self) -> None:
        definition_path = (
            self.repository_root / "tools" / "verification" / "contracts"
            / "parameter" / "cfl_ramp.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))
        summary = adapt_declared_contract(definition)

        self.assertEqual(summary.id, "parameter/cfl-ramp")
        self.assertEqual(summary.group_id, "time_stepping")
        self.assertEqual(summary.group_title, "Time stepping")
        self.assertIn("codes/ns/include/TimeStepPolicy.h", summary.implementation_paths)
        self.assertEqual(len(summary.evidence), 9)
        self.assertTrue(
            all(item.status == "linked_not_executed" for item in summary.evidence)
        )
        resolutions = resolve_contract_evidence(self.repository_root, summary)
        self.assertEqual([item.status for item in resolutions], ["resolved"] * len(summary.evidence))
        self.assertTrue(
            all(
                item.evidence.status == "linked_not_executed"
                for item in resolutions
            )
        )
        issues = validate_contract_summaries(self.repository_root, [summary])
        self.assertEqual(issues, [], "\\n".join(issues))

    def test_dual_time_index_parameter_contract_records_optional_default(self) -> None:
        path = (
            self.repository_root / "tools" / "verification" / "contracts"
            / "parameter" / "dual_time_index.json"
        )
        definition = json.loads(path.read_text(encoding="utf-8"))
        summary = adapt_declared_contract(definition)
        self.assertEqual(summary.id, "parameter/dual-time-index")
        self.assertEqual(summary.parameter_key, "idualtime")
        self.assertEqual(summary.group_id, "time_stepping")
        self.assertIn("default to zero", summary.overview)
        self.assertTrue(summary.parameter_knowledge.formulas)
        self.assertEqual(
            validate_contract_summaries(self.repository_root, [summary]), []
        )

    def test_time_step_strategy_and_global_step_have_independent_contracts(self) -> None:
        summaries, issues = load_contract_catalog(self.repository_root)
        self.assertEqual(issues, (), "\\n".join(issues))
        parameters = {item.id: item for item in summaries if item.domain == "parameter"}
        expected = {
            "parameter/time-step-strategy": "timestepModel",
            "parameter/global-time-step": "global_dt",
        }
        for contract_id, key in expected.items():
            with self.subTest(contract=contract_id):
                summary = parameters[contract_id]
                self.assertEqual(summary.parameter_key, key)
                self.assertEqual(summary.group_id, "time_stepping")
                self.assertTrue(summary.parameter_knowledge.formulas)
                if contract_id == "parameter/global-time-step":
                    self.assertIn("MG::Run", summary.architecture)
                    self.assertIn(
                        "codes/multigrid/src/Multigrid.cpp",
                        summary.implementation_paths,
                    )
                    self.assertTrue(
                        any("pdt1" in item for item in summary.limitations)
                    )
                self.assertEqual(validate_contract_summaries(self.repository_root, [summary]), [])

    def test_parameter_max_steps_contract_records_stop_condition_scope(self) -> None:
        definition_path = (
            self.repository_root / "tools" / "verification" / "contracts"
            / "parameter" / "max_steps.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))
        summary = adapt_declared_contract(definition)

        self.assertEqual(summary.id, "parameter/max-steps")
        self.assertEqual(summary.group_id, "run_control")
        self.assertEqual(summary.group_title, "Run Control")
        self.assertEqual(summary.group_order, 10)
        self.assertIn("maxSteps", summary.overview)
        self.assertIn("iexitflag", summary.overview)
        self.assertEqual(
            set(summary.implementation_paths),
            {
                "codes/global/src/Iteration.cpp",
                "codes/global/include/Iteration.h",
                "codes/global/src/Ctrl.cpp",
            },
        )
        self.assertEqual(len(summary.evidence), 1)
        self.assertEqual(summary.evidence[0].status, "declared_not_executed")
        self.assertEqual(
            resolve_contract_evidence(self.repository_root, summary)[0].status,
            "declared",
        )
        issues = validate_contract_summaries(self.repository_root, [summary])
        self.assertEqual(issues, [], "\\n".join(issues))

    def test_parameter_mach_contract_adapts_typed_formula_knowledge(self) -> None:
        definition_path = (
            self.repository_root / "tools" / "verification" / "contracts"
            / "parameter" / "mach_ref.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))
        summary = adapt_declared_contract(definition)

        self.assertIsNotNone(summary.parameter_knowledge)
        knowledge = summary.parameter_knowledge
        self.assertEqual(len(knowledge.formulas), 4)
        self.assertEqual(knowledge.formulas[0].id, "local_mach_definition")
        self.assertEqual(knowledge.formulas[2].id, "freestream_mach_definition")
        self.assertEqual(knowledge.formulas[2].symbols[0].oneflow_name, "mach_ref")
        self.assertIn("NASA Glenn", knowledge.references[2])
        self.assertTrue(any("Local Mach" in item for item in knowledge.caveats))

    def test_time_integration_parameters_adapt_typed_knowledge(self) -> None:
        summaries, issues = load_contract_catalog(self.repository_root)
        self.assertEqual(issues, (), "\\n".join(issues))
        parameters = {
            item.id: item for item in summaries if item.domain == "parameter"
        }

        selector = parameters["parameter/time_integral"]
        self.assertEqual(selector.parameter_key, "time_integral")
        self.assertEqual(
            selector.parameter_knowledge.formulas[0].id,
            "time_integral_method_dispatch",
        )
        self.assertIn("MULTI_STAGE", selector.parameter_knowledge.formulas[0].plain)
        self.assertIn("SIMPLE", selector.parameter_knowledge.formulas[0].plain)
        self.assertIn("TimeIntegral::Simple", selector.parameter_knowledge.caveats[1])

        stage_count = parameters["parameter/rk_stage"]
        self.assertEqual(stage_count.parameter_key, "rk_stage")
        formula = stage_count.parameter_knowledge.formulas[0]
        self.assertEqual(formula.id, "oneflow_rk_stage_weight_schedule")
        self.assertIn("1 / (s + 1)", formula.plain)
        self.assertIn(
            "codes/global/include/TimeIntegralPolicy.h",
            stage_count.implementation_paths,
        )
        resolutions = resolve_contract_evidence(self.repository_root, stage_count)
        self.assertEqual([item.status for item in resolutions], ["declared", "declared"])

    def test_rejects_parameter_formula_without_latex_source(self) -> None:
        definition_path = (
            self.repository_root / "tools" / "verification" / "contracts"
            / "parameter" / "mach_ref.json"
        )
        definition = json.loads(definition_path.read_text(encoding="utf-8"))
        definition["knowledge"]["formulas"][0]["latex"] = ""

        with self.assertRaisesRegex(ValueError, "latex must be non-empty"):
            adapt_declared_contract(definition)

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
            "\\\\n".join(issues),
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

    def test_loads_source_backed_parameter_contracts(self) -> None:
        summaries, issues = load_contract_catalog(self.repository_root)
        parameters = {
            item.id: item for item in summaries if item.domain == "parameter"
        }

        source_backed_parameters = {
            contract_id: item
            for contract_id, item in parameters.items()
            if item.status != "proposed"
        }
        proposed_parameters = {
            contract_id: item
            for contract_id, item in parameters.items()
            if item.status == "proposed"
        }
        self.assertEqual(
            set(source_backed_parameters),
            {
                "parameter/time_integral",
                "parameter/rk_stage",
                "parameter/mach_ref",
                "parameter/max-steps",
                "parameter/max-time",
                "parameter/cfl-ramp",
                "parameter/dual-time-index",
                "parameter/time-step-strategy",
                "parameter/global-time-step",
                "parameter/turb_cfl_ratio",
                "parameter/rprod",
                "parameter/max_vis_ratio",
                "parameter/inflow_intensity",
                "parameter/inflow_viscosity",
            },
        )
        self.assertGreater(len(parameters), len(source_backed_parameters))
        self.assertTrue(
            all(item.parameter_key for item in proposed_parameters.values()),
            "every proposed inventory entry must retain its logical key",
        )
        self.assertEqual(parameters["parameter/max-steps"].group_id, "run_control")
        self.assertEqual(parameters["parameter/max-time"].group_id, "run_control")
        self.assertEqual(parameters["parameter/cfl-ramp"].group_id, "time_stepping")
        self.assertIn(
            "time_integral_rk_stage_weight",
            parameters["parameter/rk_stage"].references,
        )
        self.assertIn(
            "codes/ns/src/NsCom.cpp",
            parameters["parameter/mach_ref"].implementation_paths,
        )
        self.assertEqual(issues, (), "\\n".join(issues))

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


    def test_parameter_examples_are_structured_and_rendered_for_all_contracts(self) -> None:
        summaries, issues = load_contract_catalog(self.repository_root)
        self.assertFalse(issues, "\\n".join(issues))
        parameter_summaries = {
            item.id: item for item in summaries if item.domain == "parameter"
        }
        source_backed_parameters = {
            contract_id: summary
            for contract_id, summary in parameter_summaries.items()
            if summary.status != "proposed"
        }
        proposed_parameters = {
            contract_id: summary
            for contract_id, summary in parameter_summaries.items()
            if summary.status == "proposed"
        }
        self.assertEqual(len(source_backed_parameters), 14)
        self.assertGreater(len(proposed_parameters), 0)
        self.assertTrue(
            all(not item.implementation_paths and not item.evidence
                for item in proposed_parameters.values()),
            "proposed inventory entries must not claim source or execution evidence",
        )
        missing_examples = [
            contract_id
            for contract_id, summary in source_backed_parameters.items()
            if not summary.examples
        ]
        self.assertEqual(
            missing_examples,
            [],
            "parameter contracts without worked examples: "
            + ", ".join(sorted(missing_examples)),
        )
        for contract_id, summary in source_backed_parameters.items():
            with self.subTest(contract_id=contract_id):
                self.assertTrue(summary.examples)
                for example in summary.examples:
                    self.assertTrue(example.title)
                    self.assertTrue(example.input)
                    self.assertTrue(example.calculation)
                    self.assertTrue(example.result)
                    self.assertTrue(example.limitations)

    def test_turbulence_parameter_contracts_preserve_audited_scope(self) -> None:
        summaries, issues = load_contract_catalog(self.repository_root)
        self.assertEqual(issues, (), "\\n".join(issues))
        parameters = {
            item.id: item for item in summaries if item.domain == "parameter"
        }

        max_viscosity = parameters["parameter/max_vis_ratio"]
        self.assertIn("TurbCappedEddyViscosity", max_viscosity.overview)
        self.assertIn(
            "codes/turb/src/TurbViscosityPolicy.cpp",
            max_viscosity.implementation_paths,
        )
        self.assertTrue(max_viscosity.parameter_knowledge)
        self.assertIn("min(", max_viscosity.parameter_knowledge.formulas[0].plain)
        evidence_resolutions = resolve_contract_evidence(
            self.repository_root, max_viscosity
        )
        self.assertEqual(
            [item.status for item in evidence_resolutions],
            ["resolved", "resolved"],
        )

        for contract_id in (
            "parameter/inflow_intensity",
            "parameter/inflow_viscosity",
        ):
            with self.subTest(contract_id=contract_id):
                self.assertIn("ITReGama", parameters[contract_id].principle)
                self.assertTrue(parameters[contract_id].examples)

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
