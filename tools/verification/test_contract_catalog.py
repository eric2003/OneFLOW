from __future__ import annotations

from dataclasses import replace
import json
import subprocess
import tempfile
import sys
import unittest
from pathlib import Path

from .contract_catalog import (
    ContractEvidence,
    ContractSummary,
    adapt_declared_contract,
    adapt_geometry_scenario,
    adapt_theory_entry,
    load_contract_catalog,
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
        self.assertEqual(len(summary.evidence), 3)
        self.assertTrue(
            all(item.status == "linked_not_executed" for item in summary.evidence)
        )
        issues = validate_contract_summaries(self.repository_root, [summary])
        self.assertEqual(issues, [], "\n".join(issues))

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

        issues = validate_contract_summaries(self.repository_root, [summary])
        self.assertEqual(issues, [], "\n".join(issues))

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
        self.assertTrue(any(item.domain == "runtime" for item in summaries))
        case_contracts = [item for item in summaries if item.domain == "case"]
        self.assertEqual(len(case_contracts), 1)
        self.assertEqual(
            case_contracts[0].id,
            "case/laminar-flat-plate-slau2",
        )
        self.assertIsNotNone(case_contracts[0].case_data)
        self.assertEqual(len({item.id for item in summaries}), len(summaries))
        self.assertEqual(issues, (), "\n".join(issues))

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
