"""Tests for source-based test inventory discovery."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.test_inventory import audit_python_test_registry, discover_tests


class TestInventoryDiscoveryTest(unittest.TestCase):
    def test_discovers_cpp_tests_and_links_contract_data(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            cpp = root / "tests" / "geometry" / "ExampleTest.cpp"
            cpp.parent.mkdir(parents=True)
            cpp.write_text(
                "TEST(ExampleSuite, FindsIdentity) {}\n"
                "TEST_F(ExampleFixture, KeepsMapping) {}\n",
                encoding="utf-8",
            )
            scenario_dir = root / "tools" / "geometry_workbench" / "scenarios"
            scenario_dir.mkdir(parents=True)
            (scenario_dir / "identity.json").write_text(
                json.dumps({
                    "id": "identity",
                    "test": {
                        "framework": "GoogleTest",
                        "suite": "ExampleSuite",
                        "name": "FindsIdentity",
                    },
                    "input": {"candidateIds": [3, 5]},
                    "expected": {"registeredIds": [5]},
                }),
                encoding="utf-8",
            )

            records = discover_tests(root)

            self.assertEqual(len(records), 2)
            record = next(item for item in records if item["name"] == "FindsIdentity")
            self.assertEqual(record["source"], "tests/geometry/ExampleTest.cpp")
            self.assertEqual(record["line"], 1)
            self.assertEqual(record["scenarioId"], "identity")
            self.assertEqual(record["input"], {"candidateIds": [3, 5]})
            self.assertEqual(record["expected"], {"registeredIds": [5]})
            self.assertIsNone(record["result"])

    def test_unlinked_test_has_no_fabricated_input_or_expected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            cpp = root / "tests" / "SampleTest.cpp"
            cpp.parent.mkdir(parents=True)
            cpp.write_text("TEST(Sample, NoFixture) {}\n", encoding="utf-8")

            record = discover_tests(root)[0]

            self.assertEqual(record["id"], "gtest:Sample.NoFixture")
            self.assertIsNone(record["input"])
            self.assertIsNone(record["expected"])
            self.assertIsNone(record["result"])
            self.assertEqual(record["discovery"], "source")

    def test_discovers_python_unittest_methods(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            gui = root / "tools" / "oneflow_gui"
            gui.mkdir(parents=True)
            (gui / "test_sample.py").write_text(
                "import unittest\n"
                "class SampleTest(unittest.TestCase):\n"
                "    def test_behavior(self):\n"
                "        self.assertTrue(True)\n",
                encoding="utf-8",
            )

            record = discover_tests(root)[0]

            self.assertEqual(record["framework"], "unittest")
            self.assertEqual(record["suite"], "SampleTest")
            self.assertEqual(record["name"], "test_behavior")
            self.assertEqual(record["line"], 3)
            self.assertIn("tools.oneflow_gui.test_sample.SampleTest.test_behavior",
                          record["execution"])

    def test_discovery_order_is_stable(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            cpp = root / "tests" / "SampleTest.cpp"
            cpp.parent.mkdir(parents=True)
            cpp.write_text(
                "TEST(Sample, Zulu) {}\nTEST(Sample, Alpha) {}\n",
                encoding="utf-8",
            )

            first = discover_tests(root)
            second = discover_tests(root)

            self.assertEqual([item["id"] for item in first],
                             [item["id"] for item in second])
            self.assertEqual(len(first), 2)



    def test_discovers_contract_catalog_tests_for_workbench(self) -> None:
        repository_root = Path(__file__).resolve().parents[2]

        records = discover_tests(repository_root)

        contract_tests = [
            item for item in records
            if item["source"] == "tools/verification/test_contract_catalog.py"
        ]
        self.assertTrue(contract_tests)
        self.assertTrue(
            any(
                item["name"] == "test_loads_complete_combined_catalog_without_reference_errors"
                for item in contract_tests
            )
        )
        self.assertTrue(
            any(
                item["name"]
                == "test_adapts_gui_process_execution_contract_and_resolves_python_tests"
                for item in contract_tests
            )
        )

    def test_discovers_workflow_contract_and_workbench_tests(self) -> None:
        repository_root = Path(__file__).resolve().parents[2]

        records = discover_tests(repository_root)
        expected = {
            (
                "tools/verification/test_contract_catalog.py",
                "test_adapts_rae2822_dependency_workflow_contract",
            ),
            (
                "tools/verification/test_contract_catalog.py",
                "test_detects_workflow_definition_drift",
            ),
            (
                "tools/oneflow_gui/test_contract_workbench_dialog.py",
                "test_workflow_domain_displays_typed_graph_and_definition_navigation",
            ),
        }
        discovered = {(item["source"], item["name"]) for item in records}

        self.assertTrue(expected.issubset(discovered), sorted(expected - discovered))

    def test_discovers_verification_execution_record_tests(self) -> None:
        repository_root = Path(__file__).resolve().parents[2]

        records = discover_tests(repository_root)

        execution_tests = [
            item for item in records
            if item["source"] == "tools/verification/test_execution_records.py"
        ]
        self.assertTrue(execution_tests)
        self.assertTrue(
            any(
                item["name"] == "test_jsonl_store_appends_and_reads_multiple_records"
                for item in execution_tests
            )
        )
        self.assertTrue(
            any(
                item["name"] == "test_freshness_requires_exact_revision_match"
                for item in execution_tests
            )
        )

    def test_discovers_workflow_execution_evidence_tests(self) -> None:
        repository_root = Path(__file__).resolve().parents[2]
        records = discover_tests(repository_root)
        workflow_tests = [
            item for item in records
            if item["source"] == "tools/oneflow_gui/test_workflow_execution_capture.py"
        ]
        names = {item["name"] for item in workflow_tests}
        self.assertIn(
            "test_captures_launched_workflow_step_with_contract_and_artifact_evidence",
            names,
        )
        self.assertIn(
            "test_does_not_record_dependency_blocked_steps_as_process_runs",
            names,
        )
        self.assertIn(
            "test_worker_preserves_source_and_staged_execution_metadata",
            names,
        )

    def test_discovers_python_execution_capture_tests(self) -> None:
        repository_root = Path(__file__).resolve().parents[2]

        records = discover_tests(repository_root)

        capture_tests = [
            item for item in records
            if item["source"] == "tools/oneflow_gui/test_execution_capture.py"
        ]
        self.assertTrue(capture_tests)
        self.assertTrue(
            any(
                item["name"] == "test_captures_completed_unittest_command_and_log_artifact"
                for item in capture_tests
            )
        )
        self.assertTrue(
            any(
                item["name"] == "test_records_failed_and_cancelled_process_outcomes"
                for item in capture_tests
            )
        )

    def test_discovers_test_worker_execution_evidence_tests(self) -> None:
        repository_root = Path(__file__).resolve().parents[2]

        records = discover_tests(repository_root)

        worker_tests = [
            item for item in records
            if item["source"] == "tools/oneflow_gui/test_test_dialog_execution_capture.py"
        ]
        self.assertTrue(worker_tests)
        self.assertTrue(
            any(
                item["name"] == "test_worker_persists_evidence_after_unittest_process_completes"
                for item in worker_tests
            )
        )
        self.assertTrue(
            any(
                item["name"] == "test_process_launch_failure_does_not_create_execution_evidence"
                for item in worker_tests
            )
        )

    def test_geometry_contracts_expose_input_expected_and_unexecuted_result(self) -> None:
        repository_root = Path(__file__).resolve().parents[2]
        records = discover_tests(repository_root)
        by_name = {record["name"]: record for record in records}

        registered = by_name[
            "RegisteredChildResolutionPreservesCandidateIndexForPartialMatch"
        ]
        self.assertEqual(registered["framework"], "GoogleTest")
        self.assertEqual(
            registered["input"]["candidateFaces"],
            [[10, 20], [20, 30]],
        )
        self.assertEqual(registered["expected"]["candidateIndices"], [1])
        self.assertIsNone(registered["result"])

        split = by_name["SplitOneSharedParentViaTablePipeline"]
        self.assertEqual(split["framework"], "GoogleTest")
        self.assertEqual(split["input"]["parents"], [0, 1])
        self.assertEqual(split["expected"]["l2g"], [[2, 3, 1], [2, 3, 1]])
        self.assertIsNone(split["result"])

        mapping = by_name[
            "RejectsMissingParentRowInsteadOfUsingSharedChildGlobalRow"
        ]
        self.assertEqual(mapping["framework"], "GoogleTest")
        self.assertEqual(mapping["input"]["childFaceIdsOfParent"], [[2], [2], []])
        self.assertEqual(mapping["expected"]["exception"], "std::logic_error")
        self.assertEqual(mapping["expected"]["outputNodeIdsUnchanged"], [[[90]]])
        self.assertIsNone(mapping["result"])


    def test_discovers_runtime_parameter_consumer_evidence_tests(self) -> None:
        repository_root = Path(__file__).resolve().parents[2]
        records = discover_tests(repository_root)
        consumer_tests = [
            item for item in records
            if item["source"] == "tools/oneflow_gui/test_parameter_consumer_catalog.py"
        ]
        names = {item["name"] for item in consumer_tests}
        self.assertIn("test_catalog_read_and_behavior_evidence_resolves_to_current_sources", names)
        self.assertIn("test_rk_stage_tracks_all_solver_control_read_sites", names)

    def test_python_test_registry_covers_discovered_test_modules(self) -> None:
        repository_root = Path(__file__).resolve().parents[2]

        audit = audit_python_test_registry(repository_root)

        self.assertTrue(audit["discovered"])
        self.assertTrue(audit["registered"])
        self.assertTrue(audit["excluded"])
        self.assertEqual(audit["unregistered"], ())
        self.assertEqual(audit["stale_exclusions"], ())

    def test_python_test_registry_reports_unregistered_test_modules(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            test_file = root / "tools" / "test_unregistered.py"
            test_file.parent.mkdir(parents=True)
            test_file.write_text(
                "import unittest\n"
                "class UnregisteredTest(unittest.TestCase):\n"
                "    def test_missing_registration(self):\n"
                "        self.assertTrue(True)\n",
                encoding="utf-8",
            )

            audit = audit_python_test_registry(root, commands=(), exclusions={})

            self.assertEqual(
                audit["unregistered"],
                ("tools/test_unregistered.py",),
            )


if __name__ == "__main__":
    unittest.main()
