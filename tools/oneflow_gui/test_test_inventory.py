"""Tests for source-based test inventory discovery."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.test_inventory import discover_tests


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


if __name__ == "__main__":
    unittest.main()
