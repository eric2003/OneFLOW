from __future__ import annotations

import json
import unittest
from pathlib import Path

from .contract_catalog import (
    ContractSummary,
    adapt_geometry_scenario,
    adapt_theory_entry,
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

        self.assertEqual(summary.id, "sst_f2_eddy")
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

        self.assertEqual(summary.id, "registered_child_resolution")
        self.assertEqual(summary.domain, "geometry")
        self.assertEqual(summary.kind, "scenario")
        self.assertEqual(summary.evidence[0].kind, "fixture")
        self.assertEqual(summary.evidence[0].status, "fixture_available")
        linked_tests = [
            item for item in summary.evidence if item.kind == "linked_test"
        ]
        self.assertEqual(len(linked_tests), 1)
        self.assertEqual(linked_tests[0].status, "linked_not_executed")
        self.assertTrue(
            any("does not execute production C++" in x for x in summary.limitations)
        )

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
        self.assertEqual(issues, [], "\\n".join(issues))

    def test_reports_missing_and_escaping_paths(self) -> None:
        summary = ContractSummary(
            id="path_check",
            title="Path check",
            domain="test",
            kind="test",
            overview="",
            principle="",
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
        self.assertTrue(any("escapes root" in issue for issue in issues))

    def test_reports_duplicate_contract_ids(self) -> None:
        summary = ContractSummary(
            id="duplicate",
            title="Duplicate",
            domain="test",
            kind="test",
            overview="",
            principle="",
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

        self.assertEqual(issues, ["duplicate contract id: duplicate"])


if __name__ == "__main__":
    unittest.main()
