"""Tests for conservative Contract verification status aggregation."""
from __future__ import annotations

import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

from tools.verification.contract_catalog import ContractEvidence, ContractSummary, EvidenceResolution
from tools.verification.contract_verification import evaluate_contract_verification
from tools.verification.execution_records import ExecutionRecord


class ContractVerificationStatusTest(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name).resolve()
        (self.root / "tests").mkdir()
        (self.root / "tests" / "sample.cpp").write_text("// source\n", encoding="utf-8")
        self.summary = ContractSummary(
            id="mesh/sample", source_id="sample", title="Sample", domain="mesh",
            kind="test", overview="", principle="", architecture="", status="implemented",
            references=(), implementation_paths=(), documentation_paths=(),
            evidence=(
                ContractEvidence("google_test", "SampleSuite.Works", "linked_not_executed", "tests/sample.cpp"),
                ContractEvidence("google_test", "SampleSuite.Rejects", "linked_not_executed", "tests/sample.cpp"),
            ),
            claims=(), limitations=(),
        )
        self.resolutions = tuple(
            EvidenceResolution(item, "resolved", str(self.root / "tests" / "sample.cpp"))
            for item in self.summary.evidence
        )

    def record(self, name, outcome="passed", *, platform="Linux", configuration="Debug", commit="current"):
        return ExecutionRecord(
            test_id="ctest:" + name, framework="CTest", test_name=name,
            source_path="tests/sample.cpp", command=("ctest", "--output-on-failure"),
            outcome=outcome, started_at="2026-11-20T10:00:00+00:00",
            finished_at="2026-11-20T10:00:01+00:00",
            exit_code=0 if outcome == "passed" else 1 if outcome == "failed" else None,
            commit=commit, platform=platform, configuration=configuration,
            contract_ids=(self.summary.id,),
        )

    def evaluate(self, records, resolutions=None, current="current"):
        return evaluate_contract_verification(
            self.summary, self.resolutions if resolutions is None else resolutions,
            records, repository_root=self.root, current_commit=current,
        )

    def test_all_exact_evidence_passed_in_context(self):
        result = self.evaluate((self.record("SampleSuite.Works"), self.record("SampleSuite.Rejects")))
        self.assertEqual(result.status, "passed")
        self.assertEqual(result.contexts[0].status, "passed")
        self.assertEqual([item.status for item in result.evidence], ["passed", "passed"])

    def test_missing_evidence_is_partially_verified_not_passed(self):
        result = self.evaluate((self.record("SampleSuite.Works"),))
        self.assertEqual(result.status, "partially_verified")
        self.assertEqual(result.contexts[0].missing, ("SampleSuite.Rejects",))

    def test_failure_is_not_hidden_by_another_pass(self):
        result = self.evaluate((
            self.record("SampleSuite.Works"),
            self.record("SampleSuite.Rejects", "failed"),
        ))
        self.assertEqual(result.status, "failed")

    def test_blocked_execution_does_not_count_as_pass(self):
        result = self.evaluate((
            self.record("SampleSuite.Works"),
            self.record("SampleSuite.Rejects", "blocked"),
        ))
        self.assertEqual(result.status, "blocked")

    def test_platform_and_configuration_contexts_are_not_merged(self):
        result = self.evaluate((
            self.record("SampleSuite.Works", platform="Linux", configuration="Debug"),
            self.record("SampleSuite.Rejects", platform="Linux", configuration="Debug"),
            self.record("SampleSuite.Works", platform="Windows", configuration="Release"),
        ))
        self.assertEqual(result.status, "partially_verified")
        self.assertEqual(
            [(item.platform, item.configuration, item.status) for item in result.contexts],
            [("Linux", "Debug", "passed"), ("Windows", "Release", "partially_verified")],
        )

    def test_stale_commit_does_not_count_as_current_pass(self):
        result = self.evaluate((
            self.record("SampleSuite.Works", commit="old"),
            self.record("SampleSuite.Rejects", commit="old"),
        ))
        self.assertEqual(result.status, "stale")
        self.assertTrue(all(item.status == "stale" for item in result.evidence))

    def test_contract_id_and_exact_test_identity_are_required(self):
        wrong_identity = self.record("OtherSuite.Works")
        wrong_contract = replace(self.record("SampleSuite.Works"), contract_ids=())
        result = self.evaluate((wrong_identity, wrong_contract))
        self.assertEqual(result.status, "not_run")

    def test_stale_declaration_is_reported(self):
        stale = tuple(
            EvidenceResolution(item, "stale", diagnostic="missing test source")
            if item.reference.endswith("Rejects") else resolution
            for item, resolution in zip(self.summary.evidence, self.resolutions)
        )
        result = self.evaluate((self.record("SampleSuite.Works"),), stale)
        self.assertEqual(result.status, "stale")
        self.assertEqual(result.evidence[1].status, "stale")

    def test_no_test_evidence_is_not_run(self):
        summary = ContractSummary(
            id="mesh/no-tests", source_id="none", title="None", domain="mesh",
            kind="test", overview="", principle="", architecture="", status="proposed",
            references=(), implementation_paths=(), documentation_paths=(),
            evidence=(ContractEvidence("fixture", "mesh/fixture.json", "fixture_available", ""),),
            claims=(), limitations=(),
        )
        result = evaluate_contract_verification(
            summary, (), (), repository_root=self.root, current_commit="current",
        )
        self.assertEqual(result.status, "not_run")
        self.assertIn("No executable", result.diagnostic)


if __name__ == "__main__":
    unittest.main()
