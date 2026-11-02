"""Tests for durable verification execution records."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from tools.verification.execution_records import (
    ExecutionArtifact,
    ExecutionRecord,
    ExecutionRecordStore,
)


def make_record(**overrides) -> ExecutionRecord:
    values = {
        "test_id": "unittest:tools.verification.test_sample.SampleTest.test_rule",
        "framework": "unittest",
        "test_name": "SampleTest.test_rule",
        "source_path": "tools/verification/test_sample.py",
        "command": ("python", "-m", "unittest", "tools.verification.test_sample"),
        "outcome": "passed",
        "started_at": "2026-11-02T01:00:00+08:00",
        "finished_at": "2026-11-02T01:00:01+08:00",
        "exit_code": 0,
        "repository": "eric2003/OneFLOW",
        "commit": "abc123",
        "branch": "hexin",
        "configuration": "Debug",
        "platform": "Windows",
        "runtime_version": "Python 3.11",
        "source_line": 42,
        "working_directory": "D:/github/OneFLOW",
        "contract_ids": ("runtime/gui-process-execution",),
        "artifacts": (
            ExecutionArtifact("log", ".oneflow/logs/test.log"),
        ),
        "diagnostics": (),
    }
    values.update(overrides)
    return ExecutionRecord(**values)


class ExecutionRecordTest(unittest.TestCase):
    def test_round_trip_preserves_execution_identity_and_evidence(self) -> None:
        record = make_record()

        restored = ExecutionRecord.from_dict(record.to_dict())

        self.assertEqual(restored, record)
        self.assertEqual(restored.contract_ids, ("runtime/gui-process-execution",))
        self.assertEqual(restored.artifacts[0].role, "log")

    def test_invocation_id_round_trips_and_legacy_records_default_empty(self) -> None:
        record = make_record(invocation_id="workflow-run-42")
        restored = ExecutionRecord.from_dict(record.to_dict())
        self.assertEqual(restored.invocation_id, "workflow-run-42")

        legacy_data = make_record().to_dict()
        legacy_data.pop("invocation_id")
        legacy = ExecutionRecord.from_dict(legacy_data)
        self.assertEqual(legacy.invocation_id, "")

    def test_missing_optional_metadata_is_preserved_as_unknown(self) -> None:
        record = make_record(commit="", branch="", platform="", runtime_version="")

        self.assertEqual(record.commit, "")
        self.assertEqual(record.freshness("abc123"), "unknown")
        self.assertEqual(record.freshness(None), "unknown")

    def test_freshness_requires_exact_revision_match(self) -> None:
        record = make_record()

        self.assertEqual(record.freshness("abc123"), "current")
        self.assertEqual(record.freshness("def456"), "stale")

    def test_rejects_finished_time_before_start(self) -> None:
        with self.assertRaisesRegex(ValueError, "finished_at"):
            make_record(finished_at="2026-11-02T00:59:59+08:00")

    def test_rejects_timezone_free_timestamps(self) -> None:
        with self.assertRaisesRegex(ValueError, "timezone"):
            make_record(started_at="2026-11-02T01:00:00")

    def test_rejects_unknown_outcome(self) -> None:
        with self.assertRaisesRegex(ValueError, "outcome"):
            make_record(outcome="green")

    def test_passed_outcome_requires_zero_exit_code(self) -> None:
        with self.assertRaisesRegex(ValueError, "exit_code 0"):
            make_record(exit_code=1)

    def test_rejects_invalid_source_line_and_artifact_checksum(self) -> None:
        with self.assertRaisesRegex(ValueError, "source_line"):
            make_record(source_line=0)
        with self.assertRaisesRegex(ValueError, "sha256"):
            ExecutionArtifact("log", "test.log", "not-a-checksum")

    def test_jsonl_store_appends_and_reads_multiple_records(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "verification" / "execution-records.jsonl"
            store = ExecutionRecordStore(path)
            first = make_record()
            second = make_record(test_id="gtest:Suite.Case", framework="GoogleTest")

            store.append(first)
            store.append(second)

            self.assertEqual(store.list_records(), (first, second))
            self.assertEqual(len(path.read_text(encoding="utf-8").splitlines()), 2)

    def test_missing_store_is_an_empty_history(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            store = ExecutionRecordStore(Path(temporary) / "missing.jsonl")

            self.assertEqual(store.list_records(), ())

    def test_corrupt_line_reports_its_location_without_silent_data_loss(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "records.jsonl"
            valid_line = json.dumps(make_record().to_dict())
            path.write_text(valid_line + "\nnot-json\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "line 2"):
                ExecutionRecordStore(path).list_records()

    def test_rejects_unsupported_schema_version(self) -> None:
        data = make_record().to_dict()
        data["schema_version"] = 2

        with self.assertRaisesRegex(ValueError, "schema"):
            ExecutionRecord.from_dict(data)


if __name__ == "__main__":
    unittest.main()
