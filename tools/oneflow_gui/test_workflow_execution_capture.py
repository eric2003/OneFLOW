"""Tests for durable evidence emitted by real workflow process runs."""

from __future__ import annotations

import hashlib
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from tools.oneflow_gui.execution_capture import capture_workflow_execution
from tools.oneflow_gui.workflow_dialog import WorkflowRunWorker
from tools.verification.execution_records import ExecutionRecordStore


class WorkflowExecutionCaptureTest(unittest.TestCase):
    def _git_result(self, arguments, **_kwargs):
        values = {
            ("remote", "get-url", "origin"): "https://github.com/eric2003/OneFLOW.git",
            ("rev-parse", "HEAD"): "abc123",
            ("branch", "--show-current"): "hexin",
        }
        value = values.get(tuple(arguments[1:]), "")
        return SimpleNamespace(stdout=value, returncode=0 if value else 1)

    def test_captures_launched_workflow_step_with_contract_and_artifact_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "runs" / "mesh" / "results" / "mesh.ofl"
            output.parent.mkdir(parents=True)
            output.write_bytes(b"mesh-data")
            log = root / "runs" / "workflow.log"
            log.write_text("mesh completed\\n", encoding="utf-8")
            records_path = root / "verification" / "records.jsonl"
            data = {
                "workflow_id": "mesh-to-solver",
                "run_id": "workflow-run-42",
                "started_at": "2026-11-02T01:00:00+08:00",
                "finished_at": "2026-11-02T01:00:03+08:00",
                "executable": str(root / "build" / "Debug" / "OneFLOW.exe"),
                "log_path": str(log),
                "steps": [
                    {
                        "type": "case", "case_id": "mesh", "status": "succeeded",
                        "command": ["OneFLOW.exe", "r", str(root / "runs" / "mesh")],
                        "return_code": 0, "started_at": "2026-11-02T01:00:00+08:00",
                        "finished_at": "2026-11-02T01:00:02+08:00",
                        "source_path": str(root / "source" / "mesh"),
                        "working_directory": str(root / "build" / "Debug"),
                        "artifacts": [str(output)],
                    },
                    {"type": "case", "case_id": "solver", "status": "skipped"},
                ],
            }
            with patch(
                "tools.oneflow_gui.execution_capture.subprocess.run",
                side_effect=self._git_result,
            ):
                records = capture_workflow_execution(
                    data, repository_root=root, records_path=records_path,
                    contract_ids=("workflow/mesh-to-solver",),
                )
            self.assertEqual(len(records), 1)
            record = records[0]
            self.assertEqual(record.outcome, "passed")
            self.assertEqual(record.framework, "OneFLOW Workflow")
            self.assertEqual(record.test_id, "workflow:mesh-to-solver:mesh")
            self.assertEqual(record.invocation_id, "workflow-run-42")
            self.assertEqual(record.contract_ids, ("workflow/mesh-to-solver",))
            self.assertEqual(record.commit, "abc123")
            self.assertEqual(record.branch, "hexin")
            self.assertEqual(record.repository, "eric2003/OneFLOW")
            self.assertEqual(record.configuration, "Debug")
            self.assertEqual(record.command, tuple(data["steps"][0]["command"]))
            self.assertEqual(
                {item.role for item in record.artifacts},
                {"workflow-log", "workflow-output"},
            )
            artifact = next(item for item in record.artifacts if item.role == "workflow-output")
            self.assertEqual(artifact.sha256, hashlib.sha256(output.read_bytes()).hexdigest())
            self.assertEqual(ExecutionRecordStore(records_path).list_records(), records)

    def test_does_not_record_dependency_blocked_steps_as_process_runs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / "records.jsonl"
            data = {
                "workflow_id": "blocked-workflow",
                "started_at": "2026-11-02T01:00:00+08:00",
                "finished_at": "2026-11-02T01:00:01+08:00",
                "steps": [{"type": "case", "case_id": "solver", "status": "skipped"}],
            }
            self.assertEqual(
                capture_workflow_execution(data, repository_root=root, records_path=path),
                (),
            )
            self.assertFalse(path.exists())

    def test_records_output_validation_failure_even_when_process_exit_code_is_zero(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            data = {
                "workflow_id": "missing-output",
                "started_at": "2026-11-02T01:00:00+08:00",
                "finished_at": "2026-11-02T01:00:01+08:00",
                "steps": [{
                    "type": "case", "case_id": "mesh", "status": "failed",
                    "message": "Declared output artifacts were not produced.",
                    "command": ["OneFLOW.exe", "r", "mesh"], "return_code": 0,
                }],
            }
            record = capture_workflow_execution(
                data, repository_root=root, records_path=root / "records.jsonl"
            )[0]
            self.assertEqual(record.outcome, "failed")
            self.assertEqual(record.exit_code, 0)

    def test_worker_preserves_source_and_staged_execution_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            base_dir = root / "cases"
            (base_dir / "mesh").mkdir(parents=True)
            (base_dir / "solver" / "script").mkdir(parents=True)
            (base_dir / "solver" / "script" / "cfd.txt").write_text(
                'string gridFileName = "grid/old.ofl";\\n', encoding="utf-8"
            )
            executable = root / "build" / "Debug" / "OneFLOW.exe"
            executable.parent.mkdir(parents=True)
            executable.write_bytes(b"")
            worker = WorkflowRunWorker(
                executable, [], "mesh-to-solver", "run-1",
                case_graph=[
                    {"id": "mesh", "directory": "mesh",
                     "outputs": [{"id": "grid", "path": "results/mesh.ofl"}]},
                    {"id": "solver", "directory": "solver", "depends_on": ["mesh"],
                     "inputs": [{"artifact": "grid", "setting": "gridFileName",
                                 "parameter_file": "script/cfd.txt"}]},
                ],
                base_directory=base_dir, run_root=root / "runs" / "run-1",
            )
            payloads = []
            worker.completed.connect(payloads.append)
            def run_process(command, **_kwargs):
                case_dir = Path(command[-1])
                if case_dir.name == "mesh":
                    output = case_dir / "results" / "mesh.ofl"
                    output.parent.mkdir(parents=True, exist_ok=True)
                    output.write_bytes(b"grid")
                return SimpleNamespace(status="succeeded", message="ok", return_code=0)
            with patch("tools.oneflow_gui.workflow_dialog.run_process", side_effect=run_process):
                worker.run()
            self.assertEqual(len(payloads), 1)
            mesh, solver = payloads[0]["steps"]
            self.assertEqual(payloads[0]["status"], "succeeded")
            self.assertEqual(mesh["source_path"], str(base_dir / "mesh"))
            self.assertEqual(mesh["working_directory"], str(executable.parent))
            self.assertEqual(len(mesh["artifacts"]), 1)
            self.assertEqual(solver["source_path"], str(base_dir / "solver"))
            self.assertEqual(solver["command"][1], "r")


if __name__ == "__main__":
    unittest.main()
