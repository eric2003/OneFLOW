"""Unit tests for the reusable OneFLOW run-plan primitives."""

from pathlib import Path
import tempfile
import unittest

from oneflow_runner import RunPlan, load_run_plan, save_run_plan


class RunPlanTests(unittest.TestCase):
    def test_command_keeps_executable_and_case_paths_as_separate_arguments(self):
        with tempfile.TemporaryDirectory(prefix="OneFLOW runner ") as temporary:
            root = Path(temporary)
            executable = root / "build output" / "OneFLOW.exe"
            executable.parent.mkdir()
            executable.write_text("", encoding="utf-8")
            case = root / "case with spaces"
            case.mkdir()

            plan = RunPlan(executable, "r", (case,))
            command = plan.command_for(case)

            self.assertEqual(command[0], str(executable.resolve()))
            self.assertEqual(command[1], "r")
            self.assertEqual(command[2], str(case.resolve()))
            self.assertEqual(len(command), 3)

    def test_run_plan_round_trip_preserves_case_order(self):
        with tempfile.TemporaryDirectory(prefix="OneFLOW plan ") as temporary:
            root = Path(temporary)
            executable = root / "OneFLOW.exe"
            executable.write_text("", encoding="utf-8")
            first_case = root / "grid case"
            second_case = root / "solver case"
            first_case.mkdir()
            second_case.mkdir()
            plan_file = root / "workflow.json"

            original = RunPlan(
                executable=executable,
                mode="r",
                case_directories=(first_case, second_case),
                stop_on_error=False,
            )
            save_run_plan(plan_file, original)
            restored = load_run_plan(plan_file)

            self.assertEqual(restored.executable, executable.resolve())
            self.assertEqual(restored.mode, "r")
            self.assertEqual(
                restored.case_directories,
                (first_case.resolve(), second_case.resolve()),
            )
            self.assertFalse(restored.stop_on_error)

    def test_validation_rejects_missing_case_directory(self):
        with tempfile.TemporaryDirectory(prefix="OneFLOW validation ") as temporary:
            root = Path(temporary)
            executable = root / "OneFLOW.exe"
            executable.write_text("", encoding="utf-8")
            plan = RunPlan(executable, "r", (root / "missing-case",))

            with self.assertRaisesRegex(ValueError, "Case directory does not exist"):
                plan.validate()


if __name__ == "__main__":
    unittest.main()
