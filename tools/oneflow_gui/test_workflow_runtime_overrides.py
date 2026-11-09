"""Tests for runtime parameter overrides."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.workflow_runtime_overrides import apply_runtime_overrides


class WorkflowRuntimeOverridesTest(unittest.TestCase):
    def test_updates_existing_setting_and_preserves_comment(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            parameter_file = Path(temporary) / "cfd.txt"
            parameter_file.write_text(
                'string gridFileName = "grid/old.ofl"; // selected grid\n'
                'int maxSteps = 50;\n',
                encoding="utf-8",
            )

            changed = apply_runtime_overrides(
                parameter_file,
                {"gridFileName": "results/rae2822_2d.ofl"},
            )

            self.assertEqual(changed, ("gridFileName",))
            content = parameter_file.read_text(encoding="utf-8")
            self.assertIn(
                'string gridFileName = "results/rae2822_2d.ofl"; // selected grid',
                content,
            )
            self.assertIn("int maxSteps = 50;", content)

    def test_injects_source_case_directory_when_setting_is_absent(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            parameter_file = Path(temporary) / "cfd.txt"
            parameter_file.write_text('string gridFileName = "grid/a.ofl";\n', encoding="utf-8")

            apply_runtime_overrides(
                parameter_file,
                {
                    "sourceGridCaseDir": "D:/runs/mesh producer",
                    "gridFileName": "results/a.ofl",
                },
            )

            content = parameter_file.read_text(encoding="utf-8")
            self.assertIn(
                'string sourceGridCaseDir = "D:/runs/mesh producer";',
                content,
            )
            self.assertIn('string gridFileName = "results/a.ofl";', content)

    def test_escapes_quotes_and_backslashes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            parameter_file = Path(temporary) / "case.txt"
            parameter_file.write_text("", encoding="utf-8")

            apply_runtime_overrides(parameter_file, {"pathSetting": 'a"b\\c'})

            self.assertIn(
                'string pathSetting = "a\\"b\\\\c";',
                parameter_file.read_text(encoding="utf-8"),
            )

    def test_preserves_gbk_parameter_file_encoding_and_comments(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            parameter_file = Path(temporary) / "cfd.txt"
            comment = bytes.fromhex("cdf8b8f1c2b7beb6").decode("gbk")
            original = f'string gridFileName = "grid/old.ofl"; // {comment}\n'
            parameter_file.write_bytes(original.encode("gbk"))

            apply_runtime_overrides(
                parameter_file,
                {"gridFileName": "results/rae2822_2d.ofl"},
            )

            content = parameter_file.read_bytes().decode("gbk")
            self.assertIn(
                f'string gridFileName = "results/rae2822_2d.ofl"; // {comment}',
                content,
            )
            self.assertIn(comment, content)

    def test_numeric_overrides_keep_int_and_real_parameter_types(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "cfd.txt"
            path.write_text("int maxSteps = 50; // steps\nReal global_dt = 0.001;\n", encoding="utf-8")
            apply_runtime_overrides(path, {"maxSteps": 2, "global_dt": 0.002})
            content = path.read_text(encoding="utf-8")
            self.assertIn("int maxSteps = 2; // steps", content)
            self.assertIn("Real global_dt = 0.002;", content)
            self.assertNotIn('"2"', content)

    def test_rejects_wrong_numeric_override_type_without_writing(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "cfd.txt"
            original = "int maxSteps = 50;\n"
            path.write_text(original, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "integer values"):
                apply_runtime_overrides(path, {"maxSteps": "2"})
            self.assertEqual(path.read_text(encoding="utf-8"), original)

    def test_rejects_invalid_setting_name_without_writing(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            parameter_file = Path(temporary) / "case.txt"
            parameter_file.write_text("int steps = 2;\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "Invalid runtime setting"):
                apply_runtime_overrides(parameter_file, {"bad name": "x"})

            self.assertEqual(parameter_file.read_text(encoding="utf-8"), "int steps = 2;\n")


if __name__ == "__main__":
    unittest.main()
