"""Tests for ordered multi-directory workflow command construction."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.workflow_execution import (
    build_multi_directory_command,
    read_multi_directory_step,
)


class WorkflowExecutionTest(unittest.TestCase):
    def test_command_preserves_directory_order_and_paths_with_spaces(self) -> None:
        with tempfile.TemporaryDirectory(prefix="oneflow workflow ") as directory:
            root = Path(directory)
            executable = root / "bin folder" / "OneFLOW.exe"
            first = root / "grid" / "two connected cubes"
            second = root / "grid" / "minimal cube"
            third = root / "plate script"

            command = build_multi_directory_command(
                executable, [first, second, third]
            )

            self.assertEqual(
                command,
                [
                    str(executable.resolve()),
                    "r",
                    str(first.resolve()),
                    str(second.resolve()),
                    str(third.resolve()),
                ],
            )

    def test_command_rejects_empty_directory_list(self) -> None:
        with self.assertRaisesRegex(ValueError, "at least one directory"):
            build_multi_directory_command(Path("OneFLOW.exe"), [])

    def test_command_rejects_blank_directory_argument(self) -> None:
        with self.assertRaisesRegex(ValueError, "must not be empty"):
            build_multi_directory_command(Path("OneFLOW.exe"), ["grid/cube", " "])

    def test_command_is_an_argument_list_not_a_shell_command(self) -> None:
        command = build_multi_directory_command(
            Path("OneFLOW.exe"), ["folder with spaces", "second-folder"]
        )
        self.assertEqual(command[1], "r")
        self.assertEqual(len(command), 4)
        self.assertEqual(command[2], str(Path("folder with spaces").resolve()))

    def test_read_multi_directory_step_preserves_order(self) -> None:
        steps = (
            {"type": "unsupported", "value": "ignored"},
            {
                "type": "multi_directory",
                "directories": [" grid/first ", "grid/second", "plate"],
            },
        )
        self.assertEqual(
            read_multi_directory_step(steps),
            ["grid/first", "grid/second", "plate"],
        )

    def test_read_multi_directory_step_rejects_missing_step(self) -> None:
        with self.assertRaisesRegex(ValueError, "exactly one"):
            read_multi_directory_step(({"type": "other"},))

    def test_read_multi_directory_step_rejects_duplicate_steps(self) -> None:
        step = {"type": "multi_directory", "directories": ["grid/cube"]}
        with self.assertRaisesRegex(ValueError, "exactly one"):
            read_multi_directory_step((step, step))

    def test_read_multi_directory_step_rejects_empty_directories(self) -> None:
        with self.assertRaisesRegex(ValueError, "at least one directory"):
            read_multi_directory_step(
                ({"type": "multi_directory", "directories": []},)
            )

    def test_read_multi_directory_step_rejects_invalid_directory_values(self) -> None:
        with self.assertRaisesRegex(ValueError, "non-empty string"):
            read_multi_directory_step(
                ({"type": "multi_directory", "directories": ["grid/cube", 7]},)
            )


if __name__ == "__main__":
    unittest.main()
