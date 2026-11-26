"""Tests for catalog-driven development Case command construction."""

from __future__ import annotations

import unittest
from pathlib import Path, PurePosixPath

from tools.oneflow_gui.case_repository import (
    build_case_command, build_python_case_command, load_cases,
)
from tools.oneflow_gui.models import CaseDefinition


class PythonCaseCommandTest(unittest.TestCase):
    def _case(self, *, requires_gmsh: bool) -> CaseDefinition:
        return CaseDefinition(
            case_id="pipe-bend" if requires_gmsh else "sample-python",
            title="Sample Python Case",
            case_directory=PurePosixPath("/repository/case"),
            runtime_config=None,
            profile=None,
            layout=None,
            execution_kind="python",
            execution_script=PurePosixPath("/repository/run_case.py"),
            requires_gmsh=requires_gmsh,
        )

    def test_required_gmsh_is_passed_to_python_case(self) -> None:
        case = self._case(requires_gmsh=True)
        command = build_python_case_command(
            case,
            PurePosixPath("/workspace/pipe-bend"),
            python_executable="python",
            gmsh_executable="/opt/gmsh/bin/gmsh",
        )
        self.assertEqual(
            command,
            (
                "python",
                "/repository/run_case.py",
                "--case-directory",
                "/workspace/pipe-bend",
                "--gmsh",
                "/opt/gmsh/bin/gmsh",
            ),
        )

    def test_missing_required_gmsh_fails_before_launch(self) -> None:
        with self.assertRaisesRegex(ValueError, "GMSH_EXECUTABLE"):
            build_python_case_command(
                self._case(requires_gmsh=True),
                Path("/workspace/pipe-bend"),
                python_executable="python",
            )

    def test_gmsh_is_not_added_to_python_cases_that_do_not_require_it(self) -> None:
        command = build_python_case_command(
            self._case(requires_gmsh=False),
            PurePosixPath("/workspace/sample"),
            python_executable="python",
        )
        self.assertNotIn("--gmsh", command)

    def test_pipe_bend_native_cpp_implementation_uses_python_runner_adapter(self) -> None:
        repository_root = Path(__file__).resolve().parents[2]
        case = next(
            item for item in load_cases(repository_root)
            if item.case_id == "stl_to_gmsh_pipe_bend_90_degree"
        )
        command = build_case_command(
            case, Path("/workspace/pipe-bend"),
            python_executable="python",
            gmsh_executable="gmsh",
            implementation_id="native_cpp_geometry",
        )
        self.assertEqual(command[0], "python")
        self.assertTrue(command[1].endswith("run_pipe_bend_cpp_case.py"))
        self.assertIn("--case-directory", command)
        self.assertIn("--gmsh", command)
        self.assertNotIn("OneFLOW.exe", command)

    def test_pipe_bend_catalog_declares_gmsh_requirement(self) -> None:
        repository_root = Path(__file__).resolve().parents[2]
        cases = load_cases(repository_root)
        pipe_bend = next(
            case for case in cases if case.case_id == "stl_to_gmsh_pipe_bend_90_degree"
        )
        self.assertTrue(pipe_bend.requires_gmsh)


if __name__ == "__main__":
    unittest.main()
