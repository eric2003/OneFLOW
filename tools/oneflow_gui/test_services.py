"""Standard-library tests for desktop services; these do not require PySide6."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.case_repository import load_cases
from tools.oneflow_gui.executable_discovery import executable_candidates
from tools.oneflow_gui.models import CaseDefinition
from tools.oneflow_gui.runner import (
    _copy_case_to_workspace,
    build_command,
    prepare_case_directory,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


class DesktopServiceTest(unittest.TestCase):
    def test_registered_cases_load_from_shared_catalog(self) -> None:
        cases = load_cases(REPOSITORY_ROOT)
        case_ids = {case.case_id for case in cases}
        self.assertTrue(
            {"minimal_rectangle", "minimal_cube", "two_connected_cubes"} <= case_ids
        )
        self.assertTrue(all(case.runtime_config.is_file() for case in cases))

    def test_executable_discovery_prefers_debug_then_release_then_install(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            debug = root / "build" / "bin" / "Debug" / "OneFLOW.exe"
            release = root / "build" / "bin" / "Release" / "OneFLOW.exe"
            install_root = root / "ProgramFiles"
            installed = install_root / "OneFLOW" / "bin" / "OneFLOW.exe"
            for path in (debug, release, installed):
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b"test")

            candidates = executable_candidates(
                root,
                {"ProgramFiles": str(install_root)},
            )
            self.assertEqual(candidates[:2], [debug.resolve(), release.resolve()])
            self.assertIn(installed.resolve(), candidates)
            self.assertLess(
                candidates.index(release.resolve()),
                candidates.index(installed.resolve()),
            )

    def test_run_command_uses_argument_list(self) -> None:
        with tempfile.TemporaryDirectory(prefix="oneflow path ") as directory:
            executable = Path(directory) / "OneFLOW.exe"
            case_directory = Path(directory) / "case with spaces"
            command = build_command(executable, case_directory)
            self.assertEqual(command, [str(executable.resolve()), "r", str(case_directory.resolve())])

    def test_workspace_copy_keeps_source_case_unchanged(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source_case"
            (source / "script").mkdir(parents=True)
            source_config = source / "script" / "control.txt"
            source_config.write_text('"simu.txt";\n', encoding="utf-8")
            (source / "script" / "simu.txt").write_text(
                'string simutask = "Grid";\n', encoding="utf-8"
            )
            case = CaseDefinition(
                case_id="sample",
                title="Sample",
                case_directory=source,
                runtime_config=source_config,
                profile=None,
                layout=None,
            )

            copied = _copy_case_to_workspace(case, root / "runs", "run-001")
            self.assertTrue((copied / "script" / "control.txt").is_file())
            self.assertEqual(source_config.read_text(encoding="utf-8"), '"simu.txt";\n')
            self.assertNotEqual(copied, source)

    def test_source_execution_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source_case"
            source.mkdir()
            case = CaseDefinition(
                case_id="sample",
                title="Sample",
                case_directory=source,
                runtime_config=source / "control.txt",
                profile=None,
                layout=None,
            )

            with self.assertRaisesRegex(ValueError, "In-place execution is disabled"):
                prepare_case_directory(
                    case, root / "runs", "run-001", "source"
                )
            self.assertFalse((root / "runs").exists())

    def test_custom_execution_copies_case_under_selected_directory(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source_case"
            source.mkdir()
            (source / "control.txt").write_text("config", encoding="utf-8")
            case = CaseDefinition(
                case_id="sample",
                title="Sample",
                case_directory=source,
                runtime_config=source / "control.txt",
                profile=None,
                layout=None,
            )

            actual = prepare_case_directory(
                case, root / "custom_runs", "run-002", "custom_workspace"
            )
            expected = root / "custom_runs" / "run-002" / "sample"
            self.assertTrue(actual.samefile(expected))
            self.assertTrue((actual / "control.txt").is_file())
            self.assertTrue((source / "control.txt").is_file())

    def test_workspace_copy_refuses_to_overwrite_existing_run(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source_case"
            source.mkdir()
            case = CaseDefinition(
                case_id="sample",
                title="Sample",
                case_directory=source,
                runtime_config=source / "control.txt",
                profile=None,
                layout=None,
            )
            destination = root / "runs" / "run-001" / "sample"
            destination.mkdir(parents=True)
            with self.assertRaises(FileExistsError):
                _copy_case_to_workspace(case, root / "runs", "run-001")


if __name__ == "__main__":
    unittest.main()
