"""Standard-library tests for desktop services; these do not require PySide6."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.case_repository import load_cases, select_cases_by_execution_kind
from tools.oneflow_gui.executable_discovery import executable_candidates
from tools.oneflow_gui.models import CaseDefinition
from tools.oneflow_gui.runner import (
    _copy_case_to_workspace,
    build_command,
    prepare_case_database_workspace,
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
        self.assertTrue(all(
            case.runtime_config is None or case.runtime_config.is_file()
            for case in cases if case.execution_kind == "oneflow"
        ))

    def test_gmsh_cases_are_available_in_the_gmsh_development_group(self) -> None:
        cases = load_cases(REPOSITORY_ROOT)
        cases_by_id = {case.case_id: case for case in cases}
        self.assertIn("two_quads", cases_by_id)
        self.assertEqual(cases_by_id["two_quads"].group, "Gmsh")
        self.assertEqual(cases_by_id["gmsh_to_cgns"].group, "Gmsh")
        self.assertEqual(cases_by_id["gmsh_to_oneflow"].group, "Gmsh")

    def test_stl_catalog_contains_native_and_python_development_cases(self) -> None:
        cases = [case for case in load_cases(REPOSITORY_ROOT) if case.group == "STL"]
        cases_by_id = {case.case_id: case for case in cases}

        self.assertEqual(len(cases), 15)
        self.assertEqual(cases_by_id["stl_ascii_basic"].execution_kind, "oneflow")
        self.assertIsNone(cases_by_id["stl_ascii_basic"].runtime_config)
        self.assertIsNone(cases_by_id["stl_ascii_basic"].execution_script)
        self.assertEqual(cases_by_id["stl_ascii_crlf"].execution_kind, "oneflow")
        self.assertEqual(cases_by_id["stl_ascii_scientific_notation"].execution_kind, "oneflow")
        self.assertEqual(cases_by_id["stl_binary_basic"].execution_kind, "oneflow")
        self.assertEqual(cases_by_id["stl_binary_header_solid"].execution_kind, "oneflow")
        self.assertEqual(cases_by_id["stl_non_planar_surface"].execution_kind, "oneflow")
        self.assertEqual(cases_by_id["stl_weld_tolerance"].execution_kind, "oneflow")
        self.assertEqual(cases_by_id["stl_ascii_scientific_notation"].workspace_relative_directory, Path("stl/ascii_scientific_notation"))
        self.assertEqual(cases_by_id["stl_python_ascii_basic"].execution_kind, "python")
        self.assertEqual(cases_by_id["stl_python_ascii_crlf"].execution_kind, "python")
        self.assertEqual(cases_by_id["stl_ascii_crlf"].workspace_relative_directory, Path("stl/ascii_crlf"))
        self.assertEqual(cases_by_id["stl_python_ascii_crlf"].workspace_relative_directory, Path("stl/python/ascii_crlf"))
        self.assertIsNone(cases_by_id["stl_ascii_crlf"].runtime_config)
        self.assertIsNone(cases_by_id["stl_ascii_crlf"].execution_script)

        python_cases = [case for case in cases if case.execution_kind == "python"]
        self.assertEqual(len(python_cases), 8)
        self.assertTrue(all(case.runtime_config is None for case in python_cases))
        self.assertTrue(all(case.execution_script is not None and case.execution_script.is_file()
                            for case in python_cases))
        self.assertTrue(all(case.case_directory.parent.name == "python" for case in python_cases))

    def test_development_case_execution_mode_separates_native_and_python_cases(self) -> None:
        cases = [case for case in load_cases(REPOSITORY_ROOT) if case.group == "STL"]
        native = select_cases_by_execution_kind(cases, "oneflow")
        python = select_cases_by_execution_kind(cases, "python")
        self.assertEqual({case.case_id for case in native}, {
            "stl_ascii_basic", "stl_ascii_crlf", "stl_ascii_scientific_notation",
            "stl_binary_basic", "stl_binary_header_solid",
            "stl_non_planar_surface", "stl_weld_tolerance"
        })
        self.assertEqual(len(python), 8)
        self.assertEqual(len(native), 7)
        self.assertTrue(all(case.execution_kind == "oneflow" for case in native))
        self.assertTrue(all(case.execution_kind == "python" for case in python))
        with self.assertRaisesRegex(ValueError, "Unsupported"):
            select_cases_by_execution_kind(cases, "mixed")

    def test_python_stl_development_cases_run_from_staged_workspaces(self) -> None:
        cases = [
            case for case in load_cases(REPOSITORY_ROOT)
            if case.group == "STL" and case.execution_kind == "python"
        ]
        with tempfile.TemporaryDirectory() as directory:
            staged_cases = prepare_case_database_workspace(cases, Path(directory) / "grid")
            for case in cases:
                staged = staged_cases[case.case_id]
                result = subprocess.run(
                    [sys.executable, str(case.execution_script), "--case-directory", str(staged)],
                    cwd=REPOSITORY_ROOT, check=True, capture_output=True, text=True,
                )
                report = json.loads(
                    (staged / "results" / "stl_report.json").read_text(encoding="utf-8")
                )
                self.assertTrue(report["passed"], result.stdout)
                self.assertFalse((case.case_directory / "results").exists())

    def test_native_stl_cases_are_staged_without_runtime_manifests(self) -> None:
        cases = [
            case for case in load_cases(REPOSITORY_ROOT)
            if case.case_id in {
                "stl_ascii_basic", "stl_ascii_crlf", "stl_ascii_scientific_notation",
                "stl_binary_basic", "stl_binary_header_solid",
            }
        ]
        self.assertEqual(len(cases), 5)
        with tempfile.TemporaryDirectory() as directory:
            staged_cases = prepare_case_database_workspace(cases, Path(directory) / "grid")
            for case in cases:
                staged = staged_cases[case.case_id]
                self.assertFalse((staged / "script" / "control.txt").exists())
                self.assertTrue((staged / "script" / "dimension.txt").is_file())
                self.assertTrue((staged / "script" / "simu.txt").is_file())
                self.assertTrue((staged / "script" / "grid.txt").is_file())
                self.assertTrue((staged / "input" / "model.stl").is_file())
                self.assertFalse((case.case_directory / "results").exists())
            crlf_fixture = (staged_cases["stl_ascii_crlf"] / "input" / "model.stl").read_bytes()
            self.assertIn(bytes((13, 10)), crlf_fixture)
            scientific_fixture = (staged_cases["stl_ascii_scientific_notation"] / "input" / "model.stl").read_text(encoding="ascii")
            self.assertIn("e+", scientific_fixture.lower())
            for case_id in ("stl_binary_basic", "stl_binary_header_solid"):
                binary_fixture = (staged_cases[case_id] / "input" / "model.stl").read_bytes()
                self.assertEqual(len(binary_fixture), 184)
                self.assertEqual(int.from_bytes(binary_fixture[80:84], "little"), 2)
            solid_header = (staged_cases["stl_binary_header_solid"] / "input" / "model.stl").read_bytes()
            self.assertTrue(solid_header.startswith(b"solid"))

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
                platform_name="nt",
            )
            self.assertEqual(candidates[:2], [debug.resolve(), release.resolve()])
            self.assertIn(installed.resolve(), candidates)
            self.assertLess(
                candidates.index(release.resolve()),
                candidates.index(installed.resolve()),
            )

    def test_debug_and_release_choices_exist_before_building(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            candidates = executable_candidates(root, {}, platform_name="nt")
            self.assertEqual(
                candidates[:2],
                [
                    (root / "build" / "bin" / "Debug" / "OneFLOW.exe").resolve(),
                    (root / "build" / "bin" / "Release" / "OneFLOW.exe").resolve(),
                ],
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

    def test_case_database_mode_refreshes_the_workspace(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source_case"
            source.mkdir()
            (source / "control.txt").write_text("fresh", encoding="utf-8")
            workspace = root / "grid-layout-test"
            destination = workspace / "sample"
            destination.mkdir(parents=True)
            (destination / "stale-result.dat").write_text("stale", encoding="utf-8")
            case = CaseDefinition(
                case_id="sample",
                title="Sample",
                case_directory=source,
                runtime_config=source / "control.txt",
                profile=None,
                layout=None,
            )

            actual = prepare_case_directory(case, workspace, "ignored-run-id", "case_database")

            self.assertTrue(actual.samefile(destination))
            self.assertEqual((actual / "control.txt").read_text(encoding="utf-8"), "fresh")
            self.assertFalse((actual / "stale-result.dat").exists())
            self.assertTrue((source / "control.txt").is_file())

    def test_case_database_refreshes_entire_grid_and_stages_selected_cases(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source_case"
            source.mkdir()
            (source / "control.txt").write_text("fresh", encoding="utf-8")
            workspace = root / "mesh-layout-test" / "grid"
            stale = workspace / "old_case"
            stale.mkdir(parents=True)
            (stale / "stale.dat").write_text("stale", encoding="utf-8")
            destination = workspace / "sample"
            destination.mkdir()
            (destination / "stale-result.dat").write_text("stale", encoding="utf-8")
            case = CaseDefinition(
                case_id="sample",
                title="Sample",
                case_directory=source,
                runtime_config=source / "control.txt",
                profile=None,
                layout=None,
            )

            staged = prepare_case_database_workspace([case], workspace)

            self.assertTrue(staged["sample"].samefile(destination))
            self.assertEqual((destination / "control.txt").read_text(encoding="utf-8"), "fresh")
            self.assertFalse((destination / "stale-result.dat").exists())
            self.assertFalse(stale.exists())
            self.assertTrue((source / "control.txt").is_file())

    def test_case_database_preserves_catalog_directory_hierarchy(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "repository" / "grid" / "stl" / "ascii_crlf"
            source.mkdir(parents=True)
            (source / "model.stl").write_text("fixture", encoding="ascii")
            workspace = root / "mesh-layout-test" / "grid"
            case = CaseDefinition(
                case_id="stl_ascii_crlf",
                title="STL ASCII CRLF",
                case_directory=source,
                runtime_config=None,
                profile=None,
                layout=None,
                workspace_relative_directory=Path("stl/ascii_crlf"),
            )

            staged = prepare_case_database_workspace([case], workspace)
            expected = workspace / "stl" / "ascii_crlf"
            self.assertTrue(staged["stl_ascii_crlf"].samefile(expected))
            self.assertEqual((expected / "model.stl").read_text(encoding="ascii"), "fixture")
            self.assertFalse((workspace / "stl_ascii_crlf").exists())

            prepared = prepare_case_directory(
                case, workspace, "unused-run-id", "case_database_prepared"
            )
            self.assertTrue(prepared.samefile(expected))

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
