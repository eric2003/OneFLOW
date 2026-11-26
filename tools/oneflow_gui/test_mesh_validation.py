"""Mesh-to-CFD validation Case catalog and graph tests."""

from __future__ import annotations

import hashlib
import json
import sys
import tempfile
import unittest
from unittest import mock
from types import SimpleNamespace
from pathlib import Path

from tools.oneflow_gui.mesh_validation import (
    build_mesh_validation_cases, cfd_template_choices, mesh_case_choices, mesh_case_groups,
    mesh_case_lifecycle_summary, mesh_case_requires_oneflow_executable,
    mesh_case_group_categories,
    mesh_validation_workspace, project_mesh_validation_workspace,
    reset_mesh_validation_workspace, create_mesh_validation_run_workspace,
    prepare_mesh_validation_run_workspace,
    _copy_existing_mesh, _graph_step_records,
    _mesh_case_command, _mesh_validation_stage_results, run_report_stage_rows,
    existing_mesh_case_paths, _validate_required_physical_groups,
    _geometry_lineage, geometry_lineage_rows, _project_required_physical_groups,
    _geometry_artifact_lineage_path, _write_geometry_artifact_lineage,
    _read_geometry_artifact_lineage, _downstream_lineage_artifacts,
    _collect_run_artifacts, _artifact_file_state, _sha256_file,
    _persist_graph_stage_outcomes,
    recover_mesh_validation_run, list_mesh_validation_runs,
)
from tools.oneflow_gui.case_repository import _case_dimension
from tools.oneflow_gui.workflow_graph_runner import (
    CaseExecutionResult, WorkflowGraphRunResult,
)


class RunArtifactInventoryTest(unittest.TestCase):
    def test_inventory_records_relative_paths_sizes_hashes_and_categories(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            files = {
                "stl_to_gmsh/input/model.stl": b"geometry",
                "stl_to_gmsh/results/model.msh": b"mesh",
                "gmsh_to_cgns/results/model.cgns": b"cgns",
                "oneflow_cfd_solve/results/residual.dat": b"residual",
                "oneflow_cfd_solve/results/solver.log": b"log",
                "gmsh_to_cgns/results/mesh-validation.json": b"{}",
                "misc/custom.bin": b"other",
            }
            for relative, payload in files.items():
                path = root / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(payload)
            report_path = root / "mesh-validation-report.json"
            report_path.write_text("{}", encoding="utf-8")
            (root / ".temporary.tmp").write_text("temporary", encoding="utf-8")

            records = _collect_run_artifacts(root, report_path)
            by_path = {item["path"]: item for item in records}

            self.assertEqual(set(by_path), set(files))
            self.assertEqual(by_path["stl_to_gmsh/input/model.stl"]["category"], "Geometry")
            self.assertEqual(by_path["stl_to_gmsh/results/model.msh"]["category"], "Volume Mesh")
            self.assertEqual(by_path["gmsh_to_cgns/results/model.cgns"]["category"], "Volume Mesh")
            self.assertEqual(by_path["oneflow_cfd_solve/results/residual.dat"]["category"], "CFD Results")
            self.assertEqual(by_path["oneflow_cfd_solve/results/solver.log"]["category"], "Logs")
            self.assertEqual(by_path["gmsh_to_cgns/results/mesh-validation.json"]["category"], "Reports")
            self.assertEqual(by_path["misc/custom.bin"]["category"], "Other")
            self.assertEqual(by_path["stl_to_gmsh/input/model.stl"]["stage"], "stl_to_gmsh")
            self.assertEqual(by_path["stl_to_gmsh/input/model.stl"]["size_bytes"], len(b"geometry"))
            self.assertEqual(
                by_path["stl_to_gmsh/input/model.stl"]["sha256"],
                hashlib.sha256(b"geometry").hexdigest(),
            )
            self.assertNotIn(str(root), by_path["stl_to_gmsh/input/model.stl"]["path"])
            self.assertNotIn("mesh-validation-report.json", by_path)

    def test_artifact_paths_report_available_missing_and_reject_escape_paths(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            present = root / "results/mesh.msh"
            present.parent.mkdir()
            present.write_bytes(b"mesh")
            self.assertEqual(
                _artifact_file_state("results/mesh.msh", root)[0], "Available"
            )
            self.assertEqual(
                _artifact_file_state("results/missing.cgns", root)[0], "Missing"
            )
            self.assertEqual(
                _artifact_file_state("../../outside.txt", root),
                ("Unsafe path", None),
            )
            self.assertEqual(
                _artifact_file_state("results/mesh.msh", None),
                ("Unresolved", None),
            )

    def test_sha256_streaming_matches_standard_digest_for_large_file(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "large.mesh"
            payload = b"OneFLOW mesh evidence" * 100_000
            path.write_bytes(payload)
            self.assertEqual(_sha256_file(path), hashlib.sha256(payload).hexdigest())

    def test_inventory_skips_symlinks_and_missing_run_directories(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "run"
            self.assertEqual(_collect_run_artifacts(root), [])
            root.mkdir()
            target = root / "target.dat"
            target.write_text("data", encoding="utf-8")
            link = root / "link.dat"
            try:
                link.symlink_to(target)
            except (OSError, NotImplementedError):
                self.skipTest("Symlinks are unavailable in this environment.")
            self.assertEqual(
                [item["path"] for item in _collect_run_artifacts(root)],
                ["target.dat"],
            )


class RecoverMeshValidationRunTest(unittest.TestCase):
    def test_graph_final_outcomes_override_preliminary_stage_log_status(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            log = root / "stl_to_gmsh/results/mesh-validation.log"
            log.parent.mkdir(parents=True)
            log.write_text(
                '$ ["gmsh", "model.geo"]\n'
                "[stl_to_gmsh] succeeded: Process returned zero.\n",
                encoding="utf-8",
            )
            plan = SimpleNamespace(steps=[
                SimpleNamespace(
                    case_id="stl_to_gmsh",
                    runtime_directory=Path("stl_to_gmsh"),
                ),
                SimpleNamespace(
                    case_id="gmsh_to_cgns",
                    runtime_directory=Path("gmsh_to_cgns"),
                ),
            ])
            result = SimpleNamespace(cases=[
                CaseExecutionResult(
                    "stl_to_gmsh", "failed",
                    "Declared output artifacts were not produced: mesh.cgns",
                ),
                CaseExecutionResult(
                    "gmsh_to_cgns", "skipped", "Blocked by unsuccessful dependencies.",
                ),
            ])

            _persist_graph_stage_outcomes(result, plan, root)
            report = recover_mesh_validation_run(root)

            self.assertEqual(report["stages"]["mesh_generation"]["status"], "FAIL")
            self.assertIn("Declared output artifacts", report["stages"]["mesh_generation"]["message"])
            self.assertEqual(report["stages"]["gmsh_to_cgns"]["status"], "NOT_RUN")
            self.assertEqual(report["steps"][1]["status"], "skipped")

    def test_stage_log_serialization_persists_terminal_outcome(self):
        from tools.oneflow_gui.mesh_validation import _format_mesh_validation_stage_log

        rendered = _format_mesh_validation_stage_log(
            "gmsh_to_cgns", ["OneFLOW.exe", "r", "case"],
            ["Converting mesh...", "Conversion complete."],
            "succeeded", "CGNS file created.",
        )
        self.assertIn('$ ["OneFLOW.exe", "r", "case"]\n', rendered)
        self.assertTrue(rendered.endswith(
            "[gmsh_to_cgns] succeeded: CGNS file created.\n"
        ))

    def test_recovers_terminal_stage_states_metrics_and_verified_lineage(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            producer_log = root / "stl_to_gmsh/results/mesh-validation.log"
            producer_log.parent.mkdir(parents=True)
            producer_log.write_text(
                '$ ["gmsh", "model.geo"]\n'
                "[stl_to_gmsh] succeeded: Mesh generated and validated.\n",
                encoding="utf-8",
            )
            cgns_log = root / "gmsh_to_cgns/results/mesh-validation.log"
            cgns_log.parent.mkdir(parents=True)
            cgns_log.write_text(
                '$ ["OneFLOW.exe", "r", "gmsh_to_cgns"]\n'
                "[gmsh_to_cgns] failed: Conversion rejected an invalid element.\n",
                encoding="utf-8",
            )
            mesh_report = root / "stl_to_gmsh/results/mesh-validation.json"
            mesh_report.write_text(
                json.dumps({
                    "nodes": 120,
                    "elements": 240,
                    "tetrahedron_elements": 180,
                    "tetrahedron_quality": {
                        "minimum_quality": 0.25,
                        "mean_quality": 0.72,
                        "minimum_volume": 0.001,
                        "total_volume": 1.5,
                        "quality_percentiles": {
                            "p05": 0.31, "p50": 0.74, "p95": 0.93,
                        },
                        "quality_below_threshold_counts": {
                            "0.1": 0, "0.2": 2, "0.3": 8,
                        },
                    },
                    "physical_groups": [{"name": "fluid"}],
                }),
                encoding="utf-8",
            )
            cgns = root / "gmsh_to_cgns/results/mesh.cgns"
            cgns.write_bytes(b"cgns")
            lineage = {
                "provenance": {"kind": "generated", "notes": "Verified source."},
                "source_representation": {"format": "STL", "sha256": "a" * 64},
            }
            _write_geometry_artifact_lineage(cgns, "CGNS volume mesh", lineage)
            (root / "mesh-validation-console.log").write_text(
                "durable console output\n", encoding="utf-8"
            )

            report = recover_mesh_validation_run(root)

            self.assertEqual(report["status"], "recovered")
            self.assertIsNone(report["exit_code"])
            self.assertEqual(report["stages"]["mesh_generation"]["status"], "PASS")
            self.assertEqual(report["stages"]["gmsh_to_cgns"]["status"], "FAIL")
            self.assertEqual(report["stages"]["cgns_to_oneflow"]["status"], "NOT_CHECKED")
            self.assertEqual(report["metrics"], {
                "nodes": 120,
                "elements": 240,
                "tetrahedron_elements": 180,
                "tetrahedron_mean_ratio_minimum": 0.25,
                "tetrahedron_mean_ratio_mean": 0.72,
                "tetrahedron_minimum_volume": 0.001,
                "tetrahedron_total_volume": 1.5,
                "tetrahedron_mean_ratio_p05": 0.31,
                "tetrahedron_mean_ratio_p50": 0.74,
                "tetrahedron_mean_ratio_p95": 0.93,
                "tetrahedron_mean_ratio_below_0_1_count": 0,
                "tetrahedron_mean_ratio_below_0_2_count": 2,
                "tetrahedron_mean_ratio_below_0_3_count": 8,
                "physical_group_count": 1,
            })
            self.assertEqual(report["geometry_lineage"]["provenance"]["kind"], "generated")
            self.assertEqual(report["recovery"]["verified_geometry_sidecars"], 1)
            self.assertTrue(any(
                item["path"] == "mesh-validation-console.log"
                for item in report["artifacts"]["files"]
            ))

    def test_older_stage_logs_recover_terminal_status_from_console_log(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            stage_log = root / "stl_to_gmsh/results/mesh-validation.log"
            stage_log.parent.mkdir(parents=True)
            stage_log.write_text(
                '$ ["gmsh", "model.geo"]\nRaw Gmsh output only.\n',
                encoding="utf-8",
            )
            (root / "mesh-validation-console.log").write_text(
                "[stl_to_gmsh] succeeded: Legacy run completed.\n",
                encoding="utf-8",
            )

            report = recover_mesh_validation_run(root)

            self.assertEqual(report["stages"]["mesh_generation"]["status"], "PASS")
            self.assertEqual(report["steps"][0]["message"], "Legacy run completed.")

    def test_run_metadata_restores_case_parameters_without_inventing_outcomes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata = {
                "schema": "oneflow.mesh-validation-run",
                "version": 1,
                "run_id": "run-123",
                "mesh_case_id": "stl_to_gmsh_sphere",
                "mesh_implementation": "python",
                "mesh_input_mode": "generated_mesh",
                "cfd_template_id": "laminar_flat_plate",
                "max_steps": 2,
                "started_at_utc": "2026-11-25T08:00:00+00:00",
            }
            (root / "mesh-validation-run.json").write_text(
                json.dumps(metadata), encoding="utf-8"
            )

            report = recover_mesh_validation_run(root)

            self.assertEqual(report["run_id"], "run-123")
            self.assertEqual(report["mesh_case_id"], "stl_to_gmsh_sphere")
            self.assertEqual(report["mesh_implementation"], "python")
            self.assertEqual(report["cfd_template_id"], "laminar_flat_plate")
            self.assertEqual(report["max_steps"], 2)
            self.assertEqual(report["recovery"]["run_metadata_status"], "readable")
            self.assertEqual(report["stages"]["mesh_generation"]["status"], "NOT_CHECKED")

    def test_artifact_presence_alone_never_claims_a_stage_pass(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            mesh = root / "stl_to_gmsh/results/model.msh"
            mesh.parent.mkdir(parents=True)
            mesh.write_bytes(b"mesh")

            report = recover_mesh_validation_run(root)

            self.assertEqual(report["status"], "recovered")
            self.assertEqual(report["stages"]["mesh_generation"]["status"], "NOT_CHECKED")
            self.assertEqual(report["stages"]["gmsh_to_cgns"]["status"], "NOT_CHECKED")
            self.assertEqual(len(report["steps"]), 4)
            self.assertTrue(all(item["status"] == "unknown" for item in report["steps"]))

    def test_recovers_existing_mesh_input_stage_name(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            log = root / "stl_to_gmsh/results/mesh-validation.log"
            log.parent.mkdir(parents=True)
            log.write_text(
                '$ ["<copy-existing-mesh>", "source.msh"]\n'
                "[stl_to_gmsh] succeeded: Existing mesh staged.\n",
                encoding="utf-8",
            )

            report = recover_mesh_validation_run(root)

            self.assertNotIn("mesh_generation", report["stages"])
            self.assertEqual(report["stages"]["existing_msh_input"]["status"], "PASS")


class GeometryArtifactLineageTest(unittest.TestCase):
    def test_lineage_is_persisted_and_verified_for_cgns_and_native_grid(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cgns = root / "mesh.cgns"
            cgns.write_bytes(b"cgns-bytes")
            lineage = {
                "provenance": {"kind": "unknown", "notes": "Historical source unknown."},
                "source_representation": {"format": "STL", "sha256": "a" * 64},
                "derived_mesh": {"format": "Gmsh MSH", "sha256": "b" * 64},
            }
            cgns_sidecar = _write_geometry_artifact_lineage(
                cgns, "CGNS volume mesh", lineage,
                {"role": "Gmsh MSH", "sha256": "b" * 64},
            )
            cgns_record = _read_geometry_artifact_lineage(cgns)
            self.assertEqual(cgns_record["geometry_lineage"], lineage)
            self.assertEqual(cgns_record["parent_artifact"]["role"], "Gmsh MSH")
            self.assertTrue(cgns_sidecar.is_file())
            self.assertTrue(cgns_sidecar.read_bytes().endswith(b"\n"))
            self.assertNotIn(b"\\n", cgns_sidecar.read_bytes())

            native = root / "mesh.oneflow.grid"
            native.write_bytes(b"native-grid")
            native_sidecar = _write_geometry_artifact_lineage(
                native, "OneFLOW native grid", cgns_record["geometry_lineage"],
                {
                    "role": "CGNS volume mesh",
                    "file_name": cgns.name,
                    "sha256": cgns_record["artifact"]["sha256"],
                },
            )
            native_record = _read_geometry_artifact_lineage(native)
            self.assertEqual(native_record["parent_artifact"]["file_name"], cgns.name)
            self.assertEqual(native_record["geometry_lineage"], lineage)
            self.assertTrue(native_sidecar.is_file())

    def test_lineage_sidecar_is_rejected_after_artifact_changes(self):
        with tempfile.TemporaryDirectory() as directory:
            artifact = Path(directory) / "mesh.cgns"
            artifact.write_bytes(b"original")
            _write_geometry_artifact_lineage(artifact, "CGNS volume mesh", {})
            artifact.write_bytes(b"replaced")
            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                _read_geometry_artifact_lineage(artifact)


class GeometryLineageTest(unittest.TestCase):
    def test_unknown_historical_source_is_preserved_and_mesh_hash_is_recorded(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            stl = root / "model.stl"
            msh = root / "model.msh"
            stl.write_bytes(b"historical-stl")
            msh.write_bytes(b"derived-mesh")
            entry = {
                "geometry_provenance": {
                    "kind": "unknown",
                    "notes": "Original STL source and generator are not documented.",
                }
            }
            lineage = _geometry_lineage(entry, stl, msh, root)
            self.assertEqual(lineage["provenance"]["kind"], "unknown")
            self.assertNotIn("asset_id", lineage)
            self.assertEqual(lineage["source_representation"]["sha256"],
                             hashlib.sha256(b"historical-stl").hexdigest())
            self.assertEqual(lineage["derived_mesh"]["sha256"],
                             hashlib.sha256(b"derived-mesh").hexdigest())
            downstream = root / "model.cgns"
            downstream.write_bytes(b"cgns")
            lineage = _geometry_lineage(
                entry, stl, msh, root,
                [{"role": "CGNS volume mesh", "format": "CGNS", "path": downstream}],
            )
            rows = geometry_lineage_rows({"geometry_lineage": lineage})
            self.assertEqual(rows[0][1], "UNKNOWN")
            self.assertEqual(rows[-2][0], "Derived volume mesh")
            self.assertEqual(rows[-1][0], "CGNS volume mesh")

    def test_matching_manifest_provides_asset_and_generator_lineage(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            stl = root / "generated.stl"
            stl.write_bytes(b"generated-stl")
            digest = hashlib.sha256(stl.read_bytes()).hexdigest()
            manifest = {
                "schema": "oneflow.stl-geometry-manifest",
                "asset_id": "geometry.sphere.unit",
                "stl_sha256": digest,
                "definition_sha256": "a" * 64,
                "generator": {"name": "oneflow_stl_geometry", "version": "1.0",
                              "source_sha256": "b" * 64},
            }
            stl.with_suffix(".geometry.json").write_text(
                json.dumps(manifest), encoding="utf-8"
            )
            lineage = _geometry_lineage({}, stl, None, root)
            self.assertEqual(lineage["asset_id"], "geometry.sphere.unit")
            self.assertEqual(lineage["provenance"]["kind"], "generated")
            self.assertEqual(lineage["generator"]["source_sha256"], "b" * 64)

    def test_geometry_lineage_rows_accepts_string_generator_metadata(self):
        rows = geometry_lineage_rows({
            "geometry_lineage": {"generator": "oneflow_stl_geometry"}
        })
        self.assertIn(
            ("Geometry generator", "RECORDED", "oneflow_stl_geometry"),
            rows,
        )

    def test_geometry_lineage_rows_formats_structured_generator_metadata(self):
        rows = geometry_lineage_rows({
            "geometry_lineage": {
                "generator": {
                    "name": "oneflow_stl_geometry",
                    "version": "1.0",
                    "source_sha256": "b" * 64,
                }
            }
        })
        self.assertIn(
            ("Geometry generator", "RECORDED",
             f"oneflow_stl_geometry 1.0 | source SHA-256: {'b' * 64}"),
            rows,
        )

    def test_manifest_with_mismatched_stl_checksum_is_not_trusted(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            stl = root / "tampered.stl"
            stl.write_bytes(b"changed")
            stl.with_suffix(".geometry.json").write_text(json.dumps({
                "schema": "oneflow.stl-geometry-manifest",
                "asset_id": "geometry.sphere.unit",
                "stl_sha256": "0" * 64,
            }), encoding="utf-8")
            lineage = _geometry_lineage({}, stl, None, root)
            self.assertNotIn("asset_id", lineage)
            self.assertEqual(lineage["manifest_status"], "invalid_or_checksum_mismatch")


class MeshValidationCatalogTest(unittest.TestCase):
    def test_project_boundary_regions_become_required_msh_physical_groups(self):
        self.assertEqual(
            _project_required_physical_groups({
                "boundary_regions": ["wall", "inflow", "outflow"]
            }),
            [
                {"name": "wall", "dimension": 2},
                {"name": "inflow", "dimension": 2},
                {"name": "outflow", "dimension": 2},
                {"name": "fluid_domain", "dimension": 3},
            ],
        )
        self.assertEqual(_project_required_physical_groups({}), [])
        with self.assertRaisesRegex(ValueError, "unique names"):
            _project_required_physical_groups({"boundary_regions": ["wall", "wall"]})
        with self.assertRaisesRegex(ValueError, "unsupported"):
            _project_required_physical_groups({"boundary_regions": ["farfield"]})

    def test_repository_catalog_routes_generic_3d_mesh_through_cfd_smoke(self):
        repository_root = Path(__file__).resolve().parents[2]
        case_id = "stl_to_gmsh_generic_3d_mesh_handoff"
        catalog = json.loads(
            (repository_root / "grid/case_catalog.json").read_text(encoding="utf-8")
        )
        entry = next(
            item for item in catalog["mesh_validation_inputs"] if item["id"] == case_id
        )
        self.assertEqual(entry["mesh_validation"], "cfd_smoke")
        self.assertEqual(entry["dimension"], 3)
        steps = build_mesh_validation_cases(
            repository_root, case_id, "laminar_flat_plate"
        )
        self.assertEqual(
            [step["id"] for step in steps],
            ["stl_to_gmsh", "gmsh_to_cgns", "cgns_to_oneflow", "oneflow_cfd_solve"],
        )
        self.assertEqual(steps[-1]["runtime_directory"], "oneflow_cfd_solve")
        self.assertEqual(steps[-1]["directory"], "test/cfd/plateuns2dslau2")


    def test_historical_cube_cases_explicitly_declare_unknown_geometry_origin(self):
        repository_root = Path(__file__).resolve().parents[2]
        catalog = json.loads(
            (repository_root / "grid/case_catalog.json").read_text(encoding="utf-8")
        )
        cube_ids = {
            "stl_to_gmsh_cube_single_hex",
            "stl_to_gmsh_python_cube_single_hex",
            "stl_to_gmsh_cube_six_tetrahedra",
            "stl_to_gmsh_python_cube_six_tetrahedra",
            "stl_to_gmsh_cube_automatic_tetrahedra",
            "stl_to_gmsh_python_cube_automatic_tetrahedra",
        }
        entries = {item["id"]: item for item in catalog["cases"]}
        self.assertTrue(cube_ids.issubset(entries))
        for case_id in cube_ids:
            provenance = entries[case_id]["geometry_provenance"]
            self.assertEqual(provenance["kind"], "unknown")
            self.assertIn("unknown", provenance["notes"].lower())
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        (self.root / "grid/stl_to_gmsh/2d/square/single_quad/cpp/script").mkdir(parents=True)
        (self.root / "grid/stl_to_gmsh/2d/square/two_triangles/cpp/script").mkdir(parents=True)
        (self.root / "grid/stl_to_gmsh/2d/square/two_triangles/python").mkdir(parents=True)
        (self.root / "grid/stl_to_gmsh/2d/square/automatic_triangles/cpp/script").mkdir(parents=True)
        (self.root / "grid/stl_to_gmsh/2d/square/automatic_triangles/python").mkdir(parents=True)
        (self.root / "grid/stl_to_gmsh/2d/circle_external_flow/cpp/script").mkdir(parents=True)
        (self.root / "grid/stl_to_gmsh/3d/cube/single_hex/cpp/script").mkdir(parents=True)
        (self.root / "grid/stl_to_gmsh/3d/sphere/automatic_tetrahedra/cpp/script").mkdir(parents=True)
        (self.root / "grid/stl_to_gmsh/3d/sphere/automatic_tetrahedra/python/script").mkdir(parents=True)
        (self.root / "grid/stl_to_gmsh/3d/cylinder_internal_flow").mkdir(parents=True)
        (self.root / "grid/stl_to_gmsh/3d/cylinder_internal_flow/oblique_cut/cpp/script").mkdir(parents=True)
        (self.root / "grid/stl_to_gmsh/3d/cylinder_internal_flow/oblique_cut/python").mkdir(parents=True)
        (self.root / "grid/stl_to_gmsh/2d/square/single_quad/python").mkdir(parents=True)
        (self.root / "test/cfd/plateuns2dslau2/script").mkdir(parents=True)
        (self.root / "grid/case_catalog.json").write_text(json.dumps({
            "schema": "oneflow.case-catalog",
            "cases": [
                {"id": "stl_to_gmsh_square", "title": "One Square QUAD_4 Cell",
                 "directory": "grid/stl_to_gmsh/2d/square/single_quad/cpp", "group": "STL to Gmsh", "dimension": 2},
                {"id": "stl_to_gmsh_square_two_triangles",
                 "title": "Square with Two TRI_3 Cells (Native C++)",
                 "directory": "grid/stl_to_gmsh/2d/square/two_triangles/cpp",
                 "group": "STL to Gmsh", "dimension": 2},
                {"id": "stl_to_gmsh_square_automatic_triangles",
                 "title": "Square with Automatic TRI_3 Refinement (Native C++)",
                 "directory": "grid/stl_to_gmsh/2d/square/automatic_triangles/cpp",
                 "group": "STL to Gmsh", "dimension": 2},
                {"id": "stl_to_gmsh_circle_external_flow",
                 "title": "Circle External Flow (Native C++)",
                 "directory": "grid/stl_to_gmsh/2d/circle_external_flow/cpp",
                 "group": "STL to Gmsh", "dimension": 2},
                {"id": "native_cube", "title": "Cube Single Hex (Native C++)",
                 "directory": "grid/stl_to_gmsh/3d/cube/single_hex/cpp",
                 "group": "STL to Gmsh", "dimension": 3},
                {"id": "python_cube", "title": "Cube Single Hex (Python Reference)",
                 "directory": "grid/stl_to_gmsh/3d/cube/single_hex/python",
                 "group": "STL to Gmsh", "dimension": 3,
                 "execution": {"kind": "python", "script": "grid/stl_to_gmsh/run_python_case.py"}},
                {"id": "native_sphere", "title": "Unit Sphere Automatic Tetrahedra (Native C++)",
                 "directory": "grid/stl_to_gmsh/3d/sphere/automatic_tetrahedra/cpp",
                 "group": "STL to Gmsh", "dimension": 3, "mesh_validation": "cfd_smoke"},
                {"id": "python_sphere", "title": "Unit Sphere Automatic Tetrahedra (Python Reference)",
                 "directory": "grid/stl_to_gmsh/3d/sphere/automatic_tetrahedra/python",
                 "group": "STL to Gmsh", "dimension": 3,
                 "execution": {"kind": "python", "script": "grid/stl_to_gmsh/run_python_case.py"},
                 "mesh_validation": "cfd_smoke"},
                {"id": "cylinder_existing_mesh",
                 "title": "3D Cylinder Internal Flow (Generated MSH)",
                 "directory": "grid/stl_to_gmsh/3d/cylinder_internal_flow",
                 "group": "STL to Gmsh", "dimension": 3,
                 "mesh_validation": "cfd_smoke",
                 "execution": {"kind": "existing_mesh",
                               "project_directory": "cylinder_internal_flow",
                               "mesh_file": "results/cylinder_internal_flow.msh"},
                 "required_physical_groups": [
                     {"name": "wall", "dimension": 2},
                     {"name": "inflow", "dimension": 2},
                     {"name": "outflow", "dimension": 2},
                     {"name": "fluid_domain", "dimension": 3},
                 ]},
                {"id": "python_square", "title": "One Square QUAD_4 Cell (Python Reference)",
                 "directory": "grid/stl_to_gmsh/2d/square/single_quad/python",
                 "group": "STL to Gmsh", "dimension": 2,
                 "execution": {"kind": "python", "script": "grid/stl_to_gmsh/run_python_case.py"}},
                {"id": "stl_to_gmsh_python_square_two_triangles",
                 "title": "Square with Two TRI_3 Cells (Python Reference)",
                 "directory": "grid/stl_to_gmsh/2d/square/two_triangles/python",
                 "group": "STL to Gmsh", "dimension": 2,
                 "execution": {"kind": "python", "script": "grid/stl_to_gmsh/run_python_case.py"}},
                {"id": "stl_to_gmsh_python_square_automatic_triangles",
                 "title": "Square with Automatic TRI_3 Refinement (Python Reference)",
                 "directory": "grid/stl_to_gmsh/2d/square/automatic_triangles/python",
                 "group": "STL to Gmsh", "dimension": 2,
                 "execution": {"kind": "python", "script": "grid/stl_to_gmsh/run_python_case.py"}},
            ],
        }), encoding="utf-8")
        catalog_path = self.root / "grid/case_catalog.json"
        catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
        cylinder_input = next(
            item for item in catalog["cases"] if item["id"] == "cylinder_existing_mesh"
        )
        catalog["cases"].remove(cylinder_input)
        catalog["mesh_validation_inputs"] = [
            cylinder_input,
            {"id": "stl_to_gmsh_cylinder_internal_flow_oblique_cut",
             "title": "Internal-flow Cylinder with Generic Oblique Cut (Native C++)",
             "directory": "grid/stl_to_gmsh/3d/cylinder_internal_flow/oblique_cut/cpp",
             "group": "STL to Gmsh", "dimension": 3, "mesh_validation": "cfd_smoke"},
            {"id": "stl_to_gmsh_python_cylinder_internal_flow_oblique_cut",
             "title": "Internal-flow Cylinder with Generic Oblique Cut (Python Reference)",
             "directory": "grid/stl_to_gmsh/3d/cylinder_internal_flow/oblique_cut/python",
             "group": "STL to Gmsh", "dimension": 3, "mesh_validation": "cfd_smoke",
             "execution": {"kind": "python",
                           "script": "grid/stl_to_gmsh/3d/cylinder_internal_flow/oblique_cut/python/run_cylinder_oblique_cut.py",
                           "output_directory_argument": "--output-directory"}},
        ]
        catalog_path.write_text(json.dumps(catalog), encoding="utf-8")
        (self.root / "grid/stl_to_gmsh/2d/square/single_quad/python/expected.json").write_text(
            json.dumps({"artifact": {"mesh_file": "square.msh"}}), encoding="utf-8"
        )
        (self.root / "grid/stl_to_gmsh/2d/circle_external_flow/cpp/script/grid.txt").write_text(
            'string targetGridFileName = "results/circle_external_flow.msh";\n',
            encoding="utf-8")
        (self.root / "grid/stl_to_gmsh/3d/cube/single_hex/cpp/script/grid.txt").write_text(
            'string sourceGridFileName = "input/model.stl";\n'
            'string targetGridFileName = "results/cube_single_hex.msh";\n', encoding="utf-8")
        (self.root / "grid/stl_to_gmsh/3d/cube/single_hex/cpp/script/layout.json").write_text(
            json.dumps({"mesh_dimension": 3}), encoding="utf-8")
        (self.root / "grid/stl_to_gmsh/3d/sphere/automatic_tetrahedra/cpp/script/grid.txt").write_text(
            'string targetGridFileName = "results/sphere_automatic_tetrahedra.msh";\n',
            encoding="utf-8")
        (self.root / "grid/stl_to_gmsh/3d/sphere/automatic_tetrahedra/cpp/script/layout.json").write_text(
            json.dumps({"mesh_dimension": 3}), encoding="utf-8")
        (self.root / "grid/stl_to_gmsh/3d/sphere/automatic_tetrahedra/python/expected.json").write_text(
            json.dumps({"artifact": {"mesh_file": "sphere_automatic_tetrahedra.msh"}}),
            encoding="utf-8")
        (self.root / "grid/stl_to_gmsh/3d/cylinder_internal_flow/expected.json").write_text(
            json.dumps({"artifact": {"mesh_file": "cylinder_internal_flow.msh"}}),
            encoding="utf-8")
        (self.root / "grid/stl_to_gmsh/3d/cylinder_internal_flow/oblique_cut/cpp/script/grid.txt").write_text(
            'string targetGridFileName = "results/cylinder_oblique_cut.msh";\n', encoding="utf-8")
        (self.root / "grid/stl_to_gmsh/3d/cylinder_internal_flow/oblique_cut/python/expected.json").write_text(
            json.dumps({"artifact": {"mesh_file": "cylinder_oblique_cut.msh"}}), encoding="utf-8")
        (self.root / "grid/stl_to_gmsh/2d/square/single_quad/cpp/script/grid.txt").write_text(
            'string targetGridFileName = "results/square.msh";\n', encoding="utf-8")
        (self.root / "grid/stl_to_gmsh/2d/square/two_triangles/cpp/script/grid.txt").write_text(
            'string targetGridFileName = "results/square_two_triangles.msh";\n', encoding="utf-8")
        (self.root / "grid/stl_to_gmsh/2d/square/two_triangles/python/expected.json").write_text(
            json.dumps({"artifact": {"mesh_file": "square_two_triangles.msh"}}),
            encoding="utf-8")
        (self.root / "grid/stl_to_gmsh/2d/square/automatic_triangles/cpp/script/grid.txt").write_text(
            'string targetGridFileName = "results/square_automatic_triangles.msh";\n',
            encoding="utf-8")
        (self.root / "grid/stl_to_gmsh/2d/square/automatic_triangles/python/expected.json").write_text(
            json.dumps({"artifact": {"mesh_file": "square_automatic_triangles.msh"}}),
            encoding="utf-8")
        (self.root / "test/cfd/plateuns2dslau2/script/cfd.txt").write_text(
            'string resFile = "results/res.dat";\n'
            'string flowFile = "restart/flow.dat";\n'
            'string visualFile = "visual/visual.dat";\n', encoding="utf-8")

    def tearDown(self):
        self.temporary.cleanup()

    def test_case_dimension_uses_explicit_metadata_then_layout_or_directory(self):
        two_d = {"directory": "grid/mesh-layout/2d/minimal"}
        self.assertEqual(_case_dimension(two_d, self.root), 2)

        layout_case = self.root / "grid/mesh-layout/2d/custom_volume_case"
        (layout_case / "script").mkdir(parents=True)
        (layout_case / "script/layout.json").write_text(
            json.dumps({"mesh_dimension": 3}), encoding="utf-8"
        )
        self.assertEqual(
            _case_dimension({"directory": "grid/mesh-layout/2d/custom_volume_case"}, layout_case), 3
        )
        self.assertEqual(
            _case_dimension(
                {"directory": "grid/mesh-layout/2d/custom_volume_case", "dimension": 2}, layout_case
            ),
            2,
        )
        self.assertIsNone(_case_dimension({"directory": "grid/unknown_case"}, self.root))

    def test_mesh_selector_contains_native_and_python_implementations(self):
        self.assertEqual(mesh_case_choices(self.root), (
            ("stl_to_gmsh_square", "One Square QUAD_4 Cell"),
            ("stl_to_gmsh_square_two_triangles", "Square with Two TRI_3 Cells (Native C++)"),
            ("stl_to_gmsh_square_automatic_triangles", "Square with Automatic TRI_3 Refinement (Native C++)"),
            ("stl_to_gmsh_circle_external_flow", "Circle External Flow (Native C++)"),
            ("native_cube", "Cube Single Hex (Native C++)"),
            ("python_cube", "Cube Single Hex (Python Reference)"),
            ("native_sphere", "Unit Sphere Automatic Tetrahedra (Native C++)"),
            ("python_sphere", "Unit Sphere Automatic Tetrahedra (Python Reference)"),
            ("python_square", "One Square QUAD_4 Cell (Python Reference)"),
            ("stl_to_gmsh_python_square_two_triangles", "Square with Two TRI_3 Cells (Python Reference)"),
            ("stl_to_gmsh_python_square_automatic_triangles", "Square with Automatic TRI_3 Refinement (Python Reference)"),
            ("cylinder_existing_mesh", "3D Cylinder Internal Flow (Generated MSH)"),
            ("stl_to_gmsh_cylinder_internal_flow_oblique_cut", "Internal-flow Cylinder with Generic Oblique Cut (Native C++)"),
            ("stl_to_gmsh_python_cylinder_internal_flow_oblique_cut", "Internal-flow Cylinder with Generic Oblique Cut (Python Reference)"),
        ))
    def test_existing_cylinder_is_grouped_as_existing_gmsh_mesh_not_stl_producer(self):
        groups = dict(mesh_case_group_categories(self.root))
        categories = dict(groups["3D"])
        self.assertIn("STL to Gmsh", categories)
        self.assertEqual(
            categories["Existing Gmsh Mesh"],
            (("cylinder_existing_mesh", "3D Cylinder Internal Flow (Generated MSH)"),),
        )

    def test_repository_catalog_exposes_3d_cylinder_existing_mesh(self):
        repository_root = Path(__file__).resolve().parents[2]
        case_id = "stl_to_gmsh_existing_mesh_3d_cylinder_internal_flow"
        choices = dict(mesh_case_choices(repository_root))
        self.assertIn(case_id, choices)

        groups = dict(mesh_case_groups(repository_root))
        self.assertIn(case_id, {item[0] for item in groups["3D"]})

        catalog = json.loads(
            (repository_root / "grid/case_catalog.json").read_text(encoding="utf-8")
        )
        entry = next(
            item for item in catalog["mesh_validation_inputs"] if item["id"] == case_id
        )
        self.assertEqual(entry["mesh_validation"], "cfd_smoke")
        steps = build_mesh_validation_cases(
            repository_root, case_id, "laminar_flat_plate"
        )
        self.assertEqual(
            [step["id"] for step in steps],
            ["stl_to_gmsh", "gmsh_to_cgns", "cgns_to_oneflow", "oneflow_cfd_solve"],
        )
        self.assertEqual(steps[3]["runtime_directory"], "oneflow_cfd_solve")
        self.assertEqual(steps[3]["directory"], "test/cfd/plateuns2dslau2")
        self.assertEqual(
            steps[0]["outputs"][0]["path"],
            "results/cylinder_internal_flow.msh",
        )

    def test_repository_catalog_exposes_unit_sphere_in_mesh_validation_3d(self):
        repository_root = Path(__file__).resolve().parents[2]
        choices = dict(mesh_case_choices(repository_root))
        cfd_smoke_ids = {
            "stl_to_gmsh_sphere_automatic_tetrahedra",
            "stl_to_gmsh_python_sphere_automatic_tetrahedra",
        }
        internal_flow_cfd_smoke_ids = {
            "stl_to_gmsh_cylinder_internal_flow_automatic_tetrahedra",
            "stl_to_gmsh_python_cylinder_internal_flow_automatic_tetrahedra",
        }
        cfd_smoke_ids |= internal_flow_cfd_smoke_ids
        expected_ids = cfd_smoke_ids
        self.assertTrue(expected_ids.issubset(choices))

        groups = dict(mesh_case_groups(repository_root))
        sphere_ids = {case_id for case_id, _title in groups["3D"]}
        self.assertTrue(expected_ids.issubset(sphere_ids))

        for case_id in sorted(cfd_smoke_ids):
            steps = build_mesh_validation_cases(
                repository_root, case_id, "laminar_flat_plate"
            )
            self.assertEqual(
                [step["id"] for step in steps],
                ["stl_to_gmsh", "gmsh_to_cgns", "cgns_to_oneflow", "oneflow_cfd_solve"],
            )
            self.assertEqual(steps[3]["runtime_directory"], "oneflow_cfd_solve")
            self.assertEqual(steps[3]["directory"], "test/cfd/plateuns2dslau2")


    def test_mesh_selector_groups_cases_by_spatial_dimension(self):
        groups = mesh_case_groups(self.root)
        self.assertEqual([dimension for dimension, _cases in groups], ["2D", "3D"])
        self.assertEqual(
            [case_id for case_id, _title in groups[0][1]],
            [
                "stl_to_gmsh_square",
                "stl_to_gmsh_square_two_triangles",
                "stl_to_gmsh_square_automatic_triangles",
                "stl_to_gmsh_circle_external_flow",
                "python_square",
                "stl_to_gmsh_python_square_two_triangles",
                "stl_to_gmsh_python_square_automatic_triangles",
            ],
        )
        self.assertEqual(
            [case_id for case_id, _title in groups[1][1]],
            [
                "native_cube", "python_cube", "native_sphere", "python_sphere",
                "cylinder_existing_mesh",
                "stl_to_gmsh_cylinder_internal_flow_oblique_cut",
                "stl_to_gmsh_python_cylinder_internal_flow_oblique_cut",
            ],
        )

    def test_mesh_producer_command_dispatches_by_implementation(self):
        python_entry = {
            "execution": {
                "kind": "python",
                "script": "grid/stl_to_gmsh/run_python_case.py",
            }
        }
        python_output_entry = {
            "execution": {
                "kind": "python",
                "script": "grid/stl_to_gmsh/3d/cylinder_internal_flow/oblique_cut/python/run_cylinder_oblique_cut.py",
                "output_directory_argument": "--output-directory",
            }
        }
        command = _mesh_case_command(
            self.root, python_entry, Path("OneFLOW"), Path("gmsh"),
            self.root / ".oneflow/staged/python_case",
        )
        self.assertEqual(command[0], sys.executable)
        self.assertEqual(command[1], str((self.root / "grid/stl_to_gmsh/run_python_case.py").resolve()))
        self.assertEqual(command[2:], [
            "--case-directory", str(self.root / ".oneflow/staged/python_case"),
            "--gmsh", "gmsh",
        ])
        staged_output = self.root / ".oneflow/staged/oblique_case"
        oblique_command = _mesh_case_command(
            self.root, python_output_entry, Path("OneFLOW"), Path("gmsh"),
            staged_output,
        )
        self.assertEqual(
            oblique_command[-2:],
            ["--output-directory", str(staged_output)],
        )
        native = _mesh_case_command(
            self.root, {}, Path("OneFLOW"), Path("gmsh"),
            self.root / ".oneflow/staged/native_case",
        )
        self.assertEqual(native, [
            "OneFLOW", "r", str(self.root / ".oneflow/staged/native_case")
        ])

    def test_python_mesh_producer_script_cannot_escape_repository(self):
        with self.assertRaisesRegex(ValueError, "inside the repository"):
            _mesh_case_command(
                self.root,
                {"execution": {"kind": "python", "script": "../outside.py"}},
                Path("OneFLOW"), Path("gmsh"), self.root / "staged",
            )

    def test_circle_external_flow_is_available_to_mesh_to_cfd_validation(self):
        cases = build_mesh_validation_cases(
            self.root, "stl_to_gmsh_circle_external_flow", "laminar_flat_plate"
        )
        self.assertEqual(len(cases), 4)
        self.assertEqual(
            cases[0]["directory"], "grid/stl_to_gmsh/2d/circle_external_flow/cpp"
        )
        self.assertEqual(
            cases[0]["outputs"][0]["path"],
            "results/circle_external_flow.msh",
        )
        self.assertEqual(cases[1]["id"], "gmsh_to_cgns")
        self.assertEqual(
            cases[1]["outputs"][0]["path"],
            "results/circle_external_flow.cgns",
        )
        self.assertEqual(cases[2]["id"], "cgns_to_oneflow")
        self.assertEqual(
            cases[2]["outputs"][0]["path"],
            "results/circle_external_flow.oneflow.grid",
        )
        self.assertEqual(cases[2]["outputs"][1], {
            "id": "bc_region_map", "path": "results/circle_external_flow_bcRegionMap.txt",
            "parameter_file": "script/grid.txt", "setting": "bcRegionMapFile",
        })
        self.assertEqual(cases[2]["outputs"][2], {
            "id": "grid_info", "path": "results/circle_external_flow.info.json",
            "parameter_file": "script/grid.txt", "setting": "gridInfoFileName",
        })
        self.assertEqual(cases[3]["id"], "oneflow_cfd_solve")
        self.assertEqual(
            cases[3]["inputs"],
            [
                {"artifact": "oneflow_grid", "setting": "gridFileName",
                 "parameter_file": "script/cfd.txt"},
                {"artifact": "bc_region_map", "setting": "bcRegionMapFile",
                 "parameter_file": "script/cfd.txt"},
            ],
        )

    def test_square_two_triangles_cases_are_available_to_mesh_to_cfd_validation(self):
        native = build_mesh_validation_cases(
            self.root, "stl_to_gmsh_square_two_triangles", "laminar_flat_plate"
        )
        python = build_mesh_validation_cases(
            self.root, "stl_to_gmsh_python_square_two_triangles", "laminar_flat_plate"
        )
        self.assertEqual(len(native), 4)
        self.assertEqual(len(python), 4)
        self.assertEqual(
            native[0]["directory"], "grid/stl_to_gmsh/2d/square/two_triangles/cpp"
        )
        self.assertEqual(
            native[0]["outputs"][0]["path"], "results/square_two_triangles.msh"
        )
        self.assertEqual(
            python[0]["directory"], "grid/stl_to_gmsh/2d/square/two_triangles/python"
        )
        self.assertEqual(
            python[0]["outputs"][0]["path"], "results/square_two_triangles.msh"
        )
        self.assertEqual(native[1]["id"], "gmsh_to_cgns")
        self.assertEqual(python[1]["id"], "gmsh_to_cgns")

    def test_automatic_triangle_cases_are_available_to_mesh_to_cfd_validation(self):
        native = build_mesh_validation_cases(
            self.root, "stl_to_gmsh_square_automatic_triangles", "laminar_flat_plate"
        )
        python = build_mesh_validation_cases(
            self.root, "stl_to_gmsh_python_square_automatic_triangles", "laminar_flat_plate"
        )
        self.assertEqual(len(native), 4)
        self.assertEqual(len(python), 4)
        self.assertEqual(native[0]["directory"],
                         "grid/stl_to_gmsh/2d/square/automatic_triangles/cpp")
        self.assertEqual(python[0]["directory"],
                         "grid/stl_to_gmsh/2d/square/automatic_triangles/python")
        self.assertEqual(native[0]["outputs"][0]["path"],
                         "results/square_automatic_triangles.msh")
        self.assertEqual(python[0]["outputs"][0]["path"],
                         "results/square_automatic_triangles.msh")
        self.assertEqual(native[1]["id"], "gmsh_to_cgns")
        self.assertEqual(python[1]["id"], "gmsh_to_cgns")

    def test_python_mesh_case_uses_the_same_artifact_handoff_contract(self):
        cases = build_mesh_validation_cases(
            self.root, "python_square", "laminar_flat_plate"
        )
        self.assertEqual(len(cases), 4)
        self.assertEqual(cases[0]["directory"], "grid/stl_to_gmsh/2d/square/single_quad/python")
        self.assertEqual(cases[0]["outputs"], [
            {"id": "gmsh_mesh", "path": "results/square.msh"}
        ])
        self.assertEqual(cases[0]["validation"], {
            "kind": "gmsh_msh4", "path": "results/square.msh"
        })
        self.assertEqual(cases[1]["inputs"][0]["artifact"], "gmsh_mesh")
        self.assertEqual(cases[1]["inputs"][0]["setting"], "sourceGridFileName")
        self.assertEqual(cases[1]["directory"], "grid/conversion/gmsh_to_cgns")
        self.assertEqual(cases[3]["inputs"][1]["artifact"], "bc_region_map")
        self.assertEqual(cases[3]["inputs"][1]["setting"], "bcRegionMapFile")



    def test_pipe_bend_runs_the_complete_three_dimensional_cfd_smoke_chain(self):
        repository_root = Path(__file__).resolve().parents[2]
        case_id = "stl_to_gmsh_pipe_bend_90_degree"
        cases = build_mesh_validation_cases(repository_root, case_id, "laminar_flat_plate")
        expected_steps = [
            "stl_to_gmsh", "gmsh_to_cgns", "cgns_to_oneflow",
            "oneflow_cfd_solve",
        ]
        self.assertEqual([item["id"] for item in cases], expected_steps)
        self.assertEqual(
            cases[0]["directory"],
            "grid/stl_to_gmsh/3d/pipe_bend/90_degree/automatic_tetrahedra",
        )
        self.assertEqual(cases[0]["outputs"][0]["path"],
                         "results/pipe_bend_90_degree.msh")
        self.assertEqual(cases[2]["outputs"], [
            {"id": "oneflow_grid", "path": "results/pipe_bend_90_degree.oneflow.grid",
             "parameter_file": "script/grid.txt", "setting": "targetGridFileName"},
            {"id": "bc_region_map", "path": "results/pipe_bend_90_degree_bcRegionMap.txt",
             "parameter_file": "script/grid.txt", "setting": "bcRegionMapFile"},
            {"id": "grid_info", "path": "results/pipe_bend_90_degree.info.json",
             "parameter_file": "script/grid.txt", "setting": "gridInfoFileName"},
        ])
        self.assertEqual(cases[3]["runtime_directory"], "oneflow_cfd_solve")
        self.assertEqual(cases[3]["directory"], "test/cfd/plateuns2dslau2")
        self.assertEqual(cases[3]["inputs"][0]["artifact"], "oneflow_grid")
        self.assertEqual(cases[3]["inputs"][1]["artifact"], "bc_region_map")
        self.assertEqual(cases[3]["inputs"][0]["setting"], "gridFileName")
        self.assertEqual(cases[3]["inputs"][1]["setting"], "bcRegionMapFile")

        lifecycle = mesh_case_lifecycle_summary(repository_root, case_id)
        self.assertEqual(
            lifecycle["verification_stages"],
            ["geometry", "surface_regions", "volume_mesh", "msh_to_cgns",
             "oneflow_ingestion", "cfd_smoke"],
        )
        stage_states = {item["id"]: item for item in lifecycle["stages"]}
        for stage_id in ("oneflow_ingestion", "cfd_smoke"):
            self.assertEqual(stage_states[stage_id]["implementation_status"], "implemented")
            self.assertTrue(stage_states[stage_id]["in_current_scope"])

    def test_native_cube_single_hex_runs_the_complete_cfd_smoke_chain(self):
        cases = build_mesh_validation_cases(self.root, "native_cube", "laminar_flat_plate")
        self.assertEqual([item["id"] for item in cases],
                         ["stl_to_gmsh", "gmsh_to_cgns", "cgns_to_oneflow",
                          "oneflow_cfd_solve"])
        self.assertEqual(cases[0]["directory"], "grid/stl_to_gmsh/3d/cube/single_hex/cpp")
        self.assertEqual(cases[0]["runtime_directory"], "stl_to_gmsh")
        self.assertEqual(cases[3]["directory"], "test/cfd/plateuns2dslau2")
        self.assertEqual(cases[3]["runtime_directory"], "oneflow_cfd_solve")
        self.assertEqual(cases[3]["inputs"][0]["artifact"], "oneflow_grid")
        self.assertEqual(cases[3]["inputs"][1]["artifact"], "bc_region_map")
        self.assertEqual(cases[0]["outputs"][0]["path"], "results/cube_single_hex.msh")
        self.assertEqual(cases[2]["outputs"], [
            {"id": "oneflow_grid", "path": "results/cube_single_hex.oneflow.grid",
             "parameter_file": "script/grid.txt", "setting": "targetGridFileName"},
            {"id": "bc_region_map", "path": "results/cube_single_hex_bcRegionMap.txt",
             "parameter_file": "script/grid.txt", "setting": "bcRegionMapFile"},
            {"id": "grid_info", "path": "results/cube_single_hex.info.json",
             "parameter_file": "script/grid.txt", "setting": "gridInfoFileName"},
        ])

    def test_existing_cylinder_mesh_runs_the_complete_3d_cfd_smoke_chain(self):
        cases = build_mesh_validation_cases(
            self.root, "cylinder_existing_mesh", "laminar_flat_plate"
        )
        self.assertEqual(
            [item["id"] for item in cases],
            ["stl_to_gmsh", "gmsh_to_cgns", "cgns_to_oneflow", "oneflow_cfd_solve"],
        )
        self.assertEqual(cases[0]["directory"],
                         "grid/stl_to_gmsh/3d/cylinder_internal_flow")
        self.assertEqual(cases[0]["outputs"][0]["path"],
                         "results/cylinder_internal_flow.msh")
        self.assertEqual(cases[3]["runtime_directory"], "oneflow_cfd_solve")
        self.assertEqual(cases[3]["directory"], "test/cfd/plateuns2dslau2")
        self.assertEqual(cases[3]["inputs"][0]["artifact"], "oneflow_grid")
        self.assertEqual(cases[3]["inputs"][1]["artifact"], "bc_region_map")

    def test_unit_sphere_cases_run_the_complete_3d_cfd_smoke_chain(self):
        choices = dict(mesh_case_choices(self.root))
        self.assertIn("native_sphere", choices)
        self.assertIn("python_sphere", choices)
        groups = dict(mesh_case_groups(self.root))
        self.assertIn("native_sphere", dict(groups["3D"]))
        self.assertIn("python_sphere", dict(groups["3D"]))

        native = build_mesh_validation_cases(
            self.root, "native_sphere", "laminar_flat_plate"
        )
        python = build_mesh_validation_cases(
            self.root, "python_sphere", "laminar_flat_plate"
        )
        expected_steps = ["stl_to_gmsh", "gmsh_to_cgns", "cgns_to_oneflow", "oneflow_cfd_solve"]
        self.assertEqual([item["id"] for item in native], expected_steps)
        self.assertEqual([item["id"] for item in python], expected_steps)
        self.assertEqual(
            native[0]["directory"],
            "grid/stl_to_gmsh/3d/sphere/automatic_tetrahedra/cpp",
        )
        self.assertEqual(
            native[0]["outputs"][0]["path"],
            "results/sphere_automatic_tetrahedra.msh",
        )
        self.assertEqual(
            python[0]["directory"],
            "grid/stl_to_gmsh/3d/sphere/automatic_tetrahedra/python",
        )
        self.assertEqual(
            python[0]["outputs"][0]["path"],
            "results/sphere_automatic_tetrahedra.msh",
        )
        self.assertEqual(
            native[2]["outputs"][0]["path"],
            "results/sphere_automatic_tetrahedra.oneflow.grid",
        )
        self.assertEqual(native[3]["runtime_directory"], "oneflow_cfd_solve")
        self.assertEqual(python[3]["runtime_directory"], "oneflow_cfd_solve")
        self.assertEqual(native[3]["directory"], "test/cfd/plateuns2dslau2")
        self.assertEqual(python[3]["directory"], "test/cfd/plateuns2dslau2")

    def test_oblique_cut_native_and_python_cases_are_registered_for_validation(self):
        native_id = "stl_to_gmsh_cylinder_internal_flow_oblique_cut"
        python_id = "stl_to_gmsh_python_cylinder_internal_flow_oblique_cut"
        choices = dict(mesh_case_choices(self.root))
        self.assertIn(native_id, choices)
        self.assertIn(python_id, choices)

        native = build_mesh_validation_cases(
            self.root, native_id, "laminar_flat_plate"
        )
        python = build_mesh_validation_cases(
            self.root, python_id, "laminar_flat_plate"
        )
        expected_steps = [
            "stl_to_gmsh", "gmsh_to_cgns", "cgns_to_oneflow",
            "oneflow_cfd_solve",
        ]
        self.assertEqual([item["id"] for item in native], expected_steps)
        self.assertEqual([item["id"] for item in python], expected_steps)
        self.assertEqual(
            native[0]["directory"],
            "grid/stl_to_gmsh/3d/cylinder_internal_flow/oblique_cut/cpp",
        )
        self.assertEqual(
            native[0]["outputs"][0]["path"],
            "results/cylinder_oblique_cut.msh",
        )
        self.assertEqual(
            python[0]["directory"],
            "grid/stl_to_gmsh/3d/cylinder_internal_flow/oblique_cut/python",
        )
        self.assertEqual(
            python[0]["outputs"][0]["path"],
            "results/cylinder_oblique_cut.msh",
        )
        self.assertEqual(native[0]["runtime_directory"], "stl_to_gmsh")
        self.assertEqual(native[3]["directory"], "test/cfd/plateuns2dslau2")

    def test_python_cube_single_hex_runs_the_complete_cfd_smoke_chain(self):
        cube_dir = self.root / "grid/stl_to_gmsh/3d/cube/single_hex/python"
        (cube_dir / "script").mkdir(parents=True, exist_ok=True)
        (cube_dir / "script/layout.json").write_text(
            json.dumps({"mesh_dimension": 3}), encoding="utf-8"
        )
        (cube_dir / "expected.json").write_text(
            json.dumps({"artifact": {"mesh_file": "cube_single_hex.msh"}}),
            encoding="utf-8",
        )
        cases = build_mesh_validation_cases(
            self.root, "python_cube", "laminar_flat_plate"
        )
        self.assertEqual([item["id"] for item in cases],
                         ["stl_to_gmsh", "gmsh_to_cgns", "cgns_to_oneflow",
                          "oneflow_cfd_solve"])
        self.assertEqual(cases[0]["validation"], {
            "kind": "gmsh_msh4", "path": "results/cube_single_hex.msh"
        })
        self.assertEqual(cases[2]["outputs"][1]["path"],
                         "results/cube_single_hex_bcRegionMap.txt")
        self.assertEqual(cases[2]["outputs"][2]["path"],
                         "results/cube_single_hex.info.json")
        self.assertIn("oneflow_cfd_solve", [item["id"] for item in cases])

    def test_first_cfd_template_is_the_existing_laminar_flat_plate_case(self):
        self.assertEqual(cfd_template_choices(), (("laminar_flat_plate", "Laminar Flat Plate"),))

    def test_selected_case_builds_conversion_chain_and_real_cfd_consumer(self):
        cases = build_mesh_validation_cases(self.root, "stl_to_gmsh_square", "laminar_flat_plate")
        self.assertEqual(len(cases), 4)
        self.assertEqual(cases[0]["id"], "stl_to_gmsh")
        self.assertEqual(cases[0]["outputs"][0]["id"], "gmsh_mesh")
        self.assertEqual(cases[0]["outputs"][0]["path"], "results/square.msh")
        self.assertEqual(
            [item["runtime_directory"] for item in cases],
            ["stl_to_gmsh", "gmsh_to_cgns", "cgns_to_oneflow", "oneflow_cfd_solve"],
        )
        self.assertEqual(cases[1]["id"], "gmsh_to_cgns")
        self.assertEqual(cases[1]["depends_on"], ["stl_to_gmsh"])
        self.assertEqual(cases[1]["directory"], "grid/conversion/gmsh_to_cgns")
        self.assertEqual(cases[1]["outputs"][0], {
            "id": "cgns_mesh", "path": "results/square.cgns",
            "parameter_file": "script/grid.txt", "setting": "targetGridFileName",
        })
        self.assertEqual(cases[2]["id"], "cgns_to_oneflow")
        self.assertEqual(cases[2]["depends_on"], ["gmsh_to_cgns"])
        self.assertEqual(cases[2]["directory"], "grid/conversion/cgns_to_oneflow")
        self.assertEqual(cases[2]["outputs"], [
            {"id": "oneflow_grid", "path": "results/square.oneflow.grid",
             "parameter_file": "script/grid.txt", "setting": "targetGridFileName"},
            {"id": "bc_region_map", "path": "results/square_bcRegionMap.txt",
             "parameter_file": "script/grid.txt", "setting": "bcRegionMapFile"},
            {"id": "grid_info", "path": "results/square.info.json",
             "parameter_file": "script/grid.txt", "setting": "gridInfoFileName"},
        ])
        self.assertEqual(cases[3]["id"], "oneflow_cfd_solve")
        self.assertEqual(cases[3]["depends_on"], ["cgns_to_oneflow"])
        self.assertEqual(cases[3]["directory"], "test/cfd/plateuns2dslau2")
        self.assertEqual(cases[3]["inputs"], [
            {"artifact": "oneflow_grid", "setting": "gridFileName",
             "parameter_file": "script/cfd.txt"},
            {"artifact": "bc_region_map", "setting": "bcRegionMapFile",
             "parameter_file": "script/cfd.txt"},
        ])

    def test_run_report_rows_are_ready_for_gui_display(self):
        rows = run_report_stage_rows({
            "stages": {
                "mesh_generation": {"status": "PASS", "message": "MSH checked."},
                "gmsh_to_cgns": {"status": "FAIL", "message": "CGNS writer failed."},
                "cfd_smoke": {"status": "NOT_RUN"},
            }
        })
        self.assertEqual(rows, (
            ("mesh_generation", "PASS", "MSH checked."),
            ("gmsh_to_cgns", "FAIL", "CGNS writer failed."),
            ("cfd_smoke", "NOT_RUN", ""),
        ))
        with self.assertRaisesRegex(ValueError, "must be an object"):
            run_report_stage_rows({"stages": []})

    def test_workflow_graph_maps_to_common_report_states(self):
        result = WorkflowGraphRunResult("failed", (
            CaseExecutionResult("stl_to_gmsh", "succeeded", "MSH generated."),
            CaseExecutionResult("gmsh_to_cgns", "failed", "CGNS writer failed."),
            CaseExecutionResult("cgns_to_oneflow", "skipped", "Blocked."),
        ))
        stages = _mesh_validation_stage_results(result)
        self.assertEqual(stages["mesh_generation"]["status"], "PASS")
        self.assertEqual(stages["gmsh_to_cgns"]["status"], "FAIL")
        self.assertEqual(stages["cgns_to_oneflow"]["status"], "NOT_RUN")
        self.assertEqual(stages["cfd_smoke"]["status"], "NOT_RUN")
        self.assertEqual(stages["mesh_quality"]["status"], "NOT_CHECKED")
        self.assertEqual(stages["cfd_accuracy"]["status"], "NOT_CHECKED")

    def test_existing_mesh_report_names_the_input_stage_without_claiming_generation(self):
        result = WorkflowGraphRunResult("succeeded", (
            CaseExecutionResult("stl_to_gmsh", "succeeded", "Staged existing mesh."),
            CaseExecutionResult("gmsh_to_cgns", "succeeded", "CGNS written."),
            CaseExecutionResult("cgns_to_oneflow", "succeeded", "Grid loaded."),
        ))
        stages = _mesh_validation_stage_results(result, existing_mesh=True)
        self.assertIn("existing_msh_input", stages)
        self.assertNotIn("mesh_generation", stages)
        self.assertEqual(stages["existing_msh_input"]["status"], "PASS")

    def test_successful_handoff_does_not_claim_cfd_or_quality_validation(self):
        result = WorkflowGraphRunResult("succeeded", (
            CaseExecutionResult("stl_to_gmsh", "succeeded", "MSH inspected."),
            CaseExecutionResult("gmsh_to_cgns", "succeeded", "CGNS written."),
            CaseExecutionResult("cgns_to_oneflow", "succeeded", "Grid loaded."),
        ))
        stages = _mesh_validation_stage_results(result)
        self.assertEqual(stages["cfd_smoke"]["status"], "NOT_RUN")
        self.assertEqual(stages["mesh_quality"]["status"], "NOT_CHECKED")
        self.assertEqual(stages["cfd_accuracy"]["status"], "NOT_CHECKED")

    def test_report_records_final_graph_status_and_blocked_steps(self):
        result = WorkflowGraphRunResult("failed", (
            CaseExecutionResult("stl_to_gmsh", "succeeded", "OneFLOW exited successfully."),
            CaseExecutionResult(
                "gmsh_to_cgns", "failed",
                "Declared output artifacts were not produced: results/mesh.cgns",
            ),
            CaseExecutionResult(
                "cgns_to_oneflow", "skipped", "Blocked by unsuccessful dependencies: gmsh_to_cgns",
            ),
            CaseExecutionResult(
                "oneflow_cfd_solve", "skipped", "Blocked by unsuccessful dependencies: cgns_to_oneflow",
            ),
        ))
        self.assertEqual(_graph_step_records(result), [
            {"case_id": "stl_to_gmsh", "status": "succeeded",
             "message": "OneFLOW exited successfully."},
            {"case_id": "gmsh_to_cgns", "status": "failed",
             "message": "Declared output artifacts were not produced: results/mesh.cgns"},
            {"case_id": "cgns_to_oneflow", "status": "skipped",
             "message": "Blocked by unsuccessful dependencies: gmsh_to_cgns"},
            {"case_id": "oneflow_cfd_solve", "status": "skipped",
             "message": "Blocked by unsuccessful dependencies: cgns_to_oneflow"},
        ])

    def test_prolate_spheroid_native_and_python_cases_are_validation_choices(self):
        repository_root = Path(__file__).resolve().parents[2]
        choices = {case_id for case_id, _title in mesh_case_choices(repository_root)}
        for case_id in (
            "stl_to_gmsh_prolate_spheroid_automatic_tetrahedra",
            "stl_to_gmsh_python_prolate_spheroid_automatic_tetrahedra",
        ):
            self.assertIn(case_id, choices)
            graph = build_mesh_validation_cases(
                repository_root, case_id, "laminar_flat_plate"
            )
            self.assertEqual(len(graph), 4)
            self.assertEqual(graph[0]["id"], "stl_to_gmsh")
            self.assertEqual(
                graph[0]["outputs"][0]["path"],
                "results/prolate_spheroid_automatic_tetrahedra.msh",
            )
            self.assertEqual(graph[-1]["id"], "oneflow_cfd_solve")

    def test_downstream_lineage_artifacts_are_collected_from_validation_cases(self):
        steps = [
            {"id": "stl_to_gmsh", "runtime_directory": "stl_to_gmsh",
             "outputs": [{"id": "gmsh_mesh", "path": "results/mesh.msh"}]},
            {"id": "gmsh_to_cgns", "runtime_directory": "gmsh_to_cgns",
             "outputs": [{"id": "cgns_mesh", "path": "results/mesh.cgns"}]},
            {"id": "cgns_to_oneflow", "runtime_directory": "cgns_to_oneflow",
             "outputs": [
                 {"id": "oneflow_grid", "path": "results/mesh.oneflow.grid"},
                 {"id": "bc_region_map", "path": "results/mesh_bcRegionMap.txt"},
             ]},
        ]
        run_root = Path("validation")
        artifacts = _downstream_lineage_artifacts(steps, run_root)
        self.assertEqual(
            artifacts,
            [
                {"role": "CGNS volume mesh", "format": "CGNS",
                 "path": run_root / "gmsh_to_cgns/results/mesh.cgns",
                 "parent": "derived_mesh"},
                {"role": "OneFLOW native grid", "format": "OneFLOW native grid",
                 "path": run_root / "cgns_to_oneflow/results/mesh.oneflow.grid",
                 "parent": "cgns_mesh"},
            ],
        )

    def test_existing_cylinder_mesh_is_a_3d_cfd_smoke_validation_entry(self):
        choices = dict(mesh_case_choices(self.root))
        self.assertIn("cylinder_existing_mesh", choices)
        groups = dict(mesh_case_groups(self.root))
        self.assertIn(
            ("cylinder_existing_mesh", "3D Cylinder Internal Flow (Generated MSH)"),
            groups["3D"],
        )

        steps = build_mesh_validation_cases(
            self.root, "cylinder_existing_mesh", "laminar_flat_plate"
        )
        self.assertEqual(
            [step["id"] for step in steps],
            ["stl_to_gmsh", "gmsh_to_cgns", "cgns_to_oneflow", "oneflow_cfd_solve"],
        )
        self.assertEqual(
            steps[0]["outputs"][0]["path"],
            "results/cylinder_internal_flow.msh",
        )
        self.assertEqual(steps[0]["validation"]["kind"], "gmsh_msh4")

    def test_cylinder_msh_requires_all_three_surface_patches_and_volume(self):
        entry = {
            "required_physical_groups": [
                {"name": "wall", "dimension": 2},
                {"name": "inflow", "dimension": 2},
                {"name": "outflow", "dimension": 2},
                {"name": "fluid_domain", "dimension": 3},
            ]
        }
        mesh_info = {
            "physical_groups": [
                {"name": "wall", "dimension": 2},
                {"name": "inflow", "dimension": 2},
                {"name": "outflow", "dimension": 2},
                {"name": "fluid_domain", "dimension": 3},
            ]
        }
        _validate_required_physical_groups(mesh_info, entry)

        mesh_info["physical_groups"].pop(1)
        with self.assertRaisesRegex(ValueError, "inflow"):
            _validate_required_physical_groups(mesh_info, entry)

        mesh_info["physical_groups"].append({"name": "inflow", "dimension": 3})
        with self.assertRaisesRegex(ValueError, "inflow"):
            _validate_required_physical_groups(mesh_info, entry)

    def test_pipe_bend_is_registered_as_a_generated_mesh_case(self):
        repository_root = Path(__file__).resolve().parents[2]
        case_id = "stl_to_gmsh_pipe_bend_90_degree"
        choices = dict(mesh_case_choices(repository_root))
        self.assertIn(case_id, choices)
        self.assertIn("90-Degree Pipe Bend", choices[case_id])
        catalog = json.loads(
            (repository_root / "grid/case_catalog.json").read_text(encoding="utf-8")
        )
        entry = next(item for item in catalog["cases"] if item["id"] == case_id)
        self.assertEqual(entry["execution"]["kind"], "python")
        self.assertEqual(
            entry["execution"]["script"],
            "grid/stl_to_gmsh/3d/pipe_bend/90_degree/automatic_tetrahedra/run_pipe_bend_case.py",
        )
        self.assertEqual(
            {(item["dimension"], item["name"]) for item in entry["required_physical_groups"]},
            {(2, "wall"), (2, "inflow"), (2, "outflow"), (3, "fluid_domain")},
        )
        categories = dict(mesh_case_group_categories(repository_root))
        self.assertTrue(any(
            item[0] == "STL to Gmsh" and case_id in dict(item[1])
            for item in categories["3D"]
        ))
        self.assertFalse(any(
            item[0] == "Existing Gmsh Mesh" and case_id in dict(item[1])
            for item in categories["3D"]
        ))
        steps = build_mesh_validation_cases(repository_root, case_id, "laminar_flat_plate")
        self.assertEqual(
            [step["id"] for step in steps],
            ["stl_to_gmsh", "gmsh_to_cgns", "cgns_to_oneflow", "oneflow_cfd_solve"],
        )
        self.assertEqual(steps[0]["outputs"][0]["path"], "results/pipe_bend_90_degree.msh")
        command = _mesh_case_command(
            repository_root, entry, Path("OneFLOW.exe"), Path("gmsh.exe"),
            Path("case-directory"),
        )
        self.assertIn("run_pipe_bend_case.py", command[1])
        self.assertEqual(command[2:4], ["--case-directory", "case-directory"])

    def test_pipe_bend_lifecycle_summary_distinguishes_current_and_future_stages(self):
        repository_root = Path(__file__).resolve().parents[2]
        summary = mesh_case_lifecycle_summary(
            repository_root, "stl_to_gmsh_pipe_bend_90_degree"
        )
        self.assertIsNotNone(summary)
        self.assertEqual(
            summary["verification_scope"],
            "geometry_regions_gmsh_volume_mesh_cgns_conversion_oneflow_ingestion_and_cfd_smoke",
        )
        stages = {stage["id"]: stage for stage in summary["stages"]}
        self.assertTrue(stages["geometry"]["in_current_scope"])
        self.assertTrue(stages["surface_regions"]["in_current_scope"])
        self.assertTrue(stages["volume_mesh"]["in_current_scope"])
        self.assertTrue(stages["msh_to_cgns"]["in_current_scope"])
        self.assertFalse(stages["mesh_quality"]["in_current_scope"])
        self.assertEqual(stages["mesh_quality"]["implementation_status"], "implemented")
        self.assertTrue(stages["oneflow_ingestion"]["in_current_scope"])
        self.assertTrue(stages["cfd_smoke"]["in_current_scope"])
        self.assertFalse(stages["cfd_validation"]["in_current_scope"])
        self.assertEqual(stages["geometry"]["implementation_status"], "implemented")
        self.assertEqual(stages["geometry"]["verification_status"], "not_run")
        self.assertEqual(stages["msh_to_cgns"]["implementation_status"], "implemented")
        self.assertEqual(stages["msh_to_cgns"]["verification_status"], "not_run")
        self.assertEqual(stages["oneflow_ingestion"]["implementation_status"], "implemented")
        self.assertEqual(stages["oneflow_ingestion"]["verification_status"], "not_run")
        self.assertEqual(stages["cfd_smoke"]["implementation_status"], "implemented")
        self.assertEqual(stages["cfd_smoke"]["verification_status"], "not_run")
        self.assertFalse(stages["cfd_validation"]["in_current_scope"])
        contract_path = (
            repository_root / "grid" / "stl_to_gmsh" / "3d" / "pipe_bend"
            / "90_degree" / "automatic_tetrahedra" / "case-lifecycle.json"
        )
        contract = json.loads(contract_path.read_text(encoding="utf-8"))
        contract_stages = {stage["id"]: stage for stage in contract["stages"]}
        self.assertEqual(
            contract_stages["oneflow_ingestion"]["depends_on"],
            ["msh_to_cgns"],
        )
        self.assertEqual(
            contract["quality_policy"]["role"], "diagnostic"
        )
        self.assertFalse(contract["quality_policy"]["blocking"])
        self.assertFalse(contract["quality_policy"]["affects_oneflow_ingestion"])
        self.assertFalse(contract["quality_policy"]["affects_cfd_smoke"])
        self.assertIn(
            {"id": "native_cpp_geometry", "kind": "cpp", "status": "implemented"},
            stages["geometry"]["implementations"],
        )
        self.assertIn(
            {"id": "native_cpp_region_export", "kind": "cpp", "status": "implemented"},
            stages["surface_regions"]["implementations"],
        )
        self.assertEqual(stages["cfd_validation"]["implementation_status"], "planned")

    def test_cases_without_lifecycle_contract_have_no_milestone_summary(self):
        repository_root = Path(__file__).resolve().parents[2]
        self.assertIsNone(
            mesh_case_lifecycle_summary(repository_root, "stl_to_gmsh_square")
        )
        legacy_steps = build_mesh_validation_cases(
            repository_root, "stl_to_gmsh_square", "laminar_flat_plate"
        )
        self.assertEqual(
            [step["id"] for step in legacy_steps],
            ["stl_to_gmsh", "gmsh_to_cgns", "cgns_to_oneflow", "oneflow_cfd_solve"],
        )
        self.assertTrue(
            mesh_case_requires_oneflow_executable(
                repository_root, "stl_to_gmsh_square"
            )
        )
        # The pipe-bend lifecycle now includes the shared native-grid
        # ingestion and bounded CFD smoke stages.
        self.assertTrue(
            mesh_case_requires_oneflow_executable(
                repository_root, "stl_to_gmsh_pipe_bend_90_degree"
            )
        )
        pipe_bend_steps = build_mesh_validation_cases(
            repository_root, "stl_to_gmsh_pipe_bend_90_degree",
            "laminar_flat_plate",
        )
        self.assertEqual(
            [step["id"] for step in pipe_bend_steps],
            ["stl_to_gmsh", "gmsh_to_cgns", "cgns_to_oneflow", "oneflow_cfd_solve"],
        )
        self.assertTrue(
            any(step["id"] == "cgns_to_oneflow" for step in pipe_bend_steps)
        )
        self.assertTrue(
            any(step["id"] == "oneflow_cfd_solve" for step in pipe_bend_steps)
        )

    def test_existing_mesh_case_resolves_only_inside_user_mesh_workspace(self):
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary) / "mesh-generation"
            paths = type("Paths", (), {"mesh_generation": workspace})()
            with mock.patch(
                "tools.oneflow_gui.mesh_validation.UserDataPaths.discover",
                return_value=paths,
            ):
                resolved = existing_mesh_case_paths(self.root, "cylinder_existing_mesh")
            expected_project = (workspace / "cylinder_internal_flow").resolve()
            self.assertEqual(
                resolved,
                (
                    expected_project / "results/cylinder_internal_flow.msh",
                    expected_project,
                ),
            )

            invalid_workspace = Path(temporary) / "invalid-mesh-generation"
            invalid_project = invalid_workspace / "invalid_cylinder"
            (invalid_project / "results").mkdir(parents=True)
            (invalid_project / "results/escaped.msh").write_text(
                "mesh", encoding="ascii"
            )
            (invalid_project / "project.json").write_text(
                json.dumps({
                    "schema": "oneflow.mesh-project", "version": 1,
                    "geometry": "cylinder_internal_flow",
                    "output_mesh": "../../outside.msh",
                }), encoding="utf-8"
            )
            invalid_paths = type(
                "Paths", (), {"mesh_generation": invalid_workspace}
            )()
            with mock.patch(
                "tools.oneflow_gui.mesh_validation.UserDataPaths.discover",
                return_value=invalid_paths,
            ):
                ignored_invalid_project = existing_mesh_case_paths(
                    Path(__file__).resolve().parents[2],
                    "stl_to_gmsh_existing_mesh_3d_cylinder_internal_flow",
                )
            self.assertEqual(
                ignored_invalid_project,
                (
                    (invalid_workspace / "cylinder_internal_flow"
                     / "results/cylinder_internal_flow.msh").resolve(),
                    (invalid_workspace / "cylinder_internal_flow").resolve(),
                ),
            )

            catalog = json.loads(
                (self.root / "grid/case_catalog.json").read_text(encoding="utf-8")
            )
            cylinder = next(
                x for x in catalog["mesh_validation_inputs"]
                if x["id"] == "cylinder_existing_mesh"
            )
            cylinder["execution"]["mesh_file"] = "../../outside.msh"
            (self.root / "grid/case_catalog.json").write_text(
                json.dumps(catalog), encoding="utf-8"
            )
            with mock.patch(
                "tools.oneflow_gui.mesh_validation.UserDataPaths.discover",
                return_value=paths,
            ):
                with self.assertRaisesRegex(ValueError, "safe relative"):
                    existing_mesh_case_paths(self.root, "cylinder_existing_mesh")

    def test_unknown_mesh_case_and_template_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "Unknown mesh-generation Case"):
            build_mesh_validation_cases(self.root, "missing", "laminar_flat_plate")
        with self.assertRaisesRegex(ValueError, "Unknown CFD template"):
            build_mesh_validation_cases(self.root, "stl_to_gmsh_square", "missing")



    def test_existing_mesh_is_copied_to_an_isolated_stage(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "project" / "mesh.msh"
            source.parent.mkdir()
            source.write_text("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n", encoding="ascii")
            destination = root / "validation" / "run-1" / "producer" / "results" / "mesh.msh"
            _copy_existing_mesh(source, destination)
            self.assertEqual(destination.read_bytes(), source.read_bytes())
            self.assertTrue(source.is_file())
            self.assertNotEqual(source.resolve(), destination.resolve())

    def test_existing_mesh_copy_rejects_missing_or_non_msh_sources(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with self.assertRaisesRegex(FileNotFoundError, "not found"):
                _copy_existing_mesh(root / "missing.msh", root / "stage" / "mesh.msh")
            invalid = root / "mesh.stl"
            invalid.write_text("solid mesh\nendsolid mesh\n", encoding="ascii")
            with self.assertRaisesRegex(ValueError, "extension"):
                _copy_existing_mesh(invalid, root / "stage" / "mesh.msh")

    def test_project_validation_workspace_flattens_legacy_run_directory(self):
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary) / "pipe_bend_90_degree" / "validation"
            legacy = workspace / "run-20261125T153212Z-a3ce16f7"
            (legacy / "gmsh_to_cgns" / "results").mkdir(parents=True)
            (legacy / "mesh-validation-report.json").write_text(
                '{"status":"SUCCEEDED"}', encoding="utf-8"
            )
            (legacy / "gmsh_to_cgns" / "results" / "mesh.cgns").write_text(
                "cgns", encoding="ascii"
            )
            prepared = prepare_mesh_validation_run_workspace(
                workspace, preserve_history=False, run_id="next-run"
            )
            self.assertTrue(prepared.samefile(workspace))
            self.assertTrue((workspace / "mesh-validation-report.json").is_file())
            self.assertTrue(
                (workspace / "gmsh_to_cgns" / "results" / "mesh.cgns").is_file()
            )
            self.assertFalse(legacy.exists())
            self.assertFalse(any(
                child.is_dir() and child.name.startswith("run-")
                for child in workspace.iterdir()
            ))

    def test_gui_validation_run_uses_stable_case_directory_without_run_subfolder(self):
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary) / "stl_to_gmsh" / "3d" / "pipe_bend" / "90_degree" / "automatic_tetrahedra"
            prepared = prepare_mesh_validation_run_workspace(
                workspace, preserve_history=True, run_id="20261125T191006Z-ac08120a"
            )
            self.assertTrue(prepared.samefile(workspace))
            self.assertTrue(workspace.is_dir())
            self.assertFalse(any(
                child.is_dir() and child.name.startswith("run-")
                for child in workspace.iterdir()
            ))
            (workspace / "mesh-validation-report.json").write_text(
                '{"status":"SUCCEEDED"}', encoding="utf-8"
            )
            prepared_again = prepare_mesh_validation_run_workspace(
                workspace, preserve_history=True, run_id="next-run"
            )
            self.assertTrue(prepared_again.samefile(workspace))
            self.assertFalse((workspace / "mesh-validation-report.json").exists())

    def test_workspace_uses_mesh_case_relative_directory(self):
        paths = type(
            "Paths", (),
            {"mesh_validation": self.root / ".oneflow" / "mesh-validation"},
        )()
        with mock.patch(
            "tools.oneflow_gui.mesh_validation.UserDataPaths.discover",
            return_value=paths,
        ):
            workspace = mesh_validation_workspace(self.root, "stl_to_gmsh_square")
        self.assertEqual(
            workspace,
            self.root / ".oneflow" / "mesh-validation" / "stl_to_gmsh" / "2d" / "square" / "single_quad" / "cpp",
        )


    def test_existing_mesh_validation_uses_project_validation_directory_without_run_stamp(self):
        with tempfile.TemporaryDirectory() as temporary:
            project = Path(temporary) / "mesh-generation" / "sphere"
            project.mkdir(parents=True)
            workspace = project_mesh_validation_workspace(project)
            self.assertEqual(workspace, (project / "validation").resolve())
            self.assertEqual(workspace.name, "validation")
            self.assertNotIn("T", workspace.name)

    def test_workspace_rejects_cases_outside_grid(self):
        catalog = json.loads(
            (self.root / "grid/case_catalog.json").read_text(encoding="utf-8")
        )
        catalog["cases"][0]["directory"] = "../outside"
        (self.root / "grid/case_catalog.json").write_text(
            json.dumps(catalog), encoding="utf-8"
        )
        with self.assertRaisesRegex(ValueError, "below grid"):
            mesh_validation_workspace(self.root, "stl_to_gmsh_square")

    def test_reset_clears_only_the_selected_case_directory(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            selected = root / "mesh-validation/stl_to_gmsh/2d/square"
            sibling = root / "mesh-validation/stl_to_gmsh/2d/square_2x2"
            selected.mkdir(parents=True)
            sibling.mkdir(parents=True)
            (selected / "stale.txt").write_text("old", encoding="utf-8")
            (sibling / "keep.txt").write_text("keep", encoding="utf-8")
            reset_mesh_validation_workspace(selected)
            self.assertFalse((selected / "stale.txt").exists())
            self.assertTrue((sibling / "keep.txt").is_file())


    def test_explicit_stage_root_uses_direct_output_and_resets_stale_artifacts(self):
        with tempfile.TemporaryDirectory() as temporary:
            stage_root = Path(temporary) / "oneflow-generic-curved-shell-validation"
            stage_root.mkdir()
            (stage_root / "stale-report.json").write_text("stale", encoding="utf-8")
            (stage_root / "old-stage").mkdir()
            run_root = prepare_mesh_validation_run_workspace(
                stage_root, preserve_history=False, run_id="ci-run"
            )
            self.assertEqual(run_root, stage_root.resolve())
            self.assertFalse((run_root / "stale-report.json").exists())
            self.assertFalse((run_root / "old-stage").exists())
            (run_root / "mesh-validation-report.json").write_text(
                "{}", encoding="utf-8"
            )
            (run_root / "gmsh_to_cgns").mkdir()
            (run_root / "cgns_to_oneflow").mkdir()
            self.assertTrue((stage_root / "mesh-validation-report.json").is_file())
            self.assertTrue((stage_root / "gmsh_to_cgns").is_dir())
            self.assertTrue((stage_root / "cgns_to_oneflow").is_dir())

    def test_each_run_gets_an_isolated_directory_and_preserves_history(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary) / "project" / "validation"
            old_run = base / "run-20261115T173711-ae8d11f2"
            old_run.mkdir(parents=True)
            (old_run / "mesh-validation-report.json").write_text(
                json.dumps({"status": "success", "mesh_case_id": "old_case"}), encoding="utf-8"
            )
            first = create_mesh_validation_run_workspace(base, "first123")
            second = create_mesh_validation_run_workspace(base, "second123")
            self.assertNotEqual(first, second)
            self.assertEqual(first.parent, base.resolve())
            self.assertTrue((old_run / "mesh-validation-report.json").is_file())

    def test_run_history_lists_reports_and_interrupted_runs(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary) / "validation"
            complete, interrupted = base / "run-complete", base / "run-interrupted"
            complete.mkdir(parents=True)
            interrupted.mkdir()
            (complete / "mesh-validation-report.json").write_text(
                json.dumps({"status": "success", "mesh_case_id": "sphere_case",
                            "started_at": "2026-11-24T10:00:00Z"}), encoding="utf-8"
            )
            (interrupted / "mesh-validation-run.json").write_text(
                json.dumps({"schema": "oneflow.mesh-validation-run", "version": 1,
                            "mesh_case_id": "cylinder_case",
                            "started_at_utc": "2026-11-24T11:00:00Z"}), encoding="utf-8"
            )
            nested = base / "stl_to_gmsh" / "3d" / "sphere" / "run-nested"
            nested.mkdir(parents=True)
            (nested / "mesh-validation-run.json").write_text(
                json.dumps({"schema": "oneflow.mesh-validation-run", "version": 1,
                            "mesh_case_id": "nested_sphere",
                            "started_at_utc": "2026-11-24T12:00:00Z"}), encoding="utf-8"
            )
            records = list_mesh_validation_runs(base)
            self.assertEqual([item["mesh_case_id"] for item in records],
                             ["nested_sphere", "cylinder_case", "sphere_case"])
            self.assertEqual(records[0]["status"], "RECOVERABLE")
            self.assertEqual(records[1]["status"], "RECOVERABLE")
            self.assertEqual(records[2]["status"], "SUCCESS")


if __name__ == "__main__":
    unittest.main()
