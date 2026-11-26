"""Regression tests for the 90-degree pipe-bend STL geometry asset."""
from __future__ import annotations

import importlib.util
import hashlib
import json
import os
import shutil
import subprocess
import math
import tempfile
import unittest
from pathlib import Path

from grid.stl_import.reader import read_stl
from grid.stl_import.validation import analyze_surface
from grid.stl_to_gmsh.converter import convert

ROOT = Path(__file__).resolve().parents[4]
GENERATOR = (
    ROOT / "grid/stl_to_gmsh/3d/pipe_bend/90_degree"
    / "geometry/python/generate_pipe_bend.py"
)
SPEC = importlib.util.spec_from_file_location("oneflow_pipe_bend_generator", GENERATOR)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class PipeBendGeometryTest(unittest.TestCase):
    def test_default_bend_is_watertight_and_has_three_semantic_patches(self):
        with tempfile.TemporaryDirectory() as temporary:
            result = MODULE.generate(Path(temporary) / "pipe_bend_90_degree.stl")
            surface = result["surface"]
            self.assertTrue(surface["watertight"])
            self.assertEqual(surface["connected_component_count"], 1)
            self.assertEqual(surface["boundary_or_nonmanifold_edge_count"], 0)
            self.assertEqual(surface["inconsistently_oriented_edge_count"], 0)
            self.assertGreater(surface["signed_volume"], 0.0)
            self.assertLess(surface["relative_volume_error"], 0.01)
            patches = result["patch_triangle_indices"]
            self.assertEqual(len(patches["inflow"]), 48)
            self.assertEqual(len(patches["outflow"]), 48)
            self.assertEqual(len(patches["wall"]), 2 * 48 * (12 + 24 + 12))
            self.assertEqual(sum(map(len, patches.values())), surface["triangle_count"])
            manifest = json.loads(result["manifest"].read_text(encoding="utf-8"))
            self.assertEqual(manifest["asset_id"], "geometry.pipe_bend.90_degree")
            self.assertEqual(manifest["definition"]["parameters"]["radius"], 0.5)
            self.assertEqual(len(manifest["definition_sha256"]), 64)
            self.assertEqual(len(manifest["generator"]["source_sha256"]), 64)
            self.assertEqual(
                manifest["intended_gmsh_physical_groups"]["surface"],
                ["wall", "inflow", "outflow"],
            )
            self.assertEqual(
                manifest["intended_gmsh_physical_groups"]["volume"],
                ["fluid_domain"],
            )
            region_manifest_path = Path(temporary) / "region_manifest.json"
            region_manifest = json.loads(region_manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(region_manifest["schema"], "oneflow.stl-surface-regions")
            self.assertEqual(
                region_manifest["closed_surface_sha256"],
                hashlib.sha256(result["stl"].read_bytes()).hexdigest(),
            )
            for name in ("wall", "inflow", "outflow"):
                record = region_manifest["regions"][name]
                region_file = Path(temporary) / record["file"]
                self.assertTrue(region_file.is_file())
                self.assertEqual(len(read_stl(region_file)), record["triangle_count"])
                self.assertEqual(
                    hashlib.sha256(region_file.read_bytes()).hexdigest(), record["sha256"]
                )

            layout = Path(temporary) / "layout.json"
            layout.write_text(json.dumps({
                "mesh_dimension": 3,
                "volume_mesh": "region_aware_closed_surface_automatic_tetrahedra",
                "region_manifest": str(region_manifest_path),
                "mesh_size_min": 0.2,
                "mesh_size_max": 0.2,
            }), encoding="utf-8")
            geo = Path(temporary) / "pipe_bend.geo"
            conversion = convert(result["stl"], layout, geo)
            self.assertEqual(conversion["geometry"], "region_aware_generic_stl")
            self.assertEqual(
                conversion["region_triangle_counts"],
                {"wall": 4608, "inflow": 48, "outflow": 48},
            )
            geometry_text = geo.read_text(encoding="utf-8")
            self.assertIn("Coherence Mesh;", geometry_text)
            self.assertIn("Surface Loop(1) = {wall[], inflow[], outflow[]};", geometry_text)
            self.assertIn('Physical Surface("inflow")', geometry_text)
            self.assertIn('Physical Surface("outflow")', geometry_text)

    def test_analytic_volume_and_turn_angle_are_parameterized(self):
        with tempfile.TemporaryDirectory() as temporary:
            params = {
                "radius": 0.25, "bend_radius": 0.8,
                "inlet_length": 1.0, "outlet_length": 2.0,
                "angle_degrees": 45.0, "circumferential_divisions": 24,
                "bend_divisions": 16, "inlet_divisions": 8, "outlet_divisions": 8,
            }
            result = MODULE.generate(Path(temporary) / "bend45.stl", **params)
            expected = math.pi * params["radius"]**2 * (
                params["inlet_length"] + params["outlet_length"]
                + params["bend_radius"] * math.radians(params["angle_degrees"])
            )
            self.assertAlmostEqual(result["surface"]["analytic_volume"], expected)
            self.assertLess(result["surface"]["relative_volume_error"], 0.02)

    def test_region_manifest_rejects_wrong_patch_geometry_even_when_counts_match(self):
        from grid.stl_to_gmsh.converter import Stl2DError

        with tempfile.TemporaryDirectory() as temporary:
            stl = Path(temporary) / "pipe_bend_90_degree.stl"
            MODULE.generate(stl)
            region_manifest_path = stl.parent / "region_manifest.json"
            manifest = json.loads(region_manifest_path.read_text(encoding="utf-8"))
            inflow = stl.parent / manifest["regions"]["inflow"]["file"]
            outflow = stl.parent / manifest["regions"]["outflow"]["file"]
            shutil.copyfile(outflow, inflow)
            manifest["regions"]["inflow"]["sha256"] = hashlib.sha256(
                inflow.read_bytes()
            ).hexdigest()
            region_manifest_path.write_text(
                json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            layout = stl.parent / "layout.json"
            layout.write_text(json.dumps({
                "mesh_dimension": 3,
                "volume_mesh": "region_aware_closed_surface_automatic_tetrahedra",
                "region_manifest": str(region_manifest_path),
                "mesh_size_min": 0.2,
                "mesh_size_max": 0.2,
            }), encoding="utf-8")
            with self.assertRaisesRegex(Stl2DError, "exact triangle partition"):
                convert(stl, layout, stl.parent / "invalid.geo")

    @unittest.skipUnless(
        os.environ.get("ONEFLOW_RUN_GMSH_TESTS") == "1"
        and (shutil.which("gmsh") or shutil.which("gmsh.exe")),
        "Set ONEFLOW_RUN_GMSH_TESTS=1 and install Gmsh for this integration test",
    )
    def test_gmsh_generates_pipe_bend_volume_with_named_physical_groups(self):
        from grid.stl_to_gmsh.generate_volume_mesh import generate

        gmsh = shutil.which("gmsh") or shutil.which("gmsh.exe")
        with tempfile.TemporaryDirectory() as temporary:
            stl = Path(temporary) / "pipe_bend_90_degree.stl"
            MODULE.generate(stl)
            report = generate(
                stl, Path(temporary) / "pipe_bend_90_degree.msh",
                gmsh, 0.25, 0.25,
            )
            self.assertEqual(report["stages"]["geometry"]["status"], "PASS")
            self.assertEqual(report["stages"]["gmsh"]["status"], "PASS")
            self.assertEqual(report["stages"]["physical_groups"]["status"], "PASS")
            groups = {
                (item["dimension"], item["name"])
                for item in report["mesh_summary"]["physical_groups"]
            }
            self.assertEqual(groups, {
                (2, "wall"), (2, "inflow"), (2, "outflow"), (3, "fluid_domain")
            })
            self.assertGreater(report["mesh_summary"]["tetrahedron_elements"], 0)

    def _assert_cpp_generator_matches_python_geometry(self, cpp_stl):
        self.assertTrue(cpp_stl.is_file())
        cpp_region_manifest_path = cpp_stl.parent / "region_manifest.json"
        cpp_geometry_manifest_path = cpp_stl.with_suffix(".geometry.json")
        self.assertTrue(cpp_region_manifest_path.is_file())
        self.assertTrue(cpp_geometry_manifest_path.is_file())
        cpp_region_manifest = json.loads(
            cpp_region_manifest_path.read_text(encoding="utf-8")
        )
        cpp_geometry_manifest = json.loads(
            cpp_geometry_manifest_path.read_text(encoding="utf-8")
        )
        self.assertEqual(
            cpp_region_manifest["schema"], "oneflow.stl-surface-regions"
        )
        self.assertEqual(
            cpp_region_manifest["closed_surface_sha256"],
            hashlib.sha256(cpp_stl.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            cpp_geometry_manifest["schema"], "oneflow.stl-geometry-manifest"
        )
        self.assertEqual(
            cpp_geometry_manifest["generator"]["name"],
            "oneflow_pipe_bend_geometry",
        )
        self.assertEqual(
            cpp_geometry_manifest["generator"]["source"],
            "generate_pipe_bend.cpp",
        )
        self.assertEqual(
            len(cpp_geometry_manifest["generator"]["source_sha256"]), 64
        )
        self.assertEqual(
            cpp_geometry_manifest["stl_sha256"],
            hashlib.sha256(cpp_stl.read_bytes()).hexdigest(),
        )
        for name in ("wall", "inflow", "outflow"):
            record = cpp_region_manifest["regions"][name]
            region_file = cpp_stl.parent / record["file"]
            self.assertTrue(region_file.is_file())
            self.assertEqual(
                hashlib.sha256(region_file.read_bytes()).hexdigest(),
                record["sha256"],
            )
            self.assertEqual(
                len(read_stl(region_file)), record["triangle_count"]
            )
        with tempfile.TemporaryDirectory() as temporary:
            python_stl = MODULE.generate(
                Path(temporary) / "pipe_bend_python.stl"
            )["stl"]
            layout = Path(temporary) / "cpp-layout.json"
            layout.write_text(json.dumps({
                "mesh_dimension": 3,
                "volume_mesh": "region_aware_closed_surface_automatic_tetrahedra",
                "region_manifest": str(cpp_region_manifest_path),
                "mesh_size_min": 0.2,
                "mesh_size_max": 0.2,
            }), encoding="utf-8")
            conversion = convert(
                cpp_stl, layout, Path(temporary) / "cpp-pipe-bend.geo"
            )
            self.assertEqual(
                conversion["geometry"], "region_aware_generic_stl"
            )
            self.assertEqual(
                set(conversion["region_triangle_counts"]),
                {"wall", "inflow", "outflow"},
            )
            cpp_metrics = analyze_surface(read_stl(cpp_stl))
            python_metrics = analyze_surface(read_stl(python_stl))
            for metrics in (cpp_metrics, python_metrics):
                self.assertEqual(metrics["status"], "valid")
                self.assertTrue(metrics["closed"])
                self.assertTrue(metrics["edge_manifold"])
                self.assertTrue(metrics["consistently_oriented"])
                self.assertEqual(metrics["connected_component_count"], 1)
            self.assertEqual(cpp_metrics["triangle_count"], python_metrics["triangle_count"])
            self.assertEqual(
                cpp_metrics["welded_vertex_count"], python_metrics["welded_vertex_count"]
            )
            self.assertAlmostEqual(
                cpp_metrics["surface_area"], python_metrics["surface_area"], delta=1.0e-4
            )
            cpp_volume = cpp_metrics["components"][0]["signed_volume"]
            python_volume = python_metrics["components"][0]["signed_volume"]
            self.assertAlmostEqual(cpp_volume, python_volume, delta=1.0e-5)
            for axis in ("x", "y", "z"):
                for index in (0, 1):
                    self.assertAlmostEqual(
                        cpp_metrics["coordinate_bounds"][axis][index],
                        python_metrics["coordinate_bounds"][axis][index],
                        delta=1.0e-6,
                    )

    def test_cpp_generator_matches_python_geometry(self):
        configured_stl = os.environ.get("ONEFLOW_PIPE_BEND_CPP_STL")
        if configured_stl:
            self._assert_cpp_generator_matches_python_geometry(Path(configured_stl))
            return
        source = (
            ROOT / "grid/stl_to_gmsh/3d/pipe_bend/90_degree"
            / "geometry/cpp/generate_pipe_bend.cpp"
        )
        compiler = next(
            (shutil.which(name) for name in ("c++", "g++", "clang++") if shutil.which(name)),
            None,
        )
        if not compiler:
            self.skipTest("A C++17 compiler or ONEFLOW_PIPE_BEND_CPP_STL is required")
        with tempfile.TemporaryDirectory() as build_directory:
            executable = Path(build_directory) / "pipe_bend_generator"
            subprocess.run(
                [compiler, "-std=c++17", str(source), "-o", str(executable)],
                check=True, capture_output=True, text=True,
            )
            cpp_stl = Path(build_directory) / "cpp" / "pipe_bend_90_degree.stl"
            cpp_stl.parent.mkdir(parents=True)
            subprocess.run(
                [str(executable), "--output", str(cpp_stl)],
                check=True, capture_output=True, text=True,
            )
            self._assert_cpp_generator_matches_python_geometry(cpp_stl)


    def test_pipe_bend_catalog_exposes_both_geometry_implementations(self):
        catalog = json.loads((ROOT / "grid/case_catalog.json").read_text(encoding="utf-8"))
        entry = next(
            item for item in catalog["cases"]
            if item["id"] == "stl_to_gmsh_pipe_bend_90_degree"
        )
        implementations = {
            item["id"]: item for item in entry["execution"]["implementations"]
        }
        self.assertEqual(set(implementations), {"python_reference", "native_cpp_geometry"})
        self.assertTrue(implementations["python_reference"]["default"])
        self.assertEqual(implementations["python_reference"]["kind"], "python")
        self.assertEqual(implementations["native_cpp_geometry"]["kind"], "cpp")
        self.assertEqual(
            implementations["native_cpp_geometry"]["compiler_standard"], "c++17"
        )
        self.assertTrue(
            (ROOT / implementations["native_cpp_geometry"]["script"]).is_file()
        )


    @unittest.skipUnless(
        os.environ.get("ONEFLOW_RUN_GMSH_TESTS") == "1"
        and (shutil.which("gmsh") or shutil.which("gmsh.exe")),
        "Set ONEFLOW_RUN_GMSH_TESTS=1 and install Gmsh for this integration test",
    )
    def test_cpp_case_runner_generates_msh4_with_shared_physical_groups(self):
        runner_path = (
            ROOT / "grid/stl_to_gmsh/3d/pipe_bend/90_degree"
            / "automatic_tetrahedra/run_pipe_bend_cpp_case.py"
        )
        runner_spec = importlib.util.spec_from_file_location(
            "oneflow_pipe_bend_cpp_case_integration_runner", runner_path
        )
        self.assertIsNotNone(runner_spec)
        self.assertIsNotNone(runner_spec.loader)
        runner = importlib.util.module_from_spec(runner_spec)
        runner_spec.loader.exec_module(runner)
        try:
            compiler = runner.discover_cpp_compiler()
        except RuntimeError as error:
            self.skipTest(str(error))
        gmsh = shutil.which("gmsh") or shutil.which("gmsh.exe")
        with tempfile.TemporaryDirectory() as temporary:
            mesh = runner.run_case(
                Path(temporary) / "case",
                gmsh=gmsh,
                mesh_size=0.25,
                compiler=compiler,
            )
            self.assertTrue(mesh.is_file())
            self.assertGreater(mesh.stat().st_size, 0)
            report = json.loads(
                mesh.with_suffix(".report.json").read_text(encoding="utf-8")
            )
            self.assertEqual(report["stages"]["geometry"]["status"], "PASS")
            self.assertEqual(report["stages"]["gmsh"]["status"], "PASS")
            self.assertEqual(report["stages"]["physical_groups"]["status"], "PASS")
            groups = {
                (item["dimension"], item["name"])
                for item in report["mesh_summary"]["physical_groups"]
            }
            self.assertEqual(
                groups,
                {(2, "wall"), (2, "inflow"), (2, "outflow"), (3, "fluid_domain")},
            )
            self.assertGreater(report["mesh_summary"]["tetrahedron_elements"], 0)

    def test_cpp_case_runner_builds_the_native_geometry_in_workspace(self):
        runner_path = (
            ROOT / "grid/stl_to_gmsh/3d/pipe_bend/90_degree"
            / "automatic_tetrahedra/run_pipe_bend_cpp_case.py"
        )
        runner_spec = importlib.util.spec_from_file_location(
            "oneflow_pipe_bend_cpp_case_runner", runner_path
        )
        self.assertIsNotNone(runner_spec)
        self.assertIsNotNone(runner_spec.loader)
        runner = importlib.util.module_from_spec(runner_spec)
        runner_spec.loader.exec_module(runner)
        try:
            compiler = runner.discover_cpp_compiler()
        except RuntimeError as error:
            self.skipTest(str(error))
        with tempfile.TemporaryDirectory() as temporary:
            workspace = Path(temporary)
            executable = runner.build_cpp_generator(
                workspace / "build", compiler=compiler
            )
            self.assertTrue(executable.is_file())
            stl = workspace / "input" / "pipe_bend_90_degree.stl"
            stl.parent.mkdir(parents=True)
            subprocess.run(
                [str(executable), "--output", str(stl)],
                cwd=ROOT, check=True, capture_output=True, text=True,
            )
            self.assertTrue(stl.is_file())
            self.assertTrue((stl.parent / "region_manifest.json").is_file())
            self.assertTrue(stl.with_suffix(".geometry.json").is_file())
            self._assert_cpp_generator_matches_python_geometry(stl)

    def test_invalid_geometry_and_source_tree_output_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(ValueError, "bend_radius"):
                MODULE.build_surface(radius=0.5, bend_radius=0.5)
            with self.assertRaisesRegex(ValueError, "outside the OneFLOW source tree"):
                MODULE.generate(ROOT / "build" / "forbidden_pipe_bend.stl")


if __name__ == "__main__":
    unittest.main()
