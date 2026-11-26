#!/usr/bin/env python3
"""Unit tests for the OneFLOW grid case catalog."""

from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from mesh_cases import load_catalog, validate_case, validate_profile
from case_lifecycle import validate_lifecycle_contract


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


class MeshCaseCatalogTest(unittest.TestCase):
    def test_repository_catalog_entries_resolve(self) -> None:
        catalog = load_catalog(REPOSITORY_ROOT / "grid" / "case_catalog.json")
        self.assertGreaterEqual(len(catalog["cases"]), 4)
        profile_ids = {
            entry["id"] for entry in catalog["cases"] if entry.get("profile")
        }
        self.assertTrue({
            "minimal_rectangle",
            "minimal_rectangle_vertical_midpoints",
            "minimal_cube",
            "two_blocks_shared_face_corners_only",
            "two_blocks_shared_edge_endpoints_only",
            "two_blocks_shared_edge_one_interior_node",
        } <= profile_ids)
        layout_directories = {
            "minimal_cube": "grid/mesh-layout/3d/minimal_cube",
            "minimal_rectangle": "grid/mesh-layout/2d/minimal_rectangle",
            "minimal_rectangle_vertical_midpoints": "grid/mesh-layout/2d/minimal_rectangle_vertical_midpoints",
            "two_blocks_shared_face_corners_only": "grid/mesh-layout/3d/two_blocks_shared_face_corners_only",
            "two_blocks_shared_edge_endpoints_only": "grid/mesh-layout/2d/two_blocks_shared_edge_endpoints_only",
            "two_blocks_shared_edge_one_interior_node": "grid/mesh-layout/2d/two_blocks_shared_edge_one_interior_node",
        }
        for entry in catalog["cases"]:
            if entry["id"] in layout_directories:
                self.assertEqual(entry["directory"], layout_directories[entry["id"]])
        errors = [
            error
            for entry in catalog["cases"]
            for error in validate_case(REPOSITORY_ROOT, entry)
        ]
        self.assertEqual(errors, [])

    def test_two_quads_is_a_gmsh_to_cgns_route_fixture_not_a_catalog_case(self) -> None:
        catalog = load_catalog(REPOSITORY_ROOT / "grid" / "case_catalog.json")
        case_ids = {case["id"] for case in catalog["cases"]}
        self.assertNotIn("two_quads", case_ids)

        conversion_root = REPOSITORY_ROOT / "grid" / "conversion"
        self.assertTrue((conversion_root / "README.md").is_file())
        for template in ("cgns_to_oneflow", "gmsh_to_cgns", "gmsh_to_oneflow"):
            self.assertTrue((conversion_root / template / "README.md").is_file())

        fixture_relative = "grid/conversion/gmsh_to_cgns/grid/two_quads.msh"
        fixture = REPOSITORY_ROOT / fixture_relative
        self.assertTrue(fixture.is_file())
        mesh_text = fixture.read_text(encoding="ascii")
        self.assertIn("$MeshFormat\n4.1 0 8", mesh_text)
        self.assertIn('1 1 "inflow"', mesh_text)
        self.assertIn('1 2 "outflow"', mesh_text)
        self.assertIn('1 3 "symmetry"', mesh_text)
        self.assertIn('1 4 "farfield"', mesh_text)
        self.assertIn('2 5 "fluid_domain"', mesh_text)

        with (REPOSITORY_ROOT / "grid" / "interoperability" / "conversion_matrix.json").open(
            "r", encoding="utf-8"
        ) as stream:
            matrix = json.load(stream)
        fixture_cases = {
            case["id"]: case for case in matrix.get("cases", [])
        }
        entry = fixture_cases["gmsh_to_cgns_two_quads"]
        self.assertEqual(entry["fixture"], fixture_relative)
        self.assertEqual(entry["routes"], ["gmsh->cgns"])
        self.assertEqual(entry["execution"], "automated_integration_test")
        self.assertEqual(entry["expected"]["nodes"], 6)
        self.assertEqual(entry["expected"]["cells"], 2)
        self.assertEqual(entry["expected"]["internal_interface_edges"], 1)

        obsolete_case_dir = REPOSITORY_ROOT / "grid" / "conversion" / "cases" / "two_quads"
        self.assertFalse(obsolete_case_dir.exists())
        self.assertFalse((REPOSITORY_ROOT / "grid" / "gmsh_to_cgns").exists())
        self.assertFalse((REPOSITORY_ROOT / "grid" / "gmsh_to_oneflow").exists())
        self.assertFalse((REPOSITORY_ROOT / "grid" / "cgns_to_oneflow").exists())

    def test_stl_ascii_basic_is_a_native_oneflow_case(self) -> None:
        catalog = load_catalog(REPOSITORY_ROOT / "grid" / "case_catalog.json")
        entry = next(case for case in catalog["cases"] if case["id"] == "stl_ascii_basic")
        self.assertEqual(entry["directory"], "grid/stl_import/ascii_basic/cpp")
        self.assertNotIn("runtimeConfig", entry)
        self.assertNotIn("execution", entry)
        case_dir = REPOSITORY_ROOT / entry["directory"]
        self.assertFalse((case_dir / "script" / "control.txt").exists())
        grid = (case_dir / "script" / "grid.txt").read_text(encoding="utf-8")
        dimension = (case_dir / "script" / "dimension.txt").read_text(encoding="utf-8")
        self.assertIn("int dimension = 2;", dimension)
        self.assertEqual(
            [path.name for path in sorted((case_dir / "script").glob("*.txt"))],
            ["dimension.txt", "grid.txt", "simu.txt"],
        )
        self.assertIn('string sourceGridType = "stl";', grid)
        self.assertIn('string targetGridType = "cgns";', grid)
        self.assertIn("Real stlWeldTolerance = 1.0e-10;", grid)
        self.assertTrue((case_dir / "input" / "model.stl").is_file())
        self.assertTrue((case_dir.parent / "README.md").is_file())
        self.assertTrue((case_dir.parent / "python" / "input" / "model.stl").is_file())
        self.assertTrue((case_dir.parent / "python" / "expected.json").is_file())
        self.assertEqual(validate_case(REPOSITORY_ROOT, entry), [])
    def test_stl_to_gmsh_case_family_has_no_legacy_flat_layout(self) -> None:
        case_root = REPOSITORY_ROOT / "grid" / "stl_to_gmsh"
        self.assertTrue((case_root / "README.md").is_file())
        for legacy_directory in ("input", "results", "script", "config", "pipeline", "python"):
            self.assertFalse((case_root / legacy_directory).exists(), legacy_directory)

        two_square_cells = case_root / "2d" / "two_square_cells"
        for implementation in ("cpp", "python"):
            case_dir = two_square_cells / implementation
            self.assertTrue((case_dir / "input" / "model.stl").is_file())
            self.assertTrue((case_dir / "script" / "layout.json").is_file())

    def test_stl_to_gmsh_cases_register_native_cpp_and_python_variants(self) -> None:
        catalog = load_catalog(REPOSITORY_ROOT / "grid" / "case_catalog.json")
        entries = {entry["id"]: entry for entry in catalog["cases"] if entry.get("group") == "STL to Gmsh"}
        self.assertEqual(set(entries), {
            "stl_to_gmsh_square",
            "stl_to_gmsh_square_two_triangles",
            "stl_to_gmsh_square_automatic_triangles",
            "stl_to_gmsh_square_2x2",
            "stl_to_gmsh_two_square_cells",
            "stl_to_gmsh_circle_external_flow",
            "stl_to_gmsh_python_square",
            "stl_to_gmsh_python_square_two_triangles",
            "stl_to_gmsh_python_square_automatic_triangles",
            "stl_to_gmsh_python_two_square_cells",
            "stl_to_gmsh_python_circle_external_flow",
            "stl_to_gmsh_python_channel_internal_flow",
            "stl_to_gmsh_cube_single_hex",
            "stl_to_gmsh_python_cube_single_hex",
            "stl_to_gmsh_cube_six_tetrahedra",
            "stl_to_gmsh_python_cube_six_tetrahedra",
            "stl_to_gmsh_cube_automatic_tetrahedra",
            "stl_to_gmsh_python_cube_automatic_tetrahedra",
            "stl_to_gmsh_sphere_automatic_tetrahedra",
            "stl_to_gmsh_python_sphere_automatic_tetrahedra",
            "stl_to_gmsh_prolate_spheroid_automatic_tetrahedra",
            "stl_to_gmsh_python_prolate_spheroid_automatic_tetrahedra",
            "stl_to_gmsh_python_cylinder_internal_flow_automatic_tetrahedra",
            "stl_to_gmsh_cylinder_internal_flow_automatic_tetrahedra",
            "stl_to_gmsh_cylinder_internal_flow_oblique_cut",
            "stl_to_gmsh_python_cylinder_internal_flow_oblique_cut",
            "stl_to_gmsh_pipe_bend_90_degree",
        })
        for case_id in (
            "stl_to_gmsh_square",
            "stl_to_gmsh_square_two_triangles",
            "stl_to_gmsh_square_automatic_triangles",
            "stl_to_gmsh_square_2x2",
            "stl_to_gmsh_two_square_cells",
            "stl_to_gmsh_circle_external_flow",
            "stl_to_gmsh_cube_single_hex",
            "stl_to_gmsh_cube_six_tetrahedra",
            "stl_to_gmsh_cube_automatic_tetrahedra",
            "stl_to_gmsh_prolate_spheroid_automatic_tetrahedra",
            "stl_to_gmsh_cylinder_internal_flow_automatic_tetrahedra",
        ):
            entry = entries[case_id]
            self.assertNotIn("execution", entry)
            self.assertTrue(entry["directory"].endswith("/cpp"))
            case_dir = REPOSITORY_ROOT / entry["directory"]
            self.assertFalse((case_dir / "script" / "control.txt").exists())
            self.assertEqual(validate_case(REPOSITORY_ROOT, entry), [])
            grid = (case_dir / "script" / "grid.txt").read_text(encoding="utf-8")
            self.assertIn('string sourceGridType = "stl";', grid)
            self.assertIn('string targetGridType = "gmsh";', grid)
            self.assertTrue((case_dir / "script" / "layout.json").is_file())
            self.assertTrue((case_dir / "input" / "model.stl").is_file())
        self.assertFalse((REPOSITORY_ROOT / "grid/stl_to_gmsh/python").exists())
        stl_case_root = REPOSITORY_ROOT / "grid/stl_to_gmsh"
        for legacy_directory in ("input", "results", "script", "config"):
            self.assertFalse((stl_case_root / legacy_directory).exists())
        for case_id in (
            "stl_to_gmsh_python_square",
            "stl_to_gmsh_python_two_square_cells",
            "stl_to_gmsh_python_circle_external_flow",
            "stl_to_gmsh_python_channel_internal_flow",
            "stl_to_gmsh_python_cube_single_hex",
            "stl_to_gmsh_python_cube_six_tetrahedra",
            "stl_to_gmsh_python_cube_automatic_tetrahedra",
            "stl_to_gmsh_python_prolate_spheroid_automatic_tetrahedra",
            "stl_to_gmsh_python_cylinder_internal_flow_automatic_tetrahedra",
        ):
            entry = entries[case_id]
            self.assertEqual(entry["execution"]["kind"], "python")
            self.assertTrue(entry["directory"].endswith("/python"))
            self.assertEqual(entry["execution"]["script"], "grid/stl_to_gmsh/run_python_case.py")
            self.assertEqual(validate_case(REPOSITORY_ROOT, entry), [])

    def test_native_cylinder_case_uses_shared_parametric_geometry_contract(self) -> None:
        catalog = load_catalog(REPOSITORY_ROOT / "grid" / "case_catalog.json")
        entry = next(
            case for case in catalog["cases"]
            if case["id"] == "stl_to_gmsh_cylinder_internal_flow_automatic_tetrahedra"
        )
        case_dir = REPOSITORY_ROOT / entry["directory"]
        layout = json.loads(
            (case_dir / "script" / "layout.json").read_text(encoding="utf-8")
        )
        manifest = json.loads(
            (case_dir / "input" / "model.geometry.json").read_text(encoding="utf-8")
        )
        self.assertEqual(layout["volume_mesh"], "cylinder_internal_flow")
        self.assertEqual(layout["mesh_dimension"], 3)
        self.assertEqual(layout["reference_geometry"]["type"], "cylinder_internal_flow")
        self.assertEqual(manifest["asset_id"], "geometry.cylinder.internal_flow")
        self.assertEqual(manifest["parameters"]["radius"], 1.0)
        self.assertEqual(manifest["parameters"]["length"], 4.0)
        self.assertEqual(
            set(manifest["patch_triangle_indices"]),
            {"wall", "inflow", "outflow"},
        )
        self.assertTrue((case_dir / "input" / "model.stl").is_file())
        self.assertNotIn("execution", entry)
        self.assertEqual(validate_case(REPOSITORY_ROOT, entry), [])

    def test_native_cylinder_fixture_checksum_and_patch_partition(self) -> None:
        attributes = (REPOSITORY_ROOT / ".gitattributes").read_text(encoding="utf-8")
        self.assertIn("*.stl -text", attributes.splitlines())
        case_dir = (
            REPOSITORY_ROOT / "grid/stl_to_gmsh/3d/cylinder_internal_flow"
            / "automatic_tetrahedra/cpp"
        )
        stl = case_dir / "input/model.stl"
        manifest = json.loads(
            stl.with_suffix(".geometry.json").read_text(encoding="utf-8")
        )
        self.assertEqual(
            hashlib.sha256(stl.read_bytes()).hexdigest(),
            manifest["stl_sha256"],
        )
        self.assertTrue(stl.read_text(encoding="ascii").startswith("solid "))
        self.assertEqual(manifest["surface"]["triangle_count"], 1632)
        patches = manifest["patch_triangle_indices"]
        self.assertEqual(len(patches["wall"]), 1536)
        self.assertEqual(len(patches["inflow"]), 48)
        self.assertEqual(len(patches["outflow"]), 48)
        indices = patches["wall"] + patches["inflow"] + patches["outflow"]
        self.assertEqual(sorted(indices), list(range(1632)))

    def test_stl_ascii_crlf_is_a_native_oneflow_case_with_crlf_fixture(self) -> None:
        catalog = load_catalog(REPOSITORY_ROOT / "grid" / "case_catalog.json")
        entry = next(
            case for case in catalog["cases"] if case["id"] == "stl_ascii_crlf"
        )
        self.assertEqual(entry["directory"], "grid/stl_import/ascii_crlf/cpp")
        self.assertNotIn("runtimeConfig", entry)
        self.assertNotIn("execution", entry)

        case_dir = REPOSITORY_ROOT / entry["directory"]
        self.assertFalse((case_dir / "script" / "control.txt").exists())
        fixture = (case_dir / "input" / "model.stl").read_bytes()
        self.assertIn(bytes((13, 10)), fixture)
        self.assertNotIn(bytes((10,)), fixture.replace(bytes((13, 10)), b""))
        self.assertIn(b"endsolid oneflow_ascii_square\r\n", fixture)
        self.assertEqual(
            [path.name for path in sorted((case_dir / "script").glob("*.txt"))],
            ["dimension.txt", "grid.txt", "simu.txt"],
        )
        self.assertEqual(validate_case(REPOSITORY_ROOT, entry), [])

    def test_stl_catalog_uses_language_neutral_case_names_and_python_namespace(self) -> None:
        catalog = load_catalog(REPOSITORY_ROOT / "grid" / "case_catalog.json")
        stl_cases = [case for case in catalog["cases"] if case.get("group") == "STL"]
        self.assertEqual(len(stl_cases), 16)
        stl_root = REPOSITORY_ROOT / "grid" / "stl_import"
        for legacy in ("input", "results", "script", "python"):
            self.assertFalse((stl_root / legacy).exists(), legacy)

        for entry in stl_cases:
            parts = entry["directory"].split("/")
            execution = entry.get("execution")
            if isinstance(execution, dict) and execution.get("kind") == "python":
                self.assertEqual(parts[:2], ["grid", "stl_import"])
                self.assertEqual(parts[-1], "python")
                self.assertEqual(len(parts), 4)
                self.assertEqual(entry["id"], f"stl_python_{parts[-2]}")
            else:
                self.assertEqual(parts[:2], ["grid", "stl_import"])
                self.assertEqual(parts[-1], "cpp")
                self.assertEqual(len(parts), 4)
                self.assertEqual(entry["id"], f"stl_{parts[-2]}")

    def test_native_stl_binary_cases_have_standard_binary_record_layout(self) -> None:
        catalog = load_catalog(REPOSITORY_ROOT / "grid" / "case_catalog.json")
        cases = {
            case["id"]: case
            for case in catalog["cases"]
            if case["id"] in {"stl_binary_basic", "stl_binary_header_solid"}
        }
        self.assertEqual(set(cases), {"stl_binary_basic", "stl_binary_header_solid"})
        for case_id, entry in cases.items():
            self.assertNotIn("execution", entry)
            self.assertNotIn("runtimeConfig", entry)
            fixture = (
                REPOSITORY_ROOT / entry["directory"] / "input" / "model.stl"
            ).read_bytes()
            self.assertEqual(len(fixture), 84 + 50 * 2)
            self.assertEqual(int.from_bytes(fixture[80:84], "little"), 2)
        solid_header = (
            REPOSITORY_ROOT / cases["stl_binary_header_solid"]["directory"]
            / "input" / "model.stl"
        ).read_bytes()
        self.assertTrue(solid_header.startswith(b"solid"))

    def test_native_stl_surface_cases_cover_non_planarity_and_welding_tolerance(self) -> None:
        catalog = load_catalog(REPOSITORY_ROOT / "grid" / "case_catalog.json")
        cases = {
            case["id"]: case
            for case in catalog["cases"]
            if case["id"] in {"stl_non_planar_surface", "stl_weld_tolerance"}
        }
        self.assertEqual(set(cases), {"stl_non_planar_surface", "stl_weld_tolerance"})
        for entry in cases.values():
            self.assertNotIn("execution", entry)
            self.assertNotIn("runtimeConfig", entry)
            case_root = REPOSITORY_ROOT / entry["directory"]
            self.assertTrue((case_root.parent / "README.md").is_file())
            for relative in (
                "script/dimension.txt", "script/simu.txt",
                "script/grid.txt", "input/model.stl",
            ):
                self.assertTrue((case_root / relative).is_file(), relative)
        non_planar = (
            REPOSITORY_ROOT / cases["stl_non_planar_surface"]["directory"]
            / "input" / "model.stl"
        ).read_text(encoding="ascii")
        self.assertIn("vertex 0 1 0.25", non_planar)
        weld_fixture = (
            REPOSITORY_ROOT / cases["stl_weld_tolerance"]["directory"]
            / "input" / "model.stl"
        ).read_text(encoding="ascii")
        self.assertIn("vertex 1.000000000051 1 0", weld_fixture)

    def test_turbulent_flat_plate_uses_requested_bounds_and_tanh_distribution(self) -> None:
        case_dir = REPOSITORY_ROOT / "grid" / "mesh-layout" / "2d" / "turbulent_flat_plate"
        with (case_dir / "case.json").open("r", encoding="utf-8") as stream:
            profile = json.load(stream)

        runtime = profile["parameters"]["runtime"]
        expected = profile["parameters"]["expected"]
        layout = (case_dir / runtime["layout"]).read_text(encoding="utf-8")
        self.assertEqual(runtime["dimension"], 2)
        self.assertIn("Point(1)= {-0.35,0,0}", layout)
        self.assertIn("Point(2)= {0,0,0}", layout)
        self.assertIn("Point(3)= {2,0,0}", layout)
        self.assertIn("Point(4)= {-0.35,1,0}", layout)
        self.assertIn("Point(5)= {0,1,0}", layout)
        self.assertIn("Point(6)= {2,1,0}", layout)
        self.assertIn("Dim(3)={81}", layout)
        self.assertIn("Dim(4)={81}", layout)
        self.assertIn("Dim(5)={81}", layout)
        expected_distributions = [
            "Ds(1)={tanh,-1,0.002}",
            "Ds(2)={tanh,0.002,-1}",
            "Ds(3)={tanh,1.0e-6,-1}",
            "Ds(4)={tanh,1.0e-6,-1}",
            "Ds(5)={tanh,1.0e-6,-1}",
            "Ds(6)={tanh,-1,0.002}",
            "Ds(7)={tanh,0.002,-1}",
        ]
        for distribution in expected_distributions:
            self.assertIn(distribution, layout)
        self.assertEqual(expected["nodes"], 5913)
        self.assertEqual(expected["cells"], 5760)
        self.assertEqual(expected["faces"], 11672)
        self.assertEqual(expected["boundaryFaces"], 304)
        self.assertEqual(expected["internalFaces"], 11368)
        self.assertEqual(expected["boundingBox"]["min"], [-0.35, 0.0, 0.0])
        self.assertEqual(expected["boundingBox"]["max"], [2.0, 1.0, 0.0])

    def test_two_blocks_shared_edge_endpoints_only_profile_describes_minimum_topology(self) -> None:
        case_dir = REPOSITORY_ROOT / "grid" / "mesh-layout" / "2d" / "two_blocks_shared_edge_endpoints_only"
        profile = json.loads((case_dir / "case.json").read_text(encoding="utf-8"))
        expected = profile["parameters"]["expected"]
        self.assertEqual(
            {key: expected[key] for key in ("zones", "nodes", "cells", "faces", "boundaryFaces", "internalFaces")},
            {"zones": 2, "nodes": 6, "cells": 2, "faces": 7, "boundaryFaces": 6, "internalFaces": 1},
        )
        layout = (case_dir / "grid" / "twoBlocksSharedEdgeEndpointsOnlyLayout.txt").read_text(encoding="utf-8")
        self.assertEqual(layout.count("Dim("), 7)
        self.assertEqual(layout.count("= {2}"), 7)

    def test_two_blocks_shared_edge_one_interior_node_profile_carries_grid_options(self) -> None:
        case_dir = REPOSITORY_ROOT / "grid" / "mesh-layout" / "2d" / "two_blocks_shared_edge_one_interior_node"
        with (case_dir / "case.json").open("r", encoding="utf-8") as stream:
            profile = json.load(stream)

        runtime = profile["parameters"]["runtime"]
        options = runtime["gridOptions"]
        self.assertEqual(runtime["dimension"], 2)
        self.assertEqual(runtime["objective"], 0)
        self.assertEqual(runtime["generationId"], 4)
        self.assertEqual(options["topology"], "unstructured")
        self.assertTrue(options["multiBlock"])
        self.assertEqual(options["fileFormat"], 1)
        self.assertEqual(options["ioRoots"]["project"], ".")
        self.assertEqual(options["ioRoots"]["meshOutput"], "case")
        self.assertEqual(options["translation"], [0.0, 0.0, 0.0])
        self.assertEqual(options["axisDirection"], 0)
        self.assertEqual(options["sourceFormat"], "plot3d")
        self.assertEqual(options["targetFormat"], "oneflow")
        self.assertEqual(options["parallelStrategy"], 1)
        self.assertEqual(profile["parameters"]["expected"]["cells"], 8)
        self.assertEqual(profile["parameters"]["expected"]["internalFaces"], 9)

        grid_parameters = (case_dir / "script" / "grid.txt").read_text(encoding="utf-8")
        self.assertIn('string topoType = "u";', grid_parameters)
        self.assertIn("int multiBlock = 1;", grid_parameters)
        self.assertIn("int fileFormat = 1;", grid_parameters)
        self.assertIn("Real gridScale = 1.0;", grid_parameters)
        self.assertIn("Real gridTrans[3] = 0, 0, 0;", grid_parameters)
        self.assertIn('string gridLayoutFileName = "grid/twoBlocksSharedEdgeOneInteriorNodeLayout.txt";', grid_parameters)
        self.assertIn('string sourceGridType = "plot3d";', grid_parameters)
        self.assertIn('string targetGridType = "oneflow";', grid_parameters)
        self.assertIn("int parallelStrategy = 1;", grid_parameters)
        self.assertIn('string projectRoot = ".";', grid_parameters)
        self.assertIn('string inputRoot = "case";', grid_parameters)
        self.assertIn('string outputRoot = "case";', grid_parameters)
        self.assertIn('string meshInputRoot = "case";', grid_parameters)
        self.assertIn('string flowInputRoot = "case";', grid_parameters)
        self.assertIn('string meshOutputRoot = "case";', grid_parameters)
        self.assertIn('string flowOutputRoot = "case";', grid_parameters)

        layout = (case_dir / runtime["layout"]).read_text(encoding="utf-8")
        self.assertIn("Boundary(2) = {-1}", layout)
        self.assertIn("Add L2F 1 2 2", layout)
        self.assertIn("Add L2F 2 4 2", layout)
        self.assertIn("Add F2B 1 1 1", layout)
        self.assertIn("Add F2B 2 1 2", layout)
        self.assertNotIn("Add F2B 1 2 2", layout)
        self.assertNotIn("Add F2B 2 4 2", layout)

    def test_invalid_grid_options_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            case_dir = Path(directory)
            (case_dir / "script").mkdir()
            (case_dir / "script" / "grid.txt").write_text("", encoding="utf-8")
            entry = {"id": "sample", "layout": "grid/layout.txt"}
            profile = {
                "schema": "oneflow.case",
                "version": 1,
                "parameters": {
                    "runtime": {
                        "dimension": 2,
                        "layout": "grid/layout.txt",
                        "parameterFiles": ["script/grid.txt"],
                        "gridOptions": {
                            "topology": "volume",
                            "multiBlock": "yes",
                            "fileFormat": "1",
                            "scale": "1.0",
                            "translation": [0.0, 1.0],
                            "ioRoots": {"output": "case"},
                        },
                    }
                },
            }
            errors = validate_profile(case_dir, entry, profile)
            self.assertTrue(any("gridOptions.topology" in error for error in errors))
            self.assertTrue(any("gridOptions.multiBlock" in error for error in errors))
            self.assertTrue(any("gridOptions.fileFormat" in error for error in errors))
            self.assertTrue(any("gridOptions.scale" in error for error in errors))
            self.assertTrue(any("gridOptions.translation" in error for error in errors))
            self.assertTrue(any("gridOptions.ioRoots.meshOutput" in error for error in errors))

    def test_native_case_without_control_manifest_uses_parameter_file_discovery(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            script_dir = root / "grid" / "sample" / "script"
            script_dir.mkdir(parents=True)
            (script_dir / "simu.txt").write_text(
                'string simutask = "Grid";\n', encoding="utf-8"
            )
            (script_dir / "grid.txt").write_text(
                "int gridObj = 1;\n", encoding="utf-8"
            )
            entry = {"id": "sample", "directory": "grid/sample"}
            self.assertEqual(validate_case(root, entry), [])

    def test_manifest_free_case_requires_parameter_files(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            (Path(directory) / "grid" / "sample" / "script").mkdir(parents=True)
            entry = {"id": "sample", "directory": "grid/sample"}
            errors = validate_case(Path(directory), entry)
            self.assertTrue(
                any("no script/*.txt parameter files found" in error for error in errors)
            )

    def test_duplicate_case_ids_are_rejected(self) -> None:
        catalog = {
            "schema": "oneflow.case-catalog",
            "version": 1,
            "cases": [
                {"id": "same", "directory": "a", "runtimeConfig": "script/control.txt"},
                {"id": "same", "directory": "b", "runtimeConfig": "script/control.txt"},
            ],
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "catalog.json"
            path.write_text(json.dumps(catalog), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Duplicate case id"):
                load_catalog(path)

    def test_pipe_bend_declares_incremental_lifecycle_not_final_completion(self) -> None:
        catalog = load_catalog(REPOSITORY_ROOT / "grid" / "case_catalog.json")
        entry = next(
            case for case in catalog["cases"]
            if case["id"] == "stl_to_gmsh_pipe_bend_90_degree"
        )
        self.assertEqual(
            entry["lifecycle_contract"],
            "stl_to_gmsh/3d/pipe_bend/90_degree/automatic_tetrahedra/case-lifecycle.json",
        )
        contract_path = REPOSITORY_ROOT / "grid" / entry["lifecycle_contract"]
        contract = json.loads(contract_path.read_text(encoding="utf-8"))
        self.assertEqual(
            validate_lifecycle_contract(
                contract, expected_case_id="stl_to_gmsh_pipe_bend_90_degree"
            ),
            [],
        )
        stages = {stage["id"]: stage for stage in contract["stages"]}
        self.assertEqual(
            contract["verification_scope"],
            "geometry_regions_gmsh_volume_mesh_cgns_conversion_oneflow_ingestion_and_cfd_smoke",
        )
        self.assertEqual(
            contract["verification_stages"],
            [
                "geometry",
                "surface_regions",
                "volume_mesh",
                "msh_to_cgns",
                "oneflow_ingestion",
                "cfd_smoke",
            ],
        )
        self.assertEqual(stages["msh_to_cgns"]["implementation_status"], "implemented")
        self.assertEqual(stages["msh_to_cgns"]["verification_status"], "not_run")
        self.assertEqual(stages["oneflow_ingestion"]["implementation_status"], "implemented")
        self.assertEqual(stages["oneflow_ingestion"]["verification_status"], "not_run")
        self.assertEqual(stages["geometry"]["implementation_status"], "implemented")
        self.assertEqual(stages["geometry"]["verification_status"], "not_run")
        geometry_implementations = {
            item["id"]: item["status"]
            for item in stages["geometry"]["implementations"]
        }
        self.assertEqual(
            geometry_implementations,
            {"python_geometry": "implemented", "native_cpp_geometry": "implemented"},
        )
        region_implementations = {
            item["id"]: item["status"]
            for item in stages["surface_regions"]["implementations"]
        }
        self.assertEqual(
            region_implementations,
            {"python_region_export": "implemented", "native_cpp_region_export": "implemented"},
        )
        self.assertEqual(stages["cfd_smoke"]["implementation_status"], "implemented")
        self.assertEqual(stages["cfd_smoke"]["verification_status"], "not_run")
        self.assertEqual(stages["cfd_validation"]["implementation_status"], "planned")
        acceptance_references = {
            check.get("reference")
            for stage_id in ("geometry", "surface_regions")
            for check in stages[stage_id]["acceptance"]
        }
        self.assertEqual(
            acceptance_references,
            {
                "grid.tests.python.stl_to_gmsh.test_stl_to_gmsh_cases."
                "StlToGmshCaseTest.test_pipe_bend_geometry_contract_emits_verified_region_assets",
                "grid.tests.python.stl_to_gmsh.test_pipe_bend_geometry."
                "PipeBendGeometryTest.test_cpp_generator_matches_python_geometry",
            },
        )

    def test_lifecycle_contract_allows_incremental_stages_but_rejects_cycles(self) -> None:
        contract = {
            "schema": "oneflow.case-lifecycle",
            "version": 1,
            "case_id": "sample",
            "contract_version": 1,
            "verification_scope": "geometry_only",
            "verification_stages": ["geometry"],
            "stages": [
                {
                    "id": "geometry",
                    "type": "geometry",
                    "implementation_status": "implemented",
                    "verification_status": "not_run",
                    "depends_on": [],
                    "implementations": [
                        {"id": "python", "kind": "python", "status": "implemented"}
                    ],
                    "inputs": [],
                    "outputs": [{"id": "surface", "format": "stl"}],
                    "acceptance": [
                        {"id": "geometry.topology", "kind": "contract", "description": "Closed surface"}
                    ],
                },
                {
                    "id": "cfd",
                    "type": "cfd_smoke",
                    "implementation_status": "planned",
                    "verification_status": "not_run",
                    "depends_on": ["geometry"],
                    "implementations": [],
                    "inputs": [{"id": "surface", "format": "stl"}],
                    "outputs": [],
                    "acceptance": [],
                },
            ],
        }
        self.assertEqual(validate_lifecycle_contract(contract), [])
        contract["stages"][0]["depends_on"] = ["cfd"]
        self.assertTrue(any(
            "cycle" in error for error in validate_lifecycle_contract(contract)
        ))

    def test_passed_lifecycle_stage_requires_implemented_dependencies(self) -> None:
        contract = {
            "schema": "oneflow.case-lifecycle",
            "version": 1,
            "case_id": "sample",
            "contract_version": 1,
            "verification_scope": "mesh",
            "verification_stages": ["geometry", "mesh"],
            "stages": [
                {
                    "id": "geometry",
                    "type": "geometry",
                    "implementation_status": "implemented",
                    "verification_status": "not_run",
                    "depends_on": [],
                    "implementations": [],
                    "inputs": [],
                    "outputs": [],
                    "acceptance": [],
                },
                {
                    "id": "mesh",
                    "type": "volume_discretization",
                    "implementation_status": "implemented",
                    "verification_status": "passed",
                    "depends_on": ["geometry"],
                    "implementations": [],
                    "inputs": [],
                    "outputs": [],
                    "acceptance": [
                        {"id": "mesh.check", "kind": "integration", "description": "Mesh exists"}
                    ],
                },
            ],
        }
        errors = validate_lifecycle_contract(contract)
        self.assertTrue(any("cannot be passed" in error for error in errors))

    def test_non_object_catalog_root_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "catalog.json"
            path.write_text("[]", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "root must be an object"):
                load_catalog(path)

    def test_parameter_include_cannot_escape_case_directory(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            case_dir = root / "grid" / "sample"
            script_dir = case_dir / "script"
            script_dir.mkdir(parents=True)
            (script_dir / "control.txt").write_text(
                '"../../outside.txt";\n', encoding="utf-8"
            )
            (root / "outside.txt").write_text("outside", encoding="utf-8")
            entry = {
                "id": "sample",
                "directory": "grid/sample",
                "runtimeConfig": "script/control.txt",
                "layout": None,
                "profile": None,
            }
            errors = validate_case(root, entry)
            self.assertTrue(any("escapes the case directory" in error for error in errors))

    def test_profile_layout_must_match_catalog(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            case_dir = Path(directory)
            (case_dir / "script").mkdir()
            (case_dir / "script" / "grid.txt").write_text("", encoding="utf-8")
            entry = {
                "id": "sample",
                "layout": "grid/catalog-layout.txt",
            }
            profile = {
                "schema": "oneflow.case",
                "version": 1,
                "parameters": {
                    "runtime": {
                        "dimension": 3,
                        "layout": "grid/profile-layout.txt",
                        "parameterFiles": ["script/grid.txt"],
                    }
                },
            }
            errors = validate_profile(case_dir, entry, profile)
            self.assertTrue(any("does not match catalog layout" in error for error in errors))
            self.assertFalse(any("parameter file not found" in error for error in errors))

    def test_profile_root_must_be_an_object(self) -> None:
        errors = validate_profile(
            Path("."),
            {"id": "sample", "layout": None},
            [],
        )
        self.assertEqual(errors, ["sample: case profile root must be an object"])


if __name__ == "__main__":
    unittest.main()
