from __future__ import annotations

import json
import hashlib
import os
import shutil
import subprocess
import tempfile
import unittest
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "grid/stl_to_gmsh/3d/pipe_bend/90_degree/geometry/python"))

from mesh_cases import load_catalog, validate_case  # noqa: E402
from grid.stl_to_gmsh.converter import convert, inspect_msh4, read_stl  # noqa: E402
from grid.stl_to_gmsh.run_cases import NATIVE_CASES, PYTHON_CASES, REFINEMENT_LEVELS, SPHERE_REFINEMENT_CASES  # noqa: E402
from generate_pipe_bend import generate as generate_pipe_bend  # noqa: E402

def _manual_gmsh():
    """Return Gmsh only when a developer explicitly opts into integration tests."""
    if os.environ.get("ONEFLOW_RUN_GMSH_TESTS") != "1":
        return None
    return shutil.which("gmsh") or shutil.which("gmsh.exe")


EXPECTED_IDS = {
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
}


class StlToGmshCaseTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        catalog = load_catalog(ROOT / "grid/case_catalog.json")
        cls.entries = {
            entry["id"]: entry
            for entry in catalog["cases"]
            if entry.get("group") == "STL to Gmsh"
        }

    def test_staging_paths_mirror_case_implementation_layout(self):
        self.assertEqual(NATIVE_CASES, {
            "2d/square/single_quad/cpp": ("2d/square/single_quad/cpp", "square.msh"),
            "2d/square/two_triangles/cpp": (
                "2d/square/two_triangles/cpp", "square_two_triangles.msh"
            ),
            "2d/square/automatic_triangles/cpp": (
                "2d/square/automatic_triangles/cpp", "square_automatic_triangles.msh"
            ),
            "2d/square_2x2/cpp": ("2d/square_2x2/cpp", "square_2x2.msh"),
            "2d/two_square_cells/cpp": ("2d/two_square_cells/cpp", "two_square_cells.msh"),
            "2d/circle_external_flow/cpp": (
                "2d/circle_external_flow/cpp", "circle_external_flow.msh"
            ),
            "3d/cube/single_hex/cpp": ("3d/cube/single_hex/cpp", "cube_single_hex.msh"),
            "3d/cube/six_tetrahedra/cpp": ("3d/cube/six_tetrahedra/cpp", "cube_six_tetrahedra.msh"),
            "3d/cube/automatic_tetrahedra/cpp": ("3d/cube/automatic_tetrahedra/cpp", "cube_automatic_tetrahedra.msh"),
            "3d/prolate_spheroid/automatic_tetrahedra/cpp": ("3d/prolate_spheroid/automatic_tetrahedra/cpp", "prolate_spheroid_automatic_tetrahedra.msh"),
    "3d/cylinder_internal_flow/automatic_tetrahedra/cpp": ("3d/cylinder_internal_flow/automatic_tetrahedra/cpp", "cylinder_internal_flow_automatic_tetrahedra.msh"),
            "3d/cylinder_internal_flow/automatic_tetrahedra/cpp": ("3d/cylinder_internal_flow/automatic_tetrahedra/cpp", "cylinder_internal_flow_automatic_tetrahedra.msh"),
        })
        self.assertEqual(PYTHON_CASES, {
            "2d/square/single_quad/python": "2d/square/single_quad/python",
            "2d/square/two_triangles/python": "2d/square/two_triangles/python",
            "2d/square/automatic_triangles/python": "2d/square/automatic_triangles/python",
            "2d/two_square_cells/python": "2d/two_square_cells/python",
            "2d/circle_external_flow/python": "2d/circle_external_flow/python",
            "2d/channel_internal_flow/python": "2d/channel_internal_flow/python",
            "3d/cube/single_hex/python": "3d/cube/single_hex/python",
            "3d/cube/six_tetrahedra/python": "3d/cube/six_tetrahedra/python",
            "3d/cube/automatic_tetrahedra/python": "3d/cube/automatic_tetrahedra/python",
            "3d/prolate_spheroid/automatic_tetrahedra/python": "3d/prolate_spheroid/automatic_tetrahedra/python",
            "3d/cylinder_internal_flow/automatic_tetrahedra/python": "3d/cylinder_internal_flow/automatic_tetrahedra/python",
        })

    def test_cube_cases_are_registered_and_staged_for_both_implementations(self):
        native = self.entries["stl_to_gmsh_cube_single_hex"]
        python = self.entries["stl_to_gmsh_python_cube_single_hex"]
        self.assertEqual(native["directory"], "grid/stl_to_gmsh/3d/cube/single_hex/cpp")
        self.assertNotIn("execution", native)
        self.assertEqual(python["directory"], "grid/stl_to_gmsh/3d/cube/single_hex/python")
        self.assertEqual(python["execution"]["kind"], "python")
        self.assertEqual(validate_case(ROOT, native), [])
        self.assertEqual(validate_case(ROOT, python), [])
        self.assertEqual(
            (ROOT / native["directory"] / "input/model.stl").read_bytes(),
            (ROOT / python["directory"] / "input/model.stl").read_bytes(),
        )
        self.assertEqual(
            json.loads((ROOT / native["directory"] / "script/layout.json").read_text()),
            json.loads((ROOT / python["directory"] / "script/layout.json").read_text()),
        )
        self.assertEqual(
            json.loads((ROOT / native["directory"] / "script/layout.json").read_text())[
                "boundary_groups"
            ],
            {
                "xmin": "inflow",
                "xmax": "outflow",
                "ymin": "wall",
                "ymax": "farfield",
                "zmin": "symmetry",
                "zmax": "symmetry",
            },
        )
        self.assertIn("3d/cube/single_hex/cpp", NATIVE_CASES)
        self.assertIn("3d/cube/single_hex/python", PYTHON_CASES)
        self.assertIn("3d/cube/six_tetrahedra/cpp", NATIVE_CASES)
        self.assertIn("3d/cube/six_tetrahedra/python", PYTHON_CASES)
        self.assertIn("3d/cube/automatic_tetrahedra/cpp", NATIVE_CASES)
        self.assertIn("3d/cube/automatic_tetrahedra/python", PYTHON_CASES)
        for entry in (native, python):
            self.assertTrue((ROOT / entry["directory"] / "input/model.stl").is_file())
            self.assertTrue((ROOT / entry["directory"] / "script/layout.json").is_file())


    def test_cube_tetrahedral_cases_are_registered_in_native_and_python_pairs(self):
        for suffix in ("six_tetrahedra", "automatic_tetrahedra"):
            native = self.entries[f"stl_to_gmsh_cube_{suffix}"]
            python = self.entries[f"stl_to_gmsh_python_cube_{suffix}"]
            self.assertEqual(native["directory"], f"grid/stl_to_gmsh/3d/cube/{suffix}/cpp")
            self.assertEqual(python["directory"], f"grid/stl_to_gmsh/3d/cube/{suffix}/python")
            self.assertNotIn("execution", native)
            self.assertEqual(python["execution"]["kind"], "python")
            self.assertEqual(validate_case(ROOT, native), [])
            self.assertEqual(validate_case(ROOT, python), [])
            native_dir, python_dir = ROOT / native["directory"], ROOT / python["directory"]
            self.assertEqual(
                (native_dir / "input/model.stl").read_bytes(),
                (python_dir / "input/model.stl").read_bytes(),
            )
            native_layout = json.loads((native_dir / "script/layout.json").read_text())
            python_layout = json.loads((python_dir / "script/layout.json").read_text())
            self.assertEqual(native_layout, python_layout)
            self.assertEqual(native_layout["mesh_dimension"], 3)
            expected = json.loads((python_dir / "expected.json").read_text())
            if suffix == "six_tetrahedra":
                self.assertEqual(expected["mesh"]["nodes"], 8)
                self.assertEqual(expected["mesh"]["triangle_elements"], 12)
                self.assertEqual(expected["mesh"]["tetrahedron_elements"], 6)
            else:
                self.assertEqual(expected["mesh"]["hexahedron_elements"], 0)
                self.assertGreater(expected["mesh_minimum"]["tetrahedron_elements"], 0)
            self.assertFalse((native_dir / "results").exists())
            self.assertFalse((python_dir / "results").exists())

    def test_unit_sphere_cases_are_registered_as_a_native_python_pair(self):
        native = self.entries["stl_to_gmsh_sphere_automatic_tetrahedra"]
        python = self.entries["stl_to_gmsh_python_sphere_automatic_tetrahedra"]
        self.assertEqual(
            native["directory"],
            "grid/stl_to_gmsh/3d/sphere/automatic_tetrahedra/cpp",
        )
        self.assertEqual(
            python["directory"],
            "grid/stl_to_gmsh/3d/sphere/automatic_tetrahedra/python",
        )
        self.assertNotIn("execution", native)
        self.assertEqual(python["execution"]["kind"], "python")
        self.assertEqual(validate_case(ROOT, native), [])
        self.assertEqual(validate_case(ROOT, python), [])
        native_dir, python_dir = ROOT / native["directory"], ROOT / python["directory"]
        self.assertEqual(
            (native_dir / "input/model.stl").read_bytes(),
            (python_dir / "input/model.stl").read_bytes(),
        )
        native_layout = json.loads(
            (native_dir / "script/layout.json").read_text(encoding="utf-8")
        )
        python_layout = json.loads(
            (python_dir / "script/layout.json").read_text(encoding="utf-8")
        )
        self.assertEqual(native_layout, python_layout)
        self.assertEqual(native_layout["volume_mesh"], "closed_surface_automatic_tetrahedra")
        self.assertEqual(native_layout["reference_geometry"]["type"], "unit_sphere")
        self.assertEqual(native_layout["reference_geometry"]["center"], [0.0, 0.0, 0.0])
        self.assertEqual(native_layout["reference_geometry"]["radius"], 1.0)
        self.assertEqual(native_layout["surface_name"], "fluid_domain")
        self.assertEqual(native_layout["boundary_name"], "wall")
        self.assertEqual(
            len(read_stl(native_dir / "input/model.stl")),
            1280,
        )
        self.assertFalse((native_dir / "results").exists())
        self.assertFalse((python_dir / "results").exists())
        self.assertEqual(REFINEMENT_LEVELS, {"coarse": 0.40, "medium": 0.20, "fine": 0.10})
        self.assertEqual(SPHERE_REFINEMENT_CASES, {
            "cpp": "3d/sphere/automatic_tetrahedra/cpp",
            "python": "3d/sphere/automatic_tetrahedra/python",
        })

    def test_prolate_spheroid_native_and_python_cases_share_one_geometry_contract(self):
        native = self.entries["stl_to_gmsh_prolate_spheroid_automatic_tetrahedra"]
        python = self.entries["stl_to_gmsh_python_prolate_spheroid_automatic_tetrahedra"]
        self.assertEqual(
            native["directory"],
            "grid/stl_to_gmsh/3d/prolate_spheroid/automatic_tetrahedra/cpp",
        )
        self.assertEqual(
            python["directory"],
            "grid/stl_to_gmsh/3d/prolate_spheroid/automatic_tetrahedra/python",
        )
        self.assertNotIn("execution", native)
        self.assertEqual(python["execution"]["kind"], "python")
        self.assertEqual(native["mesh_validation"], "cfd_smoke")
        self.assertEqual(python["mesh_validation"], "cfd_smoke")
        self.assertEqual(validate_case(ROOT, native), [])
        self.assertEqual(validate_case(ROOT, python), [])
        native_dir, python_dir = ROOT / native["directory"], ROOT / python["directory"]
        self.assertEqual(
            (native_dir / "input/model.stl").read_bytes(),
            (python_dir / "input/model.stl").read_bytes(),
        )
        self.assertEqual(
            json.loads((native_dir / "script/layout.json").read_text(encoding="utf-8")),
            json.loads((python_dir / "script/layout.json").read_text(encoding="utf-8")),
        )
        self.assertEqual(len(read_stl(native_dir / "input/model.stl")), 1280)
        layout = json.loads((native_dir / "script/layout.json").read_text(encoding="utf-8"))
        self.assertEqual(layout["reference_geometry"]["type"], "prolate_spheroid")
        self.assertEqual(layout["reference_geometry"]["semi_axes"], [3.0, 0.5, 0.5])
        self.assertEqual(layout["reference_geometry"]["aspect_ratio"], 6.0)
        self.assertEqual(layout["reference_geometry"]["volume"], 3.141592653589793)
        grid = (native_dir / "script/grid.txt").read_text(encoding="utf-8")
        self.assertIn('string targetGridFileName = "results/prolate_spheroid_automatic_tetrahedra.msh";', grid)
        self.assertIn("3d/prolate_spheroid/automatic_tetrahedra/cpp", NATIVE_CASES)
        self.assertIn("3d/prolate_spheroid/automatic_tetrahedra/python", PYTHON_CASES)
        self.assertFalse((native_dir / "results").exists())
        self.assertFalse((python_dir / "results").exists())

    def test_internal_flow_cylinder_is_a_bundled_region_aware_reference_case(self):
        case = self.entries["stl_to_gmsh_python_cylinder_internal_flow_automatic_tetrahedra"]
        self.assertEqual(
            case["directory"],
            "grid/stl_to_gmsh/3d/cylinder_internal_flow/automatic_tetrahedra/python",
        )
        self.assertEqual(case["execution"]["kind"], "python")
        self.assertEqual(validate_case(ROOT, case), [])
        case_dir = ROOT / case["directory"]
        stl = case_dir / "input/model.stl"
        manifest = json.loads(
            (case_dir / "input/model.geometry.json").read_text(encoding="utf-8")
        )
        layout = json.loads((case_dir / "script/layout.json").read_text(encoding="utf-8"))
        expected = json.loads((case_dir / "expected.json").read_text(encoding="utf-8"))
        self.assertEqual(manifest["asset_id"], "geometry.cylinder.internal_flow")
        self.assertEqual(manifest["geometry"], "cylinder_internal_flow")
        self.assertEqual(manifest["stl_sha256"], __import__("hashlib").sha256(stl.read_bytes()).hexdigest())
        self.assertTrue(manifest["surface"]["watertight"])
        self.assertEqual(len(read_stl(stl)), 1632)
        self.assertEqual(layout["volume_mesh"], "cylinder_internal_flow")
        self.assertEqual(layout["mesh_dimension"], 3)
        self.assertEqual(layout["reference_geometry"]["radius"], 1.0)
        self.assertEqual(layout["reference_geometry"]["length"], 4.0)
        self.assertEqual(
            {(item["dimension"], item["name"]) for item in expected["physical_groups"]},
            {(2, "wall"), (2, "inflow"), (2, "outflow"), (3, "fluid_domain")},
        )
        self.assertFalse((case_dir / "results").exists())
        self.assertIn(
            "3d/cylinder_internal_flow/automatic_tetrahedra/python",
            PYTHON_CASES,
        )

    def test_oblique_cut_cylinder_cases_are_registered_for_both_implementations(self):
        native_id = "stl_to_gmsh_cylinder_internal_flow_oblique_cut"
        python_id = "stl_to_gmsh_python_cylinder_internal_flow_oblique_cut"
        native = self.entries[native_id]
        python = self.entries[python_id]

        self.assertEqual(
            native["directory"],
            "grid/stl_to_gmsh/3d/cylinder_internal_flow/oblique_cut/cpp",
        )
        self.assertNotIn("execution", native)
        self.assertEqual(
            python["directory"],
            "grid/stl_to_gmsh/3d/cylinder_internal_flow/oblique_cut/python",
        )
        self.assertEqual(
            python["execution"]["script"],
            "grid/stl_to_gmsh/3d/cylinder_internal_flow/oblique_cut/python/run_cylinder_oblique_cut.py",
        )
        self.assertEqual(python["execution"]["kind"], "python")
        self.assertEqual(validate_case(ROOT, native), [])
        self.assertEqual(validate_case(ROOT, python), [])
        native_dir = ROOT / native["directory"]
        self.assertEqual(
            (native_dir / "script/simu.txt").read_text(encoding="utf-8").splitlines(),
            ['string simutask = "Grid";'],
        )
        self.assertTrue((native_dir / "script/grid.txt").is_file())
        self.assertEqual(
            (native_dir / "script/dimension.txt").read_text(encoding="utf-8").splitlines(),
            ["int dimension = 3;"],
        )
        for entry in (native, python):
            self.assertTrue((ROOT / entry["directory"]).is_dir())
            self.assertFalse((ROOT / entry["directory"] / "results").exists())

    def test_native_and_python_cases_are_registered_separately(self):
        self.assertEqual(set(self.entries), EXPECTED_IDS)
        for case_id in (
            "stl_to_gmsh_square",
            "stl_to_gmsh_square_two_triangles",
            "stl_to_gmsh_square_automatic_triangles",
            "stl_to_gmsh_square_2x2",
            "stl_to_gmsh_two_square_cells",
            "stl_to_gmsh_circle_external_flow",
        ):
            entry = self.entries[case_id]
            case = ROOT / entry["directory"]
            self.assertNotIn("execution", entry)
            self.assertTrue(entry["directory"].endswith("/cpp"))
            self.assertEqual(validate_case(ROOT, entry), [])
            self.assertEqual(
                sorted(path.name for path in (case / "script").glob("*.txt")),
                ["dimension.txt", "grid.txt", "simu.txt"],
            )
            grid = (case / "script/grid.txt").read_text(encoding="utf-8")
            self.assertIn('string sourceGridType = "stl";', grid)
            self.assertIn('string targetGridType = "gmsh";', grid)
            self.assertIn('string gmshExecutable = "gmsh";', grid)
            self.assertFalse((case / "script/control.txt").exists())
        for case_id in (
            "stl_to_gmsh_python_square",
            "stl_to_gmsh_python_square_two_triangles",
            "stl_to_gmsh_python_square_automatic_triangles",
            "stl_to_gmsh_python_two_square_cells",
            "stl_to_gmsh_python_circle_external_flow",
            "stl_to_gmsh_python_channel_internal_flow",
        ):
            entry = self.entries[case_id]
            self.assertEqual(entry["execution"]["kind"], "python")
            self.assertTrue(entry["directory"].endswith("/python"))
            self.assertEqual(entry["execution"]["script"], "grid/stl_to_gmsh/run_python_case.py")
            self.assertEqual(validate_case(ROOT, entry), [])

    def test_automatic_triangle_case_pair_is_registered_and_uses_mesh_size_only(self):
        native = self.entries["stl_to_gmsh_square_automatic_triangles"]
        python = self.entries["stl_to_gmsh_python_square_automatic_triangles"]
        self.assertEqual(validate_case(ROOT, native), [])
        self.assertEqual(validate_case(ROOT, python), [])
        native_dir, python_dir = ROOT / native["directory"], ROOT / python["directory"]
        self.assertEqual((native_dir / "input/model.stl").read_bytes(),
                         (python_dir / "input/model.stl").read_bytes())
        native_layout = json.loads((native_dir / "script/layout.json").read_text())
        python_layout = json.loads((python_dir / "script/layout.json").read_text())
        self.assertEqual(native_layout, python_layout)
        self.assertEqual(native_layout["mesh_size"], 0.25)
        self.assertNotIn("transfinite_cells", native_layout)
        self.assertNotIn("recombine_surface", native_layout)
        expected = json.loads((python_dir / "expected.json").read_text())
        self.assertEqual(expected["mesh"], {"quadrilateral_elements": 0})
        self.assertEqual(expected["mesh_minimum"], {
            "nodes": 5, "line_elements": 4, "triangle_elements": 1,
        })
        self.assertFalse((native_dir / "results").exists())
        self.assertFalse((python_dir / "results").exists())

    def test_native_circle_external_flow_uses_same_stl_and_domain_as_python_reference(self):
        native = ROOT / "grid/stl_to_gmsh/2d/circle_external_flow/cpp"
        reference = ROOT / "grid/stl_to_gmsh/2d/circle_external_flow/python"
        self.assertEqual(
            (native / "input/model.stl").read_bytes(),
            (reference / "input/model.stl").read_bytes(),
        )
        native_layout = json.loads(
            (native / "script/layout.json").read_text(encoding="utf-8")
        )
        python_layout = json.loads(
            (reference / "script/layout.json").read_text(encoding="utf-8")
        )
        self.assertEqual(native_layout, python_layout)
        self.assertEqual(native_layout["surface_name"], "fluid_domain")
        self.assertEqual(native_layout["body_boundary_name"], "wall")
        self.assertEqual(
            native_layout["domain"],
            {
                "xmin": -3, "xmax": 6, "ymin": -3, "ymax": 3,
                "boundaries": {
                    "left": "inflow", "right": "outflow",
                    "bottom": "farfield", "top": "farfield",
                },
            },
        )
        self.assertEqual(
            native_layout["quality"],
            {"minimum_triangle_quality": 0.01, "minimum_angle_degrees": 1.0},
        )
        self.assertEqual(native_layout["domain"]["boundaries"], {
            "left": "inflow", "right": "outflow",
            "bottom": "farfield", "top": "farfield",
        })
        grid = (native / "script/grid.txt").read_text(encoding="utf-8")
        self.assertIn('string targetGridFileName = "results/circle_external_flow.msh";', grid)
        self.assertFalse((native / "results").exists())

    def test_square_2x2_reuses_square_stl_and_changes_only_mesh_density(self):
        square = ROOT / "grid/stl_to_gmsh/2d/square/single_quad/cpp"
        square_2x2 = ROOT / "grid/stl_to_gmsh/2d/square_2x2/cpp"
        self.assertEqual(
            (square / "input/model.stl").read_bytes(),
            (square_2x2 / "input/model.stl").read_bytes(),
        )
        self.assertTrue((square_2x2.parent / "README.md").is_file())
        layout = json.loads(
            (square_2x2 / "script/layout.json").read_text(encoding="utf-8")
        )
        self.assertEqual(layout["transfinite_cells"], {"x": 2, "y": 2})
        grid = (square_2x2 / "script/grid.txt").read_text(encoding="utf-8")
        self.assertIn('string sourceGridFileName = "input/model.stl";', grid)
        self.assertIn(
            'string targetGridFileName = "results/square_2x2.msh";',
            grid,
        )
        self.assertFalse((square_2x2 / "results").exists())

    def test_square_two_triangles_uses_same_stl_and_disables_recombination(self):
        native = ROOT / "grid/stl_to_gmsh/2d/square/two_triangles/cpp"
        reference = ROOT / "grid/stl_to_gmsh/2d/square/two_triangles/python"
        self.assertEqual((native / "input/model.stl").read_bytes(),
                         (reference / "input/model.stl").read_bytes())
        native_layout = json.loads((native / "script/layout.json").read_text(encoding="utf-8"))
        python_layout = json.loads((reference / "script/layout.json").read_text(encoding="utf-8"))
        self.assertEqual(native_layout, python_layout)
        self.assertEqual(native_layout["transfinite_cells"], {"x": 1, "y": 1})
        self.assertFalse(native_layout["recombine_surface"])
        with tempfile.TemporaryDirectory() as temporary:
            geo = Path(temporary) / "square_two_triangles.geo"
            convert(reference / "input/model.stl", reference / "script/layout.json",
                    geo, tolerance=1.0e-10)
            text = geo.read_text(encoding="utf-8")
            self.assertIn("Transfinite Surface {1};", text)
            self.assertNotIn("Recombine Surface {1};", text)
            self.assertIn('Physical Surface("fluid_domain")', text)
        expected = json.loads((reference / "expected.json").read_text(encoding="utf-8"))
        self.assertEqual(expected["mesh"], {
            "nodes": 4, "line_elements": 4, "triangle_elements": 2,
            "quadrilateral_elements": 0,
        })
        self.assertFalse((native / "results").exists())
        self.assertFalse((reference / "results").exists())

    def test_python_square_two_triangles_generates_exactly_two_triangle_cells(self):
        gmsh = _manual_gmsh()
        if gmsh is None:
            self.skipTest("Gmsh executable is not available")
        source = ROOT / "grid/stl_to_gmsh/2d/square/two_triangles/python"
        with tempfile.TemporaryDirectory() as temporary:
            staged = Path(temporary) / "square_two_triangles"
            shutil.copytree(source, staged)
            command = [sys.executable, str(ROOT / "grid/stl_to_gmsh/run_python_case.py"),
                       "--case-directory", str(staged), "--gmsh", gmsh]
            result = subprocess.run(command, cwd=staged, text=True, capture_output=True, check=False)
            self.assertEqual(result.returncode, 0, msg=result.stdout + "\\n" + result.stderr)
            mesh = inspect_msh4(staged / "results/square_two_triangles.msh")
            self.assertEqual(mesh["nodes"], 4)
            self.assertEqual(mesh["line_elements"], 4)
            self.assertEqual(mesh["triangle_elements"], 2)
            self.assertEqual(mesh["quadrilateral_elements"], 0)
            groups = {(item["dimension"], item["name"]) for item in mesh["physical_groups"]}
            self.assertEqual(groups, {
                (1, "inflow"), (1, "outflow"), (1, "symmetry"),
                (1, "farfield"), (2, "fluid_domain"),
            })
            report = json.loads(
                (staged / "results/square_two_triangles.stl-to-gmsh-report.json")
                .read_text(encoding="utf-8")
            )
            self.assertTrue(report["validated"])
        self.assertFalse((source / "results").exists())

    def test_python_geometry_compiles_square_and_two_cell_layouts(self):
        for relative, nx, ny in (
            ("2d/square/single_quad/python", 1, 1),
            ("2d/two_square_cells/python", 2, 1),
        ):
            case = ROOT / "grid/stl_to_gmsh" / relative
            with tempfile.TemporaryDirectory() as temp:
                geo = Path(temp) / "geometry.geo"
                report = convert(
                    case / "input/model.stl",
                    case / "script/layout.json",
                    geo,
                    tolerance=1.0e-10,
                )
                self.assertEqual(report["boundary_vertices"], 4)
                layout = json.loads(
                    (case / "script/layout.json").read_text(encoding="utf-8")
                )
                self.assertEqual(layout["transfinite_cells"], {"x": nx, "y": ny})
                text = geo.read_text(encoding="utf-8")
                for name in ("inflow", "outflow", "symmetry", "farfield"):
                    self.assertIn(name, text)
                self.assertIn('Physical Surface("fluid_domain")', text)
            self.assertFalse((case / "results").exists())

    def test_automatic_triangle_mode_uses_mesh_size_without_transfinite_constraints(self):
        source = ROOT / "grid/stl_to_gmsh/2d/square/two_triangles/python"
        with tempfile.TemporaryDirectory() as temporary:
            sidecar_path = Path(temporary) / "layout.json"
            sidecar = json.loads(
                (source / "script/layout.json").read_text(encoding="utf-8")
            )
            sidecar.pop("transfinite_cells")
            sidecar.pop("recombine_surface")
            sidecar["mesh_size"] = 0.25
            sidecar_path.write_text(json.dumps(sidecar), encoding="utf-8")
            from grid.stl_to_gmsh.converter import (
                build_geo, extract_boundary, read_stl, _read_sidecar,
            )
            loop = extract_boundary(read_stl(source / "input/model.stl"), 1.0e-10)
            config, assignments = _read_sidecar(sidecar_path, loop, 1.0e-10)
            generated = build_geo(loop, assignments, config)
            self.assertIn("Mesh.MeshSizeMin = 0.25;", generated)
            self.assertIn("Mesh.MeshSizeMax = 0.25;", generated)
            self.assertIn("Mesh.MeshSizeFromPoints = 1;", generated)
            self.assertNotIn("Transfinite Curve", generated)
            self.assertNotIn("Transfinite Surface", generated)
            self.assertNotIn("Recombine Surface", generated)
            self.assertIn('Physical Surface("fluid_domain")', generated)
            sidecar["recombine_surface"] = True
            sidecar_path.write_text(json.dumps(sidecar), encoding="utf-8")
            config, assignments = _read_sidecar(sidecar_path, loop, 1.0e-10)
            with self.assertRaisesRegex(ValueError, "does not support recombine_surface"):
                build_geo(loop, assignments, config)
        self.assertFalse((source / "results").exists())

    def test_python_automatic_mesh_size_controls_progressive_refinement(self):
        gmsh = _manual_gmsh()
        if gmsh is None:
            self.skipTest("Gmsh executable is not available")
        source = ROOT / "grid/stl_to_gmsh/2d/square/two_triangles/python"
        triangle_counts = []
        node_counts = []
        expected_groups = {
            (1, "inflow"), (1, "outflow"), (1, "symmetry"),
            (1, "farfield"), (2, "fluid_domain"),
        }
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            base = json.loads((source / "script/layout.json").read_text(encoding="utf-8"))
            base.pop("transfinite_cells")
            base.pop("recombine_surface")
            for index, mesh_size in enumerate((0.5, 0.25, 0.125)):
                layout = dict(base, mesh_size=mesh_size)
                sidecar = root / f"layout_{index}.json"
                sidecar.write_text(json.dumps(layout), encoding="utf-8")
                from grid.stl_to_gmsh.converter import build_geo, extract_boundary, read_stl, _read_sidecar
                loop = extract_boundary(read_stl(source / "input/model.stl"), 1.0e-10)
                config, assignments = _read_sidecar(sidecar, loop, 1.0e-10)
                geo = root / f"automatic_{index}.geo"
                geo.write_text(build_geo(loop, assignments, config), encoding="utf-8")
                msh = root / f"automatic_{index}.msh"
                result = subprocess.run(
                    [gmsh, str(geo), "-2", "-format", "msh4", "-o", str(msh)],
                    cwd=root, text=True, capture_output=True, check=False,
                )
                self.assertEqual(result.returncode, 0, msg=result.stdout + "\n" + result.stderr)
                mesh = inspect_msh4(msh)
                self.assertGreater(mesh["nodes"], 4)
                self.assertGreater(mesh["line_elements"], 4)
                self.assertGreater(mesh["triangle_elements"], 0)
                self.assertEqual(mesh["quadrilateral_elements"], 0)
                self.assertEqual(
                    {(group["dimension"], group["name"]) for group in mesh["physical_groups"]},
                    expected_groups,
                )
                triangle_counts.append(mesh["triangle_elements"])
                node_counts.append(mesh["nodes"])
        self.assertEqual(triangle_counts, sorted(triangle_counts))
        self.assertEqual(node_counts, sorted(node_counts))
        self.assertLess(triangle_counts[0], triangle_counts[-1])
        self.assertLess(node_counts[0], node_counts[-1])

    def test_expected_mesh_counts_distinguish_the_two_examples(self):
        square = json.loads(
            (ROOT / "grid/stl_to_gmsh/2d/square/single_quad/python/expected.json").read_text()
        )
        two = json.loads(
            (ROOT / "grid/stl_to_gmsh/2d/two_square_cells/python/expected.json").read_text()
        )
        self.assertEqual(
            square["mesh"],
            {"nodes": 4, "line_elements": 4, "triangle_elements": 0, "quadrilateral_elements": 1},
        )
        self.assertEqual(
            two["mesh"],
            {"nodes": 6, "line_elements": 6, "triangle_elements": 0, "quadrilateral_elements": 2},
        )

    def test_msh_inspector_measures_triangle_quality_and_domain_bounds(self):
        with tempfile.TemporaryDirectory() as temporary:
            msh = Path(temporary) / "one_triangle.msh"
            msh.write_text(
                "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n"
                '$PhysicalNames\n1\n2 1 "fluid_domain"\n$EndPhysicalNames\n'
                "$Nodes\n1 3 1 3\n2 1 0 3\n1 2 3\n"
                "0 0 0\n1 0 0\n0 1 0\n$EndNodes\n"
                "$Elements\n1 1 1 1\n2 1 2 1\n1 1 2 3\n$EndElements\n",
                encoding="utf-8",
            )
            summary = inspect_msh4(msh)
        quality = summary["triangle_quality"]
        self.assertEqual(summary["triangle_elements"], 1)
        self.assertEqual(quality["degenerate_count"], 0)
        self.assertAlmostEqual(quality["minimum_area"], 0.5)
        self.assertAlmostEqual(quality["minimum_quality"], 3.0 ** 0.5 / 2.0)
        self.assertAlmostEqual(quality["minimum_angle_degrees"], 45.0)
        self.assertEqual(summary["coordinate_bounds"]["minimum"], [0.0, 0.0, 0.0])
        self.assertEqual(summary["coordinate_bounds"]["maximum"], [1.0, 1.0, 0.0])

    def test_python_square_generates_mesh_and_report_without_domain_quality(self):
        gmsh = _manual_gmsh()
        if gmsh is None:
            self.skipTest("Gmsh executable is not available")
        source = ROOT / "grid/stl_to_gmsh/2d/square/single_quad/python"
        with tempfile.TemporaryDirectory() as temporary:
            staged = Path(temporary) / "square"
            shutil.copytree(source, staged)
            command = [
                sys.executable,
                str(ROOT / "grid/stl_to_gmsh/run_python_case.py"),
                "--case-directory",
                str(staged),
                "--gmsh",
                gmsh,
            ]
            result = subprocess.run(
                command, cwd=staged, text=True, capture_output=True, check=False
            )
            self.assertEqual(
                result.returncode, 0,
                msg=result.stdout + "\\n" + result.stderr,
            )
            mesh = inspect_msh4(staged / "results/square.msh")
            self.assertEqual(mesh["nodes"], 4)
            self.assertEqual(mesh["line_elements"], 4)
            self.assertEqual(mesh["quadrilateral_elements"], 1)
            report = json.loads(
                (staged / "results/square.stl-to-gmsh-report.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertTrue(report["validated"])
            self.assertEqual(report["mesh"], "results/square.msh")
            self.assertEqual(report["contract"]["quality_limits"], {})
            self.assertFalse((staged / "results/model.msh").exists())
        self.assertFalse((source / "results").exists())

    def test_python_circle_external_flow_generates_valid_gmsh_mesh(self):
        gmsh = _manual_gmsh()
        if gmsh is None:
            self.skipTest("Gmsh executable is not available")
        source = ROOT / "grid/stl_to_gmsh/2d/circle_external_flow/python"
        with tempfile.TemporaryDirectory() as temporary:
            staged = Path(temporary) / "circle_external_flow"
            shutil.copytree(source, staged)
            command = [
                sys.executable,
                str(ROOT / "grid/stl_to_gmsh/run_python_case.py"),
                "--case-directory",
                str(staged),
                "--gmsh",
                gmsh,
            ]
            result = subprocess.run(
                command, cwd=staged, text=True, capture_output=True, check=False
            )
            self.assertEqual(
                result.returncode, 0,
                msg=result.stdout + "\\n" + result.stderr,
            )
            mesh = inspect_msh4(staged / "results/circle_external_flow.msh")
            self.assertGreater(mesh["nodes"], 20)
            self.assertGreater(mesh["line_elements"], 0)
            self.assertGreater(mesh["triangle_elements"], 0)
            self.assertEqual(mesh["quadrilateral_elements"], 0)
            quality = mesh["triangle_quality"]
            self.assertEqual(quality["degenerate_count"], 0)
            self.assertGreaterEqual(quality["minimum_quality"], 0.01)
            self.assertGreaterEqual(quality["minimum_angle_degrees"], 1.0)
            self.assertEqual(mesh["coordinate_bounds"]["minimum"][:2], [-3.0, -3.0])
            self.assertEqual(mesh["coordinate_bounds"]["maximum"][:2], [6.0, 3.0])
            groups = {
                (group["dimension"], group["name"])
                for group in mesh["physical_groups"]
            }
            self.assertEqual(groups, {
                (1, "wall"), (1, "inflow"), (1, "outflow"),
                (1, "farfield"), (2, "fluid_domain"),
            })
            report = json.loads(
                (staged / "results/circle_external_flow.stl-to-gmsh-report.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertTrue(report["validated"])
            self.assertEqual(report["mesh"], "results/circle_external_flow.msh")
            self.assertFalse((staged / "results/model.msh").exists())
        self.assertFalse((source / "results").exists())

    def test_external_and_internal_flow_cases_build_distinct_fluid_domains(self):
        cases_root = ROOT / "grid/stl_to_gmsh"
        external = cases_root / "2d/circle_external_flow/python"
        internal = cases_root / "2d/channel_internal_flow/python"
        with tempfile.TemporaryDirectory() as temp:
            external_geo = Path(temp) / "external.geo"
            external_report = convert(
                external / "input/model.stl",
                external / "script/layout.json",
                external_geo,
                tolerance=1.0e-10,
            )
            external_text = external_geo.read_text(encoding="utf-8")
            self.assertEqual(external_report["body_boundary_edges"], 32)
            self.assertIn("Plane Surface(1) = {1, 2};", external_text)
            self.assertIn("lc = 0.22;", external_text)
            self.assertIn("lc_outer =", external_text)
            self.assertIn("Point(33) = {-3, -3, 0, lc_outer};", external_text)
            self.assertIn('Physical Curve("wall")', external_text)
            self.assertEqual(external_text.count('Physical Curve("farfield")'), 1)
            for name in ("inflow", "outflow", "farfield", "wall"):
                self.assertIn(f'Physical Curve("{name}")', external_text)
            self.assertIn('Physical Surface("fluid_domain")', external_text)

            internal_geo = Path(temp) / "internal.geo"
            internal_report = convert(
                internal / "input/model.stl",
                internal / "script/layout.json",
                internal_geo,
                tolerance=1.0e-10,
            )
            internal_text = internal_geo.read_text(encoding="utf-8")
            self.assertEqual(internal_report["boundary_edges"], 4)
            self.assertIn("Plane Surface(1) = {1};", internal_text)
            self.assertNotIn("Plane Surface(1) = {1, 2};", internal_text)
            for name in ("Inlet", "Outlet", "UpperWall", "LowerWall"):
                self.assertIn(f'Physical Curve("{name}")', internal_text)
            self.assertIn("Transfinite Curve", internal_text)
            self.assertIn("Transfinite Surface {1};", internal_text)
            self.assertIn("Recombine Surface {1};", internal_text)

        for case in (external, internal):
            self.assertFalse((case / "results").exists())


    def test_python_cube_single_hex_generates_one_hexahedron(self):
        gmsh = _manual_gmsh()
        if gmsh is None:
            self.skipTest("Gmsh executable is not available")
        source = ROOT / "grid/stl_to_gmsh/3d/cube/single_hex/python"
        with tempfile.TemporaryDirectory() as temporary:
            staged = Path(temporary) / "cube_single_hex"
            shutil.copytree(source, staged)
            command = [
                sys.executable, str(ROOT / "grid/stl_to_gmsh/run_python_case.py"),
                "--case-directory", str(staged), "--gmsh", gmsh,
            ]
            result = subprocess.run(command, cwd=staged, text=True, capture_output=True, check=False)
            self.assertEqual(result.returncode, 0, msg=result.stdout + "\\n" + result.stderr)
            mesh = inspect_msh4(staged / "results/cube_single_hex.msh")
            self.assertEqual(mesh["nodes"], 8)
            self.assertEqual(mesh["quadrilateral_elements"], 6)
            self.assertEqual(mesh["hexahedron_elements"], 1)
            self.assertEqual(mesh["triangle_elements"], 0)
            groups = {(item["dimension"], item["name"]) for item in mesh["physical_groups"]}
            self.assertEqual(groups, {
                (2, "xmin_inflow"),
                (2, "xmax_outflow"),
                (2, "ymin_wall"),
                (2, "ymax_farfield"),
                (2, "zmin_symmetry"),
                (2, "zmax_symmetry"),
                (3, "fluid_domain"),
            })
            self.assertTrue(json.loads(
                (staged / "results/cube_single_hex.stl-to-gmsh-report.json").read_text()
            )["validated"])
        self.assertFalse((source / "results").exists())

    def test_python_prolate_spheroid_generates_valid_tetrahedral_volume_mesh(self):
        gmsh = _manual_gmsh()
        if gmsh is None:
            self.skipTest("Gmsh executable is not available")
        source = ROOT / "grid/stl_to_gmsh/3d/prolate_spheroid/automatic_tetrahedra/python"
        with tempfile.TemporaryDirectory() as temporary:
            staged = Path(temporary) / "prolate_spheroid"
            shutil.copytree(source, staged)
            command = [
                sys.executable,
                str(ROOT / "grid/stl_to_gmsh/run_python_case.py"),
                "--case-directory",
                str(staged),
                "--gmsh",
                gmsh,
            ]
            result = subprocess.run(
                command, cwd=staged, text=True, capture_output=True, check=False
            )
            self.assertEqual(
                result.returncode, 0,
                msg=result.stdout + "\n" + result.stderr,
            )
            mesh = inspect_msh4(
                staged / "results/prolate_spheroid_automatic_tetrahedra.msh"
            )
            self.assertGreater(mesh["nodes"], 0)
            self.assertGreater(mesh["triangle_elements"], 0)
            self.assertGreater(mesh["tetrahedron_elements"], 0)
            self.assertEqual(mesh["hexahedron_elements"], 0)
            bounds = mesh["coordinate_bounds"]
            expected_minimum = (-3.0, -0.5, -0.5)
            expected_maximum = (3.0, 0.5, 0.5)
            # Gmsh reconstructs parametric patches from STL triangles; the
            # generated envelope is an approximation of the STL vertex bounds.
            for actual, expected in zip(bounds["minimum"], expected_minimum):
                self.assertAlmostEqual(
                    actual, expected, delta=max(abs(expected) * 0.10, 0.02)
                )
            for actual, expected in zip(bounds["maximum"], expected_maximum):
                self.assertAlmostEqual(
                    actual, expected, delta=max(abs(expected) * 0.10, 0.02)
                )
            extents = [
                upper - lower
                for lower, upper in zip(bounds["minimum"], bounds["maximum"])
            ]
            aspect_ratio = extents[0] / max(extents[1], extents[2])
            self.assertGreater(aspect_ratio, 5.0)
            self.assertLess(aspect_ratio, 7.0)
            groups = {
                (item["dimension"], item["name"])
                for item in mesh["physical_groups"]
            }
            self.assertEqual(groups, {(2, "wall"), (3, "fluid_domain")})
            report = json.loads(
                (
                    staged / "results/prolate_spheroid_automatic_tetrahedra"
                    ".stl-to-gmsh-report.json"
                ).read_text(encoding="utf-8")
            )
            self.assertTrue(report["validated"])
        self.assertFalse((source / "results").exists())

    def test_python_internal_flow_cylinder_generates_region_aware_tetrahedral_mesh(self):
        gmsh = _manual_gmsh()
        if gmsh is None:
            self.skipTest("Gmsh executable is not available")
        source = (
            ROOT / "grid/stl_to_gmsh/3d/cylinder_internal_flow"
            / "automatic_tetrahedra/python"
        )
        with tempfile.TemporaryDirectory() as temporary:
            staged = Path(temporary) / "cylinder_internal_flow"
            shutil.copytree(source, staged)
            command = [
                sys.executable,
                str(ROOT / "grid/stl_to_gmsh/run_python_case.py"),
                "--case-directory",
                str(staged),
                "--gmsh",
                gmsh,
            ]
            result = subprocess.run(
                command, cwd=staged, text=True, capture_output=True, check=False
            )
            self.assertEqual(
                result.returncode, 0,
                msg=result.stdout + "\\n" + result.stderr,
            )
            mesh = inspect_msh4(
                staged / "results/cylinder_internal_flow.msh"
            )
            self.assertGreater(mesh["nodes"], 0)
            self.assertGreater(mesh["triangle_elements"], 0)
            self.assertGreater(mesh["tetrahedron_elements"], 0)
            self.assertEqual(mesh["quadrilateral_elements"], 0)
            self.assertEqual(mesh["hexahedron_elements"], 0)
            bounds = mesh["coordinate_bounds"]
            for actual, expected in zip(
                bounds["minimum"], (-2.0, -1.0, -1.0)
            ):
                self.assertAlmostEqual(actual, expected, delta=1.0e-8)
            for actual, expected in zip(
                bounds["maximum"], (2.0, 1.0, 1.0)
            ):
                self.assertAlmostEqual(actual, expected, delta=1.0e-8)
            groups = {
                (item["dimension"], item["name"])
                for item in mesh["physical_groups"]
            }
            self.assertEqual(
                groups,
                {
                    (2, "wall"),
                    (2, "inflow"),
                    (2, "outflow"),
                    (3, "fluid_domain"),
                },
            )
            report = json.loads(
                (
                    staged / "results/cylinder_internal_flow"
                    ".stl-to-gmsh-report.json"
                ).read_text(encoding="utf-8")
            )
            self.assertTrue(report["validated"])
            self.assertEqual(
                report["geometry_report"]["geometry"],
                "cylinder_internal_flow",
            )
        self.assertFalse((source / "results").exists())

    def test_pipe_bend_geometry_contract_emits_verified_region_assets(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "pipe_bend_90_degree.stl"
            generated = generate_pipe_bend(output)
            self.assertTrue(generated["surface"]["watertight"])
            self.assertGreater(generated["surface"]["signed_volume"], 0.0)
            self.assertLess(generated["surface"]["relative_volume_error"], 0.03)

            geometry_manifest = json.loads(
                generated["manifest"].read_text(encoding="utf-8")
            )
            surface_sha256 = hashlib.sha256(output.read_bytes()).hexdigest()
            self.assertEqual(
                geometry_manifest["schema"], "oneflow.stl-geometry-manifest"
            )
            self.assertEqual(geometry_manifest["stl_sha256"], surface_sha256)
            self.assertEqual(geometry_manifest["units"], "m")
            self.assertIn("parameters", geometry_manifest)

            region_manifest = json.loads(
                generated["region_manifest_path"].read_text(encoding="utf-8")
            )
            self.assertEqual(
                region_manifest["schema"], "oneflow.stl-surface-regions"
            )
            self.assertEqual(
                region_manifest["closed_surface_sha256"], surface_sha256
            )
            self.assertEqual(
                set(region_manifest["regions"]), {"wall", "inflow", "outflow"}
            )
            total_region_triangles = 0
            for name, record in region_manifest["regions"].items():
                region_file = output.parent / record["file"]
                self.assertTrue(region_file.is_file(), name)
                self.assertEqual(
                    hashlib.sha256(region_file.read_bytes()).hexdigest(),
                    record["sha256"],
                )
                self.assertGreater(record["triangle_count"], 0)
                total_region_triangles += record["triangle_count"]
            self.assertEqual(
                total_region_triangles, generated["surface"]["triangle_count"]
            )


if __name__ == "__main__":
    unittest.main()
