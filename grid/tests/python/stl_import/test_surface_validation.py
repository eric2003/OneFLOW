"""Tests for generic STL surface diagnostics."""
import math
import unittest
from grid.stl_import.surface import StlGeometryError
from grid.stl_import.validation import analyze_surface


def tetrahedron(scale=1.0):
    a, b = (0.,0.,0.), (scale,0.,0.)
    c, d = (0.,scale,0.), (0.,0.,scale)
    return [(a,c,b), (a,b,d), (a,d,c), (b,c,d)]


class StlSurfaceValidationTest(unittest.TestCase):
    def test_closed_tetrahedron_is_volume_mesh_ready(self):
        r = analyze_surface(tetrahedron())
        self.assertEqual(r["status"], "valid")
        self.assertTrue(r["closed"] and r["edge_manifold"] and r["consistently_oriented"])
        self.assertTrue(r["volume_mesh_ready"])
        self.assertEqual(r["components"][0]["genus"], 0)
        self.assertGreater(r["components"][0]["signed_volume"], 0)

    def test_open_surface_is_warning_not_import_failure(self):
        r = analyze_surface(tetrahedron()[:3])
        self.assertEqual(r["status"], "warning")
        self.assertFalse(r["closed"] or r["volume_mesh_ready"])
        self.assertEqual(r["boundary_edge_count"], 3)
        self.assertTrue(r["boundary"]["boundary_is_loop_like"])

    def test_disconnected_closed_shells_are_reported(self):
        first = tetrahedron()
        second = [tuple((p[0]+3.,p[1],p[2]) for p in f) for f in tetrahedron()]
        r = analyze_surface(first + second)
        self.assertEqual(r["status"], "warning")
        self.assertTrue(r["closed"])
        self.assertEqual(r["connected_component_count"], 2)
        self.assertFalse(r["volume_mesh_ready"])

    def test_duplicate_and_degenerate_facets_are_detected(self):
        f = ((0.,0.,0.),(1.,0.,0.),(0.,1.,0.))
        d = ((0.,0.,0.),(0.5,0.,0.),(1.,0.,0.))
        r = analyze_surface([f,f,d])
        self.assertEqual(r["status"], "invalid")
        self.assertEqual(r["duplicate_triangle_count"], 1)
        self.assertEqual(r["degenerate_triangle_count"], 1)

    def test_non_manifold_edge_is_detected(self):
        a,b = (0.,0.,0.),(1.,0.,0.)
        faces = [(a,b,(0.,1.,0.)),(b,a,(0.,-1.,0.)),(a,b,(0.,0.,1.))]
        r = analyze_surface(faces)
        self.assertEqual(r["status"], "invalid")
        self.assertGreater(r["non_manifold_edge_count"], 0)

    def test_inconsistent_winding_is_detected(self):
        faces = tetrahedron()
        faces[1] = tuple(reversed(faces[1]))
        r = analyze_surface(faces)
        self.assertEqual(r["status"], "invalid")
        self.assertGreater(r["orientation_conflict_edge_count"], 0)
        self.assertFalse(r["volume_mesh_ready"])

    def test_default_tolerance_scales_with_model(self):
        small, large = analyze_surface(tetrahedron(1e-3)), analyze_surface(tetrahedron(1e3))
        self.assertAlmostEqual(large["weld_tolerance"]/small["weld_tolerance"], 1e6, places=5)

    def test_explicit_tolerance_must_be_positive_finite(self):
        for tol in (0., -1., math.inf, math.nan):
            with self.assertRaises(StlGeometryError):
                analyze_surface(tetrahedron(), tol)

    def test_empty_surface_is_rejected(self):
        with self.assertRaisesRegex(StlGeometryError, "no triangles"):
            analyze_surface([])

    def test_inspect_cli_writes_real_newlines_and_keeps_open_surface_as_warning(self):
        import json
        import subprocess
        import sys
        import tempfile
        from pathlib import Path

        stl = (
            "solid open_surface\n"
            "facet normal 0 0 1\nouter loop\n"
            "vertex 0 0 0\nvertex 1 0 0\nvertex 0 1 0\n"
            "endloop\nendfacet\nendsolid open_surface\n"
        )
        with tempfile.TemporaryDirectory() as temp:
            source, output = Path(temp) / "open.stl", Path(temp) / "report.json"
            source.write_text(stl, encoding="ascii", newline="\n")
            result = subprocess.run(
                [sys.executable, str(Path(__file__).resolve().parents[3] / "stl_import" / "inspect.py"),
                 "--input", str(source), "--output", str(output)],
                capture_output=True, text=True, check=False,
            )
            self.assertEqual(result.returncode, 0, msg=result.stdout + result.stderr)
            payload = output.read_bytes()
            self.assertTrue(payload.endswith(b"\n"))
            self.assertNotIn(b"\\n", payload)
            report = json.loads(payload.decode("utf-8"))
            self.assertEqual(report["surface_quality"]["status"], "warning")


    def test_cylinder_internal_flow_stl_is_closed_without_patch_manifest(self):
        """Stage A inspects only model.stl; patch metadata is deliberately unused."""
        import json
        import subprocess
        import sys
        import tempfile
        from pathlib import Path

        repository_root = Path(__file__).resolve().parents[4]
        source = (
            repository_root
            / "grid" / "stl_to_gmsh" / "3d" / "cylinder_internal_flow"
            / "automatic_tetrahedra" / "python" / "input" / "model.stl"
        )
        inspector = repository_root / "grid" / "stl_import" / "inspect.py"
        self.assertTrue(source.is_file(), msg=f"Missing cylinder STL: {source}")

        # Pass exactly one geometry input to Stage A. Do not load model.geometry.json.
        with tempfile.TemporaryDirectory() as temp:
            report_path = Path(temp) / "stage_a_report.json"
            result = subprocess.run(
                [
                    sys.executable, str(inspector),
                    "--input", str(source),
                    "--output", str(report_path),
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(result.returncode, 0, msg=result.stdout + result.stderr)
            report = json.loads(report_path.read_text(encoding="utf-8"))

        quality = report["surface_quality"]
        self.assertEqual(report["format"], "binary")
        self.assertEqual(report["triangle_count"], 1632)
        self.assertEqual(quality["status"], "valid")
        self.assertTrue(quality["closed"])
        self.assertTrue(quality["edge_manifold"])
        self.assertTrue(quality["consistently_oriented"])
        self.assertTrue(quality["volume_mesh_ready"])
        self.assertEqual(quality["boundary_edge_count"], 0)
        self.assertEqual(quality["non_manifold_edge_count"], 0)
        self.assertEqual(quality["orientation_conflict_edge_count"], 0)
        self.assertEqual(quality["connected_component_count"], 1)
        self.assertEqual(len(quality["components"]), 1)
        self.assertEqual(quality["components"][0]["euler_characteristic"], 2)
        self.assertEqual(quality["components"][0]["genus"], 0)
        self.assertTrue(quality["components"][0]["outward_winding"])
        self.assertAlmostEqual(
            quality["components"][0]["signed_volume"], 12.5305, delta=0.002
        )


if __name__ == "__main__":
    unittest.main()
