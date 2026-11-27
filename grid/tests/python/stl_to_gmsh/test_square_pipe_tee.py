from __future__ import annotations
import importlib.util, json, tempfile, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4]
SCRIPT=ROOT/"grid/stl_to_gmsh/3d/pipe_tee/square_duct/geometry/python/generate_square_tee.py"
SPEC=importlib.util.spec_from_file_location("generate_square_tee",SCRIPT)
MODULE=importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)
class SquarePipeTeeGeometryTest(unittest.TestCase):
 def test_manual_shell_is_closed_and_preserves_three_ports(self):
  with tempfile.TemporaryDirectory() as temp:
   result=MODULE.generate(Path(temp)/"square_tee.stl","main_inlets")
   self.assertEqual(result["surface"]["invalid_edge_count"],0)
   self.assertAlmostEqual(result["surface"]["signed_volume_m3"],6.0,places=5)
   self.assertEqual(set(result["surface"]["port_subregion_triangle_counts"]),{"main_left_inflow","main_right_inflow","branch_outflow"})
   manifest=json.loads(result["region_manifest"].read_text(encoding="utf-8"))
   self.assertEqual(set(manifest["regions"]),{"wall","inflow","outflow"})
   self.assertEqual(set(manifest["subregions"]),{"main_left_inflow","main_right_inflow","branch_outflow"})
   for record in (*manifest["regions"].values(),*manifest["subregions"].values()):
    patch=Path(temp)/record["file"]; self.assertTrue(patch.is_file()); self.assertEqual(MODULE.sha(patch),record["sha256"])
 def test_branch_inlet_mode_keeps_three_port_names(self):
  with tempfile.TemporaryDirectory() as temp:
   result=MODULE.generate(Path(temp)/"square_tee.stl","branch_inlet")
   self.assertEqual(set(result["surface"]["port_subregion_triangle_counts"]),{"branch_inflow","main_left_outflow","main_right_outflow"})
if __name__=="__main__": unittest.main()
