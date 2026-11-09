"""Mesh-to-CFD validation Case catalog and graph tests."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.mesh_validation import (
    build_mesh_validation_cases, cfd_template_choices, mesh_case_choices,
)


class MeshValidationCatalogTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        (self.root / "grid/stl_to_gmsh/square/script").mkdir(parents=True)
        (self.root / "test/cfd/plateuns2dslau2/script").mkdir(parents=True)
        (self.root / "grid/case_catalog.json").write_text(json.dumps({
            "schema": "oneflow.case-catalog",
            "cases": [
                {"id": "stl_to_gmsh_square", "title": "One Square Cell",
                 "directory": "grid/stl_to_gmsh/square", "group": "STL to Gmsh"},
                {"id": "python_square", "directory": "grid/python",
                 "group": "STL to Gmsh",
                 "execution": {"kind": "python", "script": "run.py"}},
            ],
        }), encoding="utf-8")
        (self.root / "grid/stl_to_gmsh/square/script/grid.txt").write_text(
            'string targetGridFileName = "results/square.msh";\n', encoding="utf-8")
        (self.root / "test/cfd/plateuns2dslau2/script/cfd.txt").write_text(
            'string resFile = "results/res.dat";\n'
            'string flowFile = "restart/flow.dat";\n'
            'string visualFile = "visual/visual.dat";\n', encoding="utf-8")

    def tearDown(self):
        self.temporary.cleanup()

    def test_mesh_selector_contains_native_cases_only(self):
        self.assertEqual(mesh_case_choices(self.root), (("stl_to_gmsh_square", "One Square Cell"),))

    def test_first_cfd_template_is_the_existing_laminar_flat_plate_case(self):
        self.assertEqual(cfd_template_choices(), (("laminar_flat_plate", "Laminar Flat Plate"),))

    def test_selected_case_builds_conversion_chain_and_real_cfd_consumer(self):
        cases = build_mesh_validation_cases(self.root, "stl_to_gmsh_square", "laminar_flat_plate")
        self.assertEqual(len(cases), 4)
        self.assertEqual(cases[0]["outputs"][0]["path"], "results/square.msh")
        self.assertEqual(cases[1]["directory"], "grid/gmsh_to_cgns")
        self.assertEqual(cases[2]["directory"], "grid/stl_to_gmsh/pipeline/cgns_to_oneflow")
        self.assertEqual(cases[3]["directory"], "test/cfd/plateuns2dslau2")
        self.assertEqual(cases[3]["inputs"][0]["artifact"], "stl_to_gmsh_square_native")

    def test_unknown_mesh_case_and_template_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "Unknown mesh-generation Case"):
            build_mesh_validation_cases(self.root, "missing", "laminar_flat_plate")
        with self.assertRaisesRegex(ValueError, "Unknown CFD template"):
            build_mesh_validation_cases(self.root, "stl_to_gmsh_square", "missing")


if __name__ == "__main__":
    unittest.main()
