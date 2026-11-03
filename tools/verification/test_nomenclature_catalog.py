"""Tests for the source-backed CFD nomenclature catalog."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from tools.verification.nomenclature_catalog import (
    NomenclatureCatalog,
    validate_nomenclature_catalog,
)


class NomenclatureCatalogTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalog = NomenclatureCatalog.load()

    def test_all_aliases_have_resolvable_sources_and_unique_concepts(self):
        self.assertEqual(validate_nomenclature_catalog(self.catalog), [])
        self.assertEqual(len({entry.id for entry in self.catalog.entries}), len(self.catalog.entries))

    def test_cgns_static_pressure_is_a_standardized_flow_quantity(self):
        matches = self.catalog.find_alias("Pressure", system="CGNS")
        self.assertEqual(len(matches), 1)
        entry, alias = matches[0]
        self.assertEqual(entry.id, "fluid.static_pressure")
        self.assertEqual(alias.relation, "exact")
        self.assertIn("DataClass", alias.caveat)

    def test_openfoam_p_is_intentionally_ambiguous_without_solver_context(self):
        matches = self.catalog.find_alias("p", system="OpenFOAM")
        self.assertEqual(
            {entry.id for entry, _ in matches},
            {"fluid.static_pressure", "fluid.kinematic_pressure"},
        )
        self.assertTrue(any(
            "kinematic" in alias.applicability for _, alias in matches
        ))

    def test_su2_freestream_values_are_parameters_not_solution_fields(self):
        for name, expected in (
            ("FREESTREAM_PRESSURE", "fluid.static_pressure"),
            ("FREESTREAM_DENSITY", "fluid.density"),
            ("FREESTREAM_VELOCITY", "fluid.velocity"),
        ):
            matches = self.catalog.find_alias(name, system="SU2")
            self.assertEqual(len(matches), 1)
            entry, alias = matches[0]
            self.assertEqual(entry.id, expected)
            self.assertEqual(alias.role, "input_parameter")
            self.assertEqual(alias.relation, "related")

    def test_fluent_gauge_and_absolute_pressure_are_not_merged(self):
        gauge = self.catalog.find_alias("Static Pressure", system="ANSYS Fluent")[0][0]
        absolute = self.catalog.find_alias("Absolute Pressure", system="ANSYS Fluent")[0][0]
        self.assertEqual(gauge.id, "fluid.static_pressure")
        self.assertEqual(absolute.id, "fluid.absolute_pressure")
        self.assertIn("Operating Pressure", self.catalog.find_alias(
            "Static Pressure", system="ANSYS Fluent"
        )[0][1].caveat)

    def test_total_pressure_names_map_to_stagnation_pressure_concept(self):
        for system, name in (
            ("CGNS", "PressureStagnation"),
            ("ANSYS Fluent", "Total Pressure"),
            ("OpenFOAM", "total(p)"),
        ):
            self.assertEqual(
                self.catalog.find_alias(name, system=system)[0][0].id,
                "fluid.total_pressure",
            )

    def test_velocity_and_mesh_coordinate_names_are_searchable(self):
        self.assertEqual(self.catalog.find_alias("U", system="OpenFOAM")[0][0].id, "fluid.velocity")
        self.assertEqual(self.catalog.find_alias("CoordinateX", system="CGNS")[0][0].id, "mesh.coordinate")
        self.assertEqual(
            [entry.id for entry in self.catalog.search("stagnation")],
            ["fluid.total_pressure"],
        )

    def test_missing_source_reference_is_rejected(self):
        payload = {
            "schema_version": 1,
            "sources": [{"id": "s", "name": "source", "url": "https://example.invalid"}],
            "concepts": [{
                "id": "x", "label": "x", "symbol": "x", "definition": "x",
                "dimensions": {"L": 1}, "si_unit": "m",
                "aliases": [{
                    "system": "test", "name": "x", "role": "field", "relation": "exact",
                    "applicability": "always", "caveat": "", "source_ids": ["missing"]
                }]
            }]
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "catalog.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "unknown source IDs"):
                NomenclatureCatalog.load(path)


if __name__ == "__main__":
    unittest.main()
