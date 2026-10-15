#!/usr/bin/env python3
"""Unit tests for the OneFLOW grid case catalog."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from mesh_cases import load_catalog, validate_case, validate_profile


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


class MeshCaseCatalogTest(unittest.TestCase):
    def test_repository_catalog_entries_resolve(self) -> None:
        catalog = load_catalog(REPOSITORY_ROOT / "grid" / "case_catalog.json")
        self.assertGreaterEqual(len(catalog["cases"]), 4)
        profile_ids = {
            entry["id"] for entry in catalog["cases"] if entry.get("profile")
        }
        self.assertTrue({"minimal_rectangle", "minimal_cube", "two_connected_cubes"} <= profile_ids)
        errors = [
            error
            for entry in catalog["cases"]
            for error in validate_case(REPOSITORY_ROOT, entry)
        ]
        self.assertEqual(errors, [])

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
