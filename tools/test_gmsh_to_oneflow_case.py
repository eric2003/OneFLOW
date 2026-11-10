#!/usr/bin/env python3
"""Preflight regression for the checked-in Gmsh-to-OneFLOW example case."""

from __future__ import annotations

import subprocess
import sys
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


class GmshToOneFlowCaseTest(unittest.TestCase):
    def test_case_is_catalogued_and_passes_preflight(self) -> None:
        case_dir = REPOSITORY_ROOT / "grid" / "conversion" / "gmsh_to_oneflow"
        validator = (
            REPOSITORY_ROOT
            / "grid"
            / "mesh_workflow"
            / "tools"
            / "validate_case.py"
        )

        result = subprocess.run(
            [sys.executable, str(validator), str(case_dir)],
            cwd=REPOSITORY_ROOT,
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(
            result.returncode,
            0,
            msg=result.stdout + result.stderr,
        )
        self.assertIn("Source: gmsh / grid/single_quad.msh", result.stdout)
        self.assertIn(
            "Target: oneflow / results/single_quad.oneflow.grid",
            result.stdout,
        )

        import json

        catalog_path = REPOSITORY_ROOT / "grid" / "case_catalog.json"
        catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
        entry = next(
            (item for item in catalog["cases"] if item["id"] == "gmsh_to_oneflow"),
            None,
        )
        self.assertIsNotNone(entry)
        self.assertEqual(entry["directory"], "grid/conversion/gmsh_to_oneflow")
        self.assertEqual(entry["runtimeConfig"], "script/control.txt")


if __name__ == "__main__":
    unittest.main()
