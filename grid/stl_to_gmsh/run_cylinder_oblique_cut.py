#!/usr/bin/env python3
"""Compatibility entry point for the oblique-cut cylinder Python reference case."""
from pathlib import Path
import runpy

if __name__ == "__main__":
    target = (
        Path(__file__).resolve().parent / "3d" / "cylinder_internal_flow"
        / "oblique_cut" / "python" / "run_cylinder_oblique_cut.py"
    )
    runpy.run_path(str(target), run_name="__main__")
