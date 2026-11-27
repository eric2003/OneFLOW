"""Import adapter for the OECD/NEA-Vattenfall T-junction STL generator.

The benchmark is grouped under a human-readable 3d directory; this adapter
keeps its reusable generator importable from Python modules.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

_GENERATOR_PATH = (
    Path(__file__).resolve().parent
    / "3d" / "pipe_tee" / "oecd_nea_vattenfall"
    / "geometry" / "python" / "generate_pipe_tee.py"
)
_SPEC = importlib.util.spec_from_file_location(
    "oneflow_pipe_tee_geometry_generator", _GENERATOR_PATH
)
if _SPEC is None or _SPEC.loader is None:
    raise ImportError(f"Cannot load T-junction geometry generator: {_GENERATOR_PATH}")
_MODULE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MODULE)

generate = _MODULE.generate
export_surface = _MODULE.export_surface
