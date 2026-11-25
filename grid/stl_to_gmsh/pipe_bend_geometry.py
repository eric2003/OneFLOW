"""Importable facade for the registered 90-degree pipe-bend generator."""
from __future__ import annotations

import importlib.util
from pathlib import Path

_GENERATOR = (
    Path(__file__).resolve().parent
    / "3d/pipe_bend/90_degree/geometry/python/generate_pipe_bend.py"
)
_SPEC = importlib.util.spec_from_file_location("oneflow_pipe_bend_geometry_impl", _GENERATOR)
if _SPEC is None or _SPEC.loader is None:
    raise ImportError(f"Cannot load pipe-bend geometry generator: {_GENERATOR}")
_MODULE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MODULE)

build_surface = _MODULE.build_surface
inspect_surface = _MODULE.inspect_surface
generate = _MODULE.generate
main = _MODULE.main

__all__ = ["build_surface", "inspect_surface", "generate", "main"]
