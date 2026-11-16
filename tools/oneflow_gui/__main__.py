"""Launch the OneFLOW desktop application as a Python module."""

from __future__ import annotations

import os
import sys
from pathlib import Path


def _configure_bytecode_cache() -> None:
    """Keep Python bytecode caches in user data, never in the source tree."""
    cache_root = Path(
        os.environ.get("LOCALAPPDATA", str(Path.home() / ".oneflow"))
    ) / "OneFLOW" / "pycache"
    cache_root.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("PYTHONPYCACHEPREFIX", str(cache_root))
    sys.pycache_prefix = os.environ["PYTHONPYCACHEPREFIX"]


if __name__ == "__main__":
    _configure_bytecode_cache()
    from .main_window import run_desktop_app

    raise SystemExit(run_desktop_app())
