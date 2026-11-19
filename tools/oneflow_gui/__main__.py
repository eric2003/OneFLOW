"""Launch the OneFLOW desktop application as a Python module."""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Mapping


def _bytecode_cache_root(
    *,
    platform: str | None = None,
    environ: Mapping[str, str] | None = None,
    home: Path | None = None,
) -> Path:
    """Return the platform-standard cache location for OneFLOW bytecode."""
    platform = sys.platform if platform is None else platform
    environ = os.environ if environ is None else environ
    home = Path.home() if home is None else Path(home)

    if platform == "win32":
        base = Path(environ.get("LOCALAPPDATA", home / "AppData" / "Local"))
        return base / "OneFLOW" / "pycache"

    if platform == "darwin":
        return home / "Library" / "Caches" / "OneFLOW" / "pycache"

    xdg_cache_home = environ.get("XDG_CACHE_HOME")
    if xdg_cache_home:
        configured_cache = Path(xdg_cache_home).expanduser()
        if configured_cache.is_absolute():
            base = configured_cache
        else:
            base = home / ".cache"
    else:
        base = home / ".cache"
    return base / "oneflow" / "pycache"


def _configure_bytecode_cache() -> None:
    """Keep Python bytecode caches outside the source tree."""
    configured = os.environ.get("PYTHONPYCACHEPREFIX")
    cache_root = (
        Path(configured).expanduser()
        if configured
        else _bytecode_cache_root()
    ).resolve()
    cache_root.mkdir(parents=True, exist_ok=True)
    os.environ["PYTHONPYCACHEPREFIX"] = str(cache_root)
    sys.pycache_prefix = str(cache_root)


if __name__ == "__main__":
    _configure_bytecode_cache()
    from .main_window import run_desktop_app

    raise SystemExit(run_desktop_app())
