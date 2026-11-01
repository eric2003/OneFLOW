"""Table-level geometry contract workbench (stdlib core; optional plots)."""

from __future__ import annotations

from typing import Any

__all__ = ["diff_tables", "run_pipeline"]


def __getattr__(name: str) -> Any:
    """Load the optional scenario runner only when its public API is requested.

    Importing lightweight scenario modules must not pull in the table runner
    or its optional analytic dependencies.
    """
    if name in __all__:
        from . import pipeline_stub

        return getattr(pipeline_stub, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
