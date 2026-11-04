"""File-independent resolution of logical parameter overrides for Case construction.

This module deliberately knows nothing about cfd.txt (or any other legacy
filename). It resolves keys through ParameterModel, then routes each update to
the document that owns the effective declaration.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Mapping

from .parameter_document import ParameterDocument
from .parameter_model import ParameterModel


@dataclass(frozen=True, slots=True)
class ParameterOverride:
    """One requested logical value, independent of its physical file."""

    key: str
    value: Any
    source: str = "user"


@dataclass(frozen=True, slots=True)
class ResolvedParameterChange:
    """A resolved change with both logical and physical provenance."""

    key: str
    previous_value: Any
    requested_value: Any
    source_path: str
    declared_type: str
    override_source: str


@dataclass(frozen=True, slots=True)
class ParameterResolutionPlan:
    """Validated, side-effect-free plan ready for preview or application."""

    changes: tuple[ResolvedParameterChange, ...]
    _document_updates: tuple[tuple[ParameterDocument, tuple[tuple[str, Any], ...]], ...]

    @property
    def is_empty(self) -> bool:
        return not self.changes

    def preview(self) -> tuple[str, ...]:
        """Return concise human-readable lines suitable for GUI confirmation."""
        return tuple(
            f"{change.key}: {change.previous_value!r} -> "
            f"{change.requested_value!r} [{change.source_path}]"
            for change in self.changes
        )

    def apply(self) -> tuple[ResolvedParameterChange, ...]:
        """Apply only after all target documents pass an external-edit preflight."""
        for document, _ in self._document_updates:
            if document.path.read_bytes() != document.original_bytes:
                raise RuntimeError(
                    f"{document.path.name} changed outside the Workbench; "
                    "reload before applying the Case parameter plan."
                )
        for document, updates in self._document_updates:
            document.save(dict(updates))
        return self.changes


def resolve_parameter_overrides(
    documents: Iterable[ParameterDocument],
    overrides: Iterable[ParameterOverride] | Mapping[str, Any],
) -> ParameterResolutionPlan:
    """Resolve logical keys to effective declarations without modifying files.

    ParameterModel defines the last-declaration-wins behavior. The resolver
    therefore targets the effective declaration even when duplicate keys occur
    across files. Missing keys are rejected instead of silently creating a new
    declaration in an arbitrary file.
    """
    document_list = tuple(documents)
    model = ParameterModel.from_documents(document_list)
    requested = (
        tuple(ParameterOverride(key, value) for key, value in overrides.items())
        if isinstance(overrides, Mapping)
        else tuple(overrides)
    )

    seen: set[str] = set()
    changes: list[ResolvedParameterChange] = []
    updates_by_path: dict[str, tuple[ParameterDocument, dict[str, Any]]] = {}

    documents_by_path = {str(document.path): document for document in document_list}
    for override in requested:
        normalized = override.key.casefold()
        if normalized in seen:
            raise ValueError(f"Duplicate override for parameter {override.key!r}.")
        seen.add(normalized)

        entry = model.get(override.key)
        if entry is None:
            raise KeyError(
                f"Parameter {override.key!r} is not declared in the Case manifest."
            )
        path_key = str(entry.source_path)
        document = documents_by_path.get(path_key)
        if document is None:
            raise RuntimeError(
                f"No loaded parameter document owns {override.key!r} at {path_key}."
            )
        bucket = updates_by_path.setdefault(path_key, (document, {}))
        bucket[1][entry.name] = override.value
        changes.append(ResolvedParameterChange(
            key=entry.name,
            previous_value=entry.value,
            requested_value=override.value,
            source_path=path_key,
            declared_type=entry.declared_type,
            override_source=override.source,
        ))

    ordered_updates = tuple(
        (document, tuple(values.items()))
        for document, values in updates_by_path.values()
    )
    return ParameterResolutionPlan(tuple(changes), ordered_updates)
