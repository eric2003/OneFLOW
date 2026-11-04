"""Versioned, declarative Case profiles built on logical parameter overrides.

Profiles describe intent and provenance. They do not imply that arbitrary
combinations are numerically valid or that a profile has passed runtime checks.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from .parameter_document import ParameterDocument
from .parameter_resolution import (
    ParameterOverride,
    ParameterResolutionPlan,
    resolve_parameter_overrides,
)


@dataclass(frozen=True, slots=True)
class CaseProfile:
    profile_id: str
    title: str
    source_case: str
    run_purpose: str
    maturity: str
    description: str
    overrides: tuple[ParameterOverride, ...]
    limitations: tuple[str, ...]


def load_case_profile(path: str | Path) -> CaseProfile:
    """Load and validate one versioned profile definition."""
    source = Path(path)
    try:
        payload = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ValueError(f"Cannot read Case profile {source}: {error}") from error
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        raise ValueError("Unsupported Case profile schema; expected schema_version 1.")

    required = ("profile_id", "title", "source_case", "run_purpose", "maturity", "description")
    values = {key: payload.get(key) for key in required}
    for key, value in values.items():
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"Case profile field {key!r} must be a non-empty string.")
    if not re.fullmatch(r"[a-z0-9][a-z0-9._-]*", values["profile_id"]):
        raise ValueError("Case profile profile_id must be a lowercase identifier.")
    if values["maturity"] not in {"proposed", "case_derived", "runtime_checked", "numerically_verified"}:
        raise ValueError(f"Unknown Case profile maturity: {values['maturity']}")

    raw_overrides = payload.get("overrides", [])
    if not isinstance(raw_overrides, list):
        raise ValueError("Case profile overrides must be a list.")
    overrides = []
    seen: set[str] = set()
    for item in raw_overrides:
        if not isinstance(item, dict):
            raise ValueError("Each Case profile override must be an object.")
        key, value = item.get("key"), item.get("value")
        if not isinstance(key, str) or not key.strip():
            raise ValueError("Each Case profile override requires a logical key.")
        if isinstance(value, (dict, list)) or value is None or isinstance(value, bool):
            raise ValueError(f"Unsupported scalar override for {key!r}.")
        normalized = key.casefold()
        if normalized in seen:
            raise ValueError(f"Duplicate Case profile override: {key!r}.")
        seen.add(normalized)
        override_source = item.get("source", values["profile_id"])
        if not isinstance(override_source, str) or not override_source.strip():
            raise ValueError(f"Override source for {key!r} must be a non-empty string.")
        overrides.append(ParameterOverride(key, value, override_source))

    limitations = payload.get("limitations", [])
    if not isinstance(limitations, list) or not all(
        isinstance(item, str) and item.strip() for item in limitations
    ):
        raise ValueError("Case profile limitations must be a list of non-empty strings.")

    return CaseProfile(
        profile_id=values["profile_id"],
        title=values["title"],
        source_case=values["source_case"].replace("\\", "/"),
        run_purpose=values["run_purpose"],
        maturity=values["maturity"],
        description=values["description"],
        overrides=tuple(overrides),
        limitations=tuple(limitations),
    )


def resolve_case_profile(
    profile: CaseProfile,
    documents: Iterable[ParameterDocument],
) -> ParameterResolutionPlan:
    """Resolve profile values without writing files or assuming file names."""
    return resolve_parameter_overrides(documents, profile.overrides)
