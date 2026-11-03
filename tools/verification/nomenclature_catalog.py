"""Source-backed CFD terminology cross-reference catalog."""
from __future__ import annotations
import json
from dataclasses import dataclass
from pathlib import Path
import re

_DEFAULT_DATA = Path(__file__).with_name("nomenclature") / "cfd_core.json"

@dataclass(frozen=True, slots=True)
class NomenclatureAlias:
    system: str
    name: str
    role: str
    relation: str
    applicability: str
    caveat: str
    source_ids: tuple[str, ...]

@dataclass(frozen=True, slots=True)
class NomenclatureEntry:
    id: str
    label: str
    symbol: str
    definition: str
    dimensions: tuple[tuple[str, int], ...]
    si_unit: str
    aliases: tuple[NomenclatureAlias, ...]

class NomenclatureCatalog:
    """Read-only semantic index; ambiguous names return multiple matches."""
    def __init__(self, entries, sources):
        self.entries = entries
        self.sources = sources

    @classmethod
    def load(cls, path=None):
        payload = json.loads(Path(path or _DEFAULT_DATA).read_text(encoding="utf-8"))
        if payload.get("schema_version") != 1:
            raise ValueError("Unsupported nomenclature schema version")
        sources = {item["id"]: item for item in payload["sources"]}
        if len(sources) != len(payload["sources"]):
            raise ValueError("Source IDs must be unique")
        entries = []
        seen = set()
        for item in payload["concepts"]:
            if item["id"] in seen:
                raise ValueError("Concept IDs must be unique")
            seen.add(item["id"])
            aliases = []
            for raw in item.get("aliases", []):
                missing = set(raw.get("source_ids", ())) - sources.keys()
                if missing:
                    raise ValueError(f"{item['id']}: unknown source IDs: {sorted(missing)}")
                if raw.get("relation") not in {"exact", "conditional", "related"}:
                    raise ValueError(f"{item['id']}: invalid alias relation")
                aliases.append(NomenclatureAlias(**{**raw, "source_ids": tuple(raw["source_ids"])}))
            entries.append(NomenclatureEntry(
                item["id"], item["label"], item["symbol"], item["definition"],
                tuple(sorted((str(k), int(v)) for k, v in item["dimensions"].items())),
                item["si_unit"], tuple(aliases)))
        return cls(tuple(entries), sources)

    def find_alias(self, name, *, system=None, role=None):
        key = _normalize(name)
        return tuple((entry, alias) for entry in self.entries for alias in entry.aliases
            if _normalize(alias.name) == key
            and (system is None or alias.system.casefold() == system.casefold())
            and (role is None or alias.role.casefold() == role.casefold()))

    def search(self, text):
        key = _normalize(text)
        if not key:
            return ()
        return tuple(entry for entry in self.entries if any(
            key in _normalize(value) for value in
            [entry.id, entry.label, entry.symbol, entry.definition,
             *(alias.name for alias in entry.aliases)]))

def _normalize(value):
    return re.sub(r"\s+", " ", str(value).strip()).casefold()

def validate_nomenclature_catalog(catalog):
    """Check dimensions and source-backed alias records; repeated names across concepts are valid."""
    issues = []
    for entry in catalog.entries:
        if not entry.dimensions:
            issues.append(f"{entry.id}: dimensions are missing")
        for alias in entry.aliases:
            if not alias.source_ids:
                issues.append(f"{entry.id}: alias {alias.name!r} has no source")
    return issues
