"""User-owned import and registration of non-STL geometry assets."""

from __future__ import annotations

import hashlib
import json
import re
import shutil
from pathlib import Path

from .user_data_paths import UserDataPaths

_FORMATS = {
    ".step": ("cad_brep", "STEP", "imported"),
    ".stp": ("cad_brep", "STEP", "imported"),
    ".iges": ("cad_brep", "IGES", "imported"),
    ".igs": ("cad_brep", "IGES", "imported"),
    ".brep": ("cad_brep", "BREP", "imported"),
    ".obj": ("surface_mesh", "OBJ", "imported"),
    ".ply": ("surface_mesh", "PLY", "imported"),
    ".xyz": ("point_cloud", "XYZ point cloud", "measured"),
    ".pts": ("point_cloud", "PTS point cloud", "measured"),
    ".pcd": ("point_cloud", "PCD point cloud", "measured"),
    ".implicit.json": ("implicit_field", "OneFLOW implicit definition JSON", "implicit"),
}


def supported_geometry_imports() -> tuple[str, ...]:
    """Return metadata-registration suffixes, not parser capabilities."""
    return tuple(sorted(_FORMATS))


def import_geometry_asset(source_file: Path) -> dict:
    """Copy a supported non-STL file into user data and register its provenance."""
    source = Path(source_file).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(f"Geometry source file does not exist: {source}")
    suffix = next((item for item in _FORMATS if source.name.lower().endswith(item)), None)
    if suffix is None:
        raise ValueError(
            "Unsupported geometry import format. Supported: "
            + ", ".join(supported_geometry_imports())
        )
    role, format_name, definition_kind = _FORMATS[suffix]
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    safe_suffix = re.sub(r"[^a-z0-9]+", "_", suffix.lstrip("."))
    asset_id = f"geometry.imported.{safe_suffix}.{digest[:20]}"
    user_root = UserDataPaths.discover().root
    workspace = user_root / "geometry-library" / "assets" / asset_id
    manifest_path = workspace / "asset.json"
    if manifest_path.is_file():
        return json.loads(manifest_path.read_text(encoding="utf-8"))
    workspace.mkdir(parents=True, exist_ok=True)
    destination = workspace / ("source" + suffix)
    shutil.copy2(source, destination)
    if hashlib.sha256(destination.read_bytes()).hexdigest() != digest:
        destination.unlink(missing_ok=True)
        raise OSError("Geometry import copy failed its SHA-256 verification.")
    asset = {
        "schema": "oneflow.geometry-asset",
        "version": 1,
        "id": asset_id,
        "name": source.name[:-len(suffix)],
        "description": f"User-imported {format_name} geometry representation.",
        "dimension": 3,
        "tags": ["imported", format_name.lower().replace(" ", "-")],
        "definition": {
            "kind": definition_kind,
            "notes": "Imported representation registered by content hash; geometric parsing and meshing adapters are not implied.",
        },
        "provenance": {
            "kind": "user_supplied",
            "notes": "User-selected source file; original authoring application and upstream history are not independently verified.",
            "source_file_name": source.name,
            "source_sha256": digest,
        },
        "representations": [{
            "role": role,
            "format": format_name,
            "path": destination.relative_to(user_root).as_posix(),
            "status": "external",
            "sha256": digest,
            "notes": "Stored in the user geometry library. Import/preview/meshing support depends on a future format adapter.",
        }],
        "metadata": {"registration_status": "registered_not_parsed"},
    }
    manifest_path.write_text(
        json.dumps(asset, indent=2, sort_keys=True) + "\n",
        encoding="utf-8", newline="\n",
    )
    return asset


def load_user_geometry_assets() -> list[dict]:
    """Load valid user asset manifests from the managed geometry library."""
    root = UserDataPaths.discover().root / "geometry-library" / "assets"
    if not root.is_dir():
        return []
    assets = []
    for manifest in sorted(root.glob("*/asset.json")):
        try:
            value = json.loads(manifest.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if (isinstance(value, dict)
                and value.get("schema") == "oneflow.geometry-asset"
                and value.get("id")):
            assets.append(value)
    return assets
