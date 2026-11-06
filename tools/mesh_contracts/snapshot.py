"""Validation for the format-neutral MeshSnapshot interchange contract.

This module validates interchange manifests; it is not the production C++
mesh container and does not claim to execute solver topology algorithms.
"""

from __future__ import annotations

import math
from typing import Any


class MeshSnapshotError(ValueError):
    """Raised when a MeshSnapshot violates its declared contract."""


_ELEMENT_DIMENSION = {
    "line2": 1,
    "triangle3": 2,
    "quad4": 2,
    "tetra4": 3,
    "hexa8": 3,
    "prism6": 3,
    "pyramid5": 3,
}

_ELEMENT_ARITY = {
    "line2": 2,
    "triangle3": 3,
    "quad4": 4,
    "tetra4": 4,
    "hexa8": 8,
    "prism6": 6,
    "pyramid5": 5,
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise MeshSnapshotError(message)


def validate_snapshot(snapshot: dict[str, Any]) -> None:
    """Validate IDs, coordinates, cell connectivity, groups and provenance."""
    _require(isinstance(snapshot, dict), "snapshot must be an object")
    _require(snapshot.get("schema") == "oneflow.mesh-snapshot", "unsupported snapshot schema")
    _require(snapshot.get("schemaVersion") == 1, "unsupported snapshot schemaVersion")
    dimension = snapshot.get("dimension")
    _require(type(dimension) is int and 1 <= dimension <= 3, "dimension must be an integer from 1 to 3")

    nodes = snapshot.get("nodes")
    cells = snapshot.get("cells")
    faces = snapshot.get("faces")
    groups = snapshot.get("groups")
    _require(isinstance(nodes, list) and nodes, "nodes must be a non-empty array")
    _require(isinstance(cells, list), "cells must be an array")
    _require(isinstance(faces, list), "faces must be an array")
    _require(isinstance(groups, list), "groups must be an array")

    node_ids: set[int] = set()
    for index, node in enumerate(nodes):
        _require(isinstance(node, dict), f"nodes[{index}] must be an object")
        node_id = node.get("id")
        _require(type(node_id) is int and node_id > 0, f"nodes[{index}].id must be a positive integer")
        _require(node_id not in node_ids, f"duplicate node id: {node_id}")
        node_ids.add(node_id)
        coordinates = node.get("coordinates")
        _require(isinstance(coordinates, list) and len(coordinates) == dimension,
                 f"nodes[{index}].coordinates must contain {dimension} values")
        _require(all(type(value) in (int, float) and math.isfinite(value) for value in coordinates),
                 f"nodes[{index}].coordinates must be finite numbers")

    cell_ids: set[int] = set()
    for index, cell in enumerate(cells):
        _require(isinstance(cell, dict), f"cells[{index}] must be an object")
        cell_id = cell.get("id")
        _require(type(cell_id) is int and cell_id > 0, f"cells[{index}].id must be a positive integer")
        _require(cell_id not in cell_ids, f"duplicate cell id: {cell_id}")
        cell_ids.add(cell_id)
        cell_type = cell.get("type")
        _require(isinstance(cell_type, str) and cell_type in _ELEMENT_ARITY,
                 f"cells[{index}].type is unsupported: {cell_type}")
        _require(_ELEMENT_DIMENSION[cell_type] == dimension,
                 f"cells[{index}].type dimension does not match mesh dimension")
        node_refs = cell.get("nodeIds")
        _require(isinstance(node_refs, list) and len(node_refs) == _ELEMENT_ARITY[cell_type],
                 f"cells[{index}].nodeIds has invalid arity for {cell_type}")
        _require(all(type(ref) is int and ref in node_ids for ref in node_refs),
                 f"cells[{index}].nodeIds contains an unknown node ID")
        _require(len(set(node_refs)) == len(node_refs), f"cells[{index}].nodeIds contains duplicate node IDs")

    face_ids: set[int] = set()
    for index, face in enumerate(faces):
        _require(isinstance(face, dict), f"faces[{index}] must be an object")
        face_id = face.get("id")
        _require(type(face_id) is int and face_id > 0, f"faces[{index}].id must be a positive integer")
        _require(face_id not in face_ids, f"duplicate face id: {face_id}")
        face_ids.add(face_id)
        node_refs = face.get("nodeIds")
        _require(isinstance(node_refs, list) and len(node_refs) >= 2, f"faces[{index}].nodeIds must contain at least two node IDs")
        _require(all(type(ref) is int and ref in node_ids for ref in node_refs), f"faces[{index}].nodeIds contains an unknown node ID")
        _require(len(set(node_refs)) == len(node_refs), f"faces[{index}].nodeIds contains duplicate node IDs")
        owner = face.get("ownerCellId")
        neighbour = face.get("neighbourCellId")
        _require(type(owner) is int and owner in cell_ids, f"faces[{index}].ownerCellId must reference an existing cell")
        _require(neighbour is None or (type(neighbour) is int and neighbour in cell_ids), f"faces[{index}].neighbourCellId must be null or reference an existing cell")
        _require(neighbour != owner, f"faces[{index}] owner and neighbour must differ")

    group_ids: set[str] = set()
    for index, group in enumerate(groups):
        _require(isinstance(group, dict), f"groups[{index}] must be an object")
        group_id = group.get("id")
        _require(isinstance(group_id, str) and group_id.strip(), f"groups[{index}].id must be non-empty")
        _require(group_id not in group_ids, f"duplicate group id: {group_id}")
        group_ids.add(group_id)
        _require(group.get("kind") in ("boundary", "region", "interface"),
                 f"groups[{index}].kind must be boundary, region or interface")
        group_dimension = group.get("dimension")
        _require(type(group_dimension) is int and 0 <= group_dimension <= dimension,
                 f"groups[{index}].dimension is outside the mesh dimension")
        if group["kind"] == "region":
            _require(group_dimension == dimension,
                     f"groups[{index}] region dimension must equal mesh dimension")
        else:
            _require(group_dimension == dimension - 1,
                     f"groups[{index}] boundary/interface dimension must be mesh dimension minus one")
        member_kind = group.get("memberKind")
        _require(member_kind in ("cell", "face"), f"groups[{index}].memberKind must be cell or face")
        if group["kind"] == "region":
            _require(member_kind == "cell", f"groups[{index}] region members must be cells")
        else:
            _require(member_kind == "face", f"groups[{index}] boundary/interface members must be faces")
        members = group.get("memberIds")
        _require(isinstance(members, list), f"groups[{index}].memberIds must be an array")
        _require(all(type(member) is int for member in members), f"groups[{index}].memberIds must contain integer IDs")
        _require(len(set(members)) == len(members), f"groups[{index}].memberIds contains duplicates")
        available_ids = cell_ids if member_kind == "cell" else face_ids
        _require(set(members).issubset(available_ids), f"groups[{index}].memberIds contains an unknown {member_kind} ID")
        _require(isinstance(group.get("name"), str) and group["name"].strip(),
                 f"groups[{index}].name must be non-empty")

    provenance = snapshot.get("provenance")
    _require(isinstance(provenance, dict), "provenance must be an object")
    _require(isinstance(provenance.get("sourceFormat"), str) and provenance["sourceFormat"].strip(),
             "provenance.sourceFormat must be non-empty")
