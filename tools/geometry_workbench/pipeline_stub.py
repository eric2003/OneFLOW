"""Table-level scenario runners; native C++ binding remains separate."""

from __future__ import annotations

from copy import deepcopy
from typing import Any, Dict, List, Optional


def _run_registered_child_resolution(scenario: Dict[str, Any]) -> Dict[str, Any]:
    """Run the small lookup contract from input, independently of expected data."""
    input_data = scenario.get("input") or {}
    registered_faces = input_data.get("registeredFaces") or []
    candidates = input_data.get("candidateFaces") or []
    face_ids_by_nodes: dict[tuple[int, ...], int] = {}
    for face in registered_faces:
        node_key = tuple(sorted(face["nodes"]))
        if node_key in face_ids_by_nodes:
            raise ValueError("duplicate registered face identity")
        face_ids_by_nodes[node_key] = int(face["globalFaceId"])

    children: list[dict[str, int]] = []
    seen_global_ids: set[int] = set()
    unregistered: list[int] = []
    for candidate_index, nodes in enumerate(candidates):
        global_face_id = face_ids_by_nodes.get(tuple(sorted(nodes)))
        if global_face_id is None:
            unregistered.append(candidate_index)
            continue
        if global_face_id in seen_global_ids:
            raise ValueError("duplicate child face for parent")
        seen_global_ids.add(global_face_id)
        children.append({
            "globalFaceId": global_face_id,
            "candidateIndex": candidate_index,
        })

    return {
        "registeredChildren": children,
        "unregisteredCandidateIndices": unregistered,
        "registeredFaceCountAfterResolution": len(registered_faces),
    }


def run_pipeline(scenario: Dict[str, Any], backend: str = "stub") -> Dict[str, Any]:
    if backend == "native":
        raise NotImplementedError(
            "native backend not wired yet; bind RunSymmetricInterfaceSplitMappingPipeline"
        )
    if backend != "stub":
        raise ValueError(f"unknown backend: {backend}")
    if scenario.get("id") == "registered_child_resolution":
        return _run_registered_child_resolution(scenario)
    expected = scenario.get("expected")
    if not expected:
        raise ValueError("scenario has no expected tables for stub backend")
    # The legacy multi-zone fixture is a display fixture, not an algorithm runner.
    return deepcopy(expected)


def diff_tables(
    actual: Dict[str, Any],
    expected: Dict[str, Any],
    keys: Optional[List[str]] = None,
) -> List[str]:
    if keys is None:
        keys = sorted(set(expected) | set(actual))
    lines: List[str] = []
    for key in keys:
        if key not in expected:
            lines.append(f"{key}: unexpected key in actual")
            continue
        if key not in actual:
            lines.append(f"{key}: missing in actual")
            continue
        if actual[key] != expected[key]:
            lines.append(f"{key}: actual={actual[key]!r} expected={expected[key]!r}")
    return lines
