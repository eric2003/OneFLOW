"""Explicit topology entities used before serializing legacy grid layouts."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True, slots=True)
class LayoutPoint:
    id: int
    x: float
    y: float
    z: float


@dataclass(frozen=True, slots=True)
class LayoutLine:
    id: int
    start_point: int
    end_point: int
    point_count: int
    distribution: Any
    boundary_type: int


@dataclass(frozen=True, slots=True)
class LayoutFace:
    id: int
    line_by_side: tuple[tuple[int, int], ...]


@dataclass(frozen=True, slots=True)
class LayoutBlock:
    id: int
    face_by_position: tuple[tuple[int, int], ...]


@dataclass(frozen=True, slots=True)
class LayoutTopology:
    points: tuple[LayoutPoint, ...]
    lines: tuple[LayoutLine, ...]
    faces: tuple[LayoutFace, ...]
    blocks: tuple[LayoutBlock, ...]

    def validate(self) -> None:
        point_ids = self._unique_ids(self.points, "point")
        line_ids = self._unique_ids(self.lines, "line")
        face_ids = self._unique_ids(self.faces, "face")
        self._unique_ids(self.blocks, "block")

        for line in self.lines:
            if line.start_point not in point_ids or line.end_point not in point_ids:
                raise ValueError(f"Line {line.id} references an unknown point.")
            if line.point_count < 2:
                raise ValueError(f"Line {line.id} must contain at least two points.")

        for face in self.faces:
            positions = [position for position, _ in face.line_by_side]
            if any(position < 1 for position in positions) or len(set(positions)) != len(positions):
                raise ValueError(f"Face {face.id} has invalid or duplicate side positions.")
            for _, line_id in face.line_by_side:
                if line_id not in line_ids:
                    raise ValueError(f"Face {face.id} references unknown line {line_id}.")

        for block in self.blocks:
            positions = [position for position, _ in block.face_by_position]
            if any(position < 1 for position in positions) or len(set(positions)) != len(positions):
                raise ValueError(f"Block {block.id} has invalid or duplicate face positions.")
            for _, face_id in block.face_by_position:
                if face_id not in face_ids:
                    raise ValueError(f"Block {block.id} references unknown face {face_id}.")

    @staticmethod
    def _unique_ids(items, entity_name: str) -> set[int]:
        identifiers = [item.id for item in items]
        if any(identifier < 1 for identifier in identifiers):
            raise ValueError(f"{entity_name.title()} identifiers must be positive.")
        if len(set(identifiers)) != len(identifiers):
            raise ValueError(f"Duplicate {entity_name} identifier.")
        return set(identifiers)


__all__ = [
    "LayoutBlock",
    "LayoutFace",
    "LayoutLine",
    "LayoutPoint",
    "LayoutTopology",
]
