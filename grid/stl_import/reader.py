"""STL file-format reader. Parsing preserves facet records and does not infer topology."""
from __future__ import annotations
import math
import re
import struct
from pathlib import Path

Point = tuple[float, float, float]
Triangle = tuple[Point, Point, Point]

class StlReadError(ValueError):
    """Invalid or unsupported STL file representation."""

def _is_binary_stl_data(data: bytes) -> bool:
    if len(data) < 84:
        return False
    count = struct.unpack_from("<I", data, 80)[0]
    return 84 + count * 50 == len(data)

def is_binary_stl(path: Path) -> bool:
    return _is_binary_stl_data(path.read_bytes())

def read_stl(path: Path) -> list[Triangle]:
    data = path.read_bytes()
    if _is_binary_stl_data(data):
        count = struct.unpack_from("<I", data, 80)[0]
        result, offset = [], 84
        for _ in range(count):
            row = struct.unpack_from("<12fH", data, offset)
            result.append(tuple(tuple(float(v) for v in row[i:i+3]) for i in (3, 6, 9)))
            offset += 50
        if not all(math.isfinite(v) for tri in result for p in tri for v in p):
            raise StlReadError("binary STL contains non-finite coordinates")
        return result
    try:
        content = data.decode("ascii")
    except UnicodeDecodeError as exc:
        raise StlReadError("input is neither binary STL nor ASCII STL") from exc
    rows = re.findall(r"\bvertex\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)", content, re.I)
    if not rows or len(rows) % 3:
        raise StlReadError("ASCII STL must contain complete groups of three vertices")
    points = [tuple(float(v) for v in row) for row in rows]
    if not all(math.isfinite(v) for p in points for v in p):
        raise StlReadError("ASCII STL contains non-finite coordinates")
    return [tuple(points[i:i+3]) for i in range(0, len(points), 3)]
