"""Optional topology sketch; editing stays in the host, not the plot library."""

from __future__ import annotations

from typing import Any, Dict, List, Sequence, Tuple


def _edges_from_result(result: Dict[str, Any]) -> List[Tuple[str, str, str]]:
    g_zid = result.get("gI2Zid") or []
    g_local = result.get("g2l") or []
    edges: List[Tuple[str, str, str]] = []
    for gid, zones in enumerate(g_zid):
        locals_row = g_local[gid] if gid < len(g_local) else []
        for slot, zid in enumerate(zones):
            local = locals_row[slot] if slot < len(locals_row) else "?"
            edges.append((f"Z{zid}", f"G{gid}", str(local)))
    return edges


def render_ascii(result: Dict[str, Any], orphan_ids: Sequence[int] = ()) -> str:
    lines = ["zone --local--> global face"]
    for zone, gface, local in _edges_from_result(result):
        mark = ""
        if gface[1:].isdigit() and int(gface[1:]) in orphan_ids:
            mark = " (orphan-row)"
        lines.append(f"  {zone} local {local} -> {gface}{mark}")
    l2g = result.get("l2g") or []
    lines.append("l2g (zone -> global ids):")
    for zid, row in enumerate(l2g):
        lines.append(f"  Z{zid}: {row}")
    return "\n".join(lines)
