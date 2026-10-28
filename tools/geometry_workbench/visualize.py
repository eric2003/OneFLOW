"""Optional topology sketch; editing stays in the host, not the plot library."""

from __future__ import annotations

from pathlib import Path
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


def save_topology_png(
    result: Dict[str, Any],
    path: Path,
    orphan_ids: Sequence[int] = (),
) -> Path:
    """Write a static PNG using matplotlib Agg (no GUI window)."""
    from pathlib import Path as PathType

    path = PathType(path)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    edges = _edges_from_result(result)
    zone_y: Dict[str, float] = {}
    glob_y: Dict[str, float] = {}
    for zone, gface, _ in edges:
        zone_y.setdefault(zone, float(len(zone_y)))
        glob_y.setdefault(gface, float(len(glob_y)))

    fig, ax = plt.subplots(figsize=(6.0, 3.5))
    for zone, gface, local in edges:
        ax.plot([0, 1], [zone_y[zone], glob_y[gface]], color="#888", lw=1)
        ax.text(0.5, 0.5 * (zone_y[zone] + glob_y[gface]), local, fontsize=8, ha="center")
    for name, y in zone_y.items():
        ax.scatter([0], [y], s=200, c="#4C78A8", zorder=3)
        ax.text(-0.05, y, name, ha="right", va="center")
    for name, y in glob_y.items():
        label = f"{name}*" if name[1:].isdigit() and int(name[1:]) in orphan_ids else name
        ax.scatter([1], [y], s=200, c="#F58518", zorder=3)
        ax.text(1.05, y, label, ha="left", va="center")
    ax.set_axis_off()
    ax.set_title("Interface dual-map (static export)")
    path = path.resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path
