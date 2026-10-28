"""Optional topology sketch for dual-map results.

Interactivity tiers (pick what you have installed):

1. ASCII always works (no deps).
2. matplotlib: zoom/pan toolbar only; weak for click-edit.
3. plotly: browser zoom/hover/pan; better for exploration, still not a full GUI.

Editable tables belong in the workbench host (CLI/REPL/future panel), not in
the plot library.
"""

from __future__ import annotations

from typing import Any, Dict, List, Sequence, Tuple


def _edges_from_result(result: Dict[str, Any]) -> List[Tuple[str, str, str]]:
    """Build (zone_label, global_face_label, local_id) edges from gI2Zid/g2l."""
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
        mark = " (orphan-row)" if gface[1:].isdigit() and int(gface[1:]) in orphan_ids else ""
        lines.append(f"  {zone} local {local} -> {gface}{mark}")
    l2g = result.get("l2g") or []
    lines.append("l2g (zone -> global ids):")
    for zid, row in enumerate(l2g):
        lines.append(f"  Z{zid}: {row}")
    return "\n".join(lines)


def try_render_plotly(result: Dict[str, Any], orphan_ids: Sequence[int] = ()) -> bool:
    try:
        import plotly.graph_objects as go
    except ImportError:
        return False

    edges = _edges_from_result(result)
    # Simple bipartite layout: zones on x=0, globals on x=1
    zone_y: Dict[str, float] = {}
    glob_y: Dict[str, float] = {}
    for zone, gface, _ in edges:
        zone_y.setdefault(zone, float(len(zone_y)))
        glob_y.setdefault(gface, float(len(glob_y)))

    edge_x: List[Any] = []
    edge_y: List[Any] = []
    for zone, gface, _ in edges:
        edge_x += [0.0, 1.0, None]
        edge_y += [zone_y[zone], glob_y[gface], None]

    node_x = [0.0] * len(zone_y) + [1.0] * len(glob_y)
    node_y = list(zone_y.values()) + list(glob_y.values())
    node_text = list(zone_y.keys()) + [
        f"{g}*" if g[1:].isdigit() and int(g[1:]) in orphan_ids else g for g in glob_y.keys()
    ]

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=edge_x, y=edge_y, mode="lines", line=dict(width=1, color="#888"), hoverinfo="none"))
    fig.add_trace(
        go.Scatter(
            x=node_x,
            y=node_y,
            mode="markers+text",
            text=node_text,
            textposition="top center",
            marker=dict(size=18, color="#4C78A8"),
            hoverinfo="text",
        )
    )
    fig.update_layout(
        title="Interface dual-map (zone <-> global face)",
        showlegend=False,
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        margin=dict(l=20, r=20, t=40, b=20),
    )
    fig.show()
    return True


def try_render_matplotlib(result: Dict[str, Any], orphan_ids: Sequence[int] = ()) -> bool:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return False

    edges = _edges_from_result(result)
    zone_y: Dict[str, float] = {}
    glob_y: Dict[str, float] = {}
    for zone, gface, _ in edges:
        zone_y.setdefault(zone, float(len(zone_y)))
        glob_y.setdefault(gface, float(len(glob_y)))

    fig, ax = plt.subplots(figsize=(6, 3.5))
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
    ax.set_title("Interface dual-map (matplotlib: zoom/pan only)")
    plt.tight_layout()
    plt.show()
    return True


def render(result: Dict[str, Any], orphan_ids: Sequence[int] = (), prefer: str = "auto") -> str:
    """Return ASCII always; optionally open an interactive window/browser."""
    text = render_ascii(result, orphan_ids)
    if prefer == "none":
        return text
    if prefer in ("auto", "plotly") and try_render_plotly(result, orphan_ids):
        return text
    if prefer in ("auto", "matplotlib") and try_render_matplotlib(result, orphan_ids):
        return text
    return text
