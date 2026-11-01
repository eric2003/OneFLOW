"""Best-effort formula preview (mathtext is not full LaTeX)."""
from __future__ import annotations
import re
from pathlib import Path

def sanitize_for_mathtext(latex_line: str) -> str:
    s = latex_line.strip()
    for pat, rep in [
        (r"\\mathrm\{([^}]*)\}", r"\1"),
        (r"\\text\{([^}]*)\}", r"\1"),
        (r"\\left", ""), (r"\\right", ""),
        (r"\\Bigg?", ""), (r"\\bigg?", ""),
        (r"\\,", " "), (r"\\quad", " "),
        (r"\\dfrac", r"\\frac"), (r"\\tfrac", r"\\frac"),
    ]:
        s = re.sub(pat, rep, s)
    return s

def render_latex_preview_png(latex_body: str, out_png: Path) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    lines = [ln.strip() for ln in latex_body.splitlines() if ln.strip()] or ["(empty)"]
    fig, ax = plt.subplots(figsize=(8.0, max(1.8, 0.55 * len(lines) + 0.8)))
    ax.axis("off")
    y = 0.92
    fallback = False
    for ln in lines[:10]:
        shown = False
        for cand in (ln, sanitize_for_mathtext(ln)):
            if not cand:
                continue
            try:
                ax.text(0.02, y, f"${cand}$", fontsize=13, transform=ax.transAxes, va="top")
                shown = True
                break
            except Exception:
                continue
        if not shown:
            fallback = True
            ax.text(0.02, y, ln, fontsize=11, transform=ax.transAxes, va="top", family="monospace")
        y -= 0.09
    if fallback:
        ax.text(0.02, 0.04, "Note: some lines as plain text (mathtext != full LaTeX).",
                fontsize=8, transform=ax.transAxes, color="#666")
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=120, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_png
