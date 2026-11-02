"""Formula preview for Contract Workbench using matplotlib.

Rendering order:
1. Full LaTeX via matplotlib usetex (if a TeX engine is available).
2. High-quality matplotlib mathtext (no external TeX required).
3. Plain monospace fallback for unsupported lines.

mathtext is not full LaTeX; sanitize common macros before mathtext.
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path


def _has_tex_engine() -> bool:
    return any(shutil.which(name) for name in ("latex", "pdflatex", "xelatex", "lualatex"))


def sanitize_for_mathtext(latex_line: str) -> str:
    s = latex_line.strip().rstrip("\\")
    if not s:
        return s
    replacements = [
        (r"\\mathrm\{([^}]*)\}", r"\1"),
        (r"\\operatorname\{([^}]*)\}", r"\1"),
        (r"\\text\{([^}]*)\}", r"\\mathrm{\1}"),
        (r"\\left", ""),
        (r"\\right", ""),
        (r"\\Bigg?", ""),
        (r"\\bigg?", ""),
        (r"\\,", r"\,"),
        (r"\\;", r"\;"),
        (r"\\quad", r"\quad"),
        (r"\\dfrac", r"\frac"),
        (r"\\tfrac", r"\frac"),
        (r"\\longrightarrow", r"\rightarrow"),
        (r"\\leftarrow", r"\leftarrow"),
        (r"\\,", " "),
    ]
    for pat, rep in replacements:
        s = re.sub(pat, rep, s)
    # Display-style fractions look better in preview
    if r"\frac" in s and not s.lstrip().startswith(r"\displaystyle"):
        s = r"\displaystyle " + s
    return s


def _prepare_lines(latex_body: str) -> list[str]:
    lines = []
    for ln in latex_body.splitlines():
        ln = ln.strip()
        if not ln:
            continue
        lines.append(ln)
    return lines or [r"\text{(empty formula)}"]


def _render_with_usetex(lines: list[str], out_png: Path) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # One stacked equation block
    body = " \\\\\n".join(lines)
    tex = (
        r"\begin{align*}" + "\n" + body + "\n" + r"\end{align*}"
    )

    plt.rcParams.update(
        {
            "text.usetex": True,
            "font.family": "serif",
            "font.size": 14,
            "figure.dpi": 160,
        }
    )
    fig = plt.figure(figsize=(9.0, max(1.6, 0.55 * len(lines) + 0.9)))
    fig.patch.set_facecolor("white")
    ax = fig.add_axes([0.02, 0.05, 0.96, 0.9])
    ax.set_axis_off()
    ax.set_facecolor("white")
    try:
        ax.text(
            0.5,
            0.5,
            f"${tex}$" if "align" not in tex else f"${body}$",
            ha="center",
            va="center",
            fontsize=16,
            color="#111111",
            transform=ax.transAxes,
        )
        # usetex often prefers non-align in single text; use line-by-line
    except Exception:
        plt.close(fig)
        raise

    # Prefer reliable line-by-line usetex
    plt.close(fig)
    fig = plt.figure(figsize=(9.0, max(1.8, 0.7 * len(lines) + 1.0)))
    fig.patch.set_facecolor("white")
    ax = fig.add_axes([0.04, 0.08, 0.92, 0.84])
    ax.set_axis_off()
    y = 0.9
    for ln in lines:
        ax.text(
            0.5,
            y,
            f"${ln}$",
            ha="center",
            va="top",
            fontsize=16,
            color="#111111",
            transform=ax.transAxes,
        )
        y -= 0.85 / max(len(lines), 1)
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=160, facecolor="white", bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)
    # reset usetex so other plots are not forced
    plt.rcParams.update({"text.usetex": False})
    return out_png


def _render_with_mathtext(lines: list[str], out_png: Path, *, dpi: int = 160, fontsize: float = 18) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "text.usetex": False,
            "mathtext.fontset": "cm",  # Computer Modern-like, closer to LaTeX
            "mathtext.default": "it",
            "font.size": 14,
            "figure.dpi": 160,
        }
    )

    n = len(lines)
    fig_h = max(2.0, 0.85 * n + 1.0)
    fig = plt.figure(figsize=(9.5, fig_h))
    fig.patch.set_facecolor("white")
    ax = fig.add_axes([0.05, 0.08, 0.9, 0.84])
    ax.set_axis_off()
    ax.set_facecolor("white")

    # Light card background
    ax.add_patch(
        plt.Rectangle(
            (0.0, 0.0),
            1.0,
            1.0,
            transform=ax.transAxes,
            facecolor="#fafafa",
            edgecolor="#e0e0e0",
            linewidth=1.0,
            zorder=0,
        )
    )

    y = 0.88
    dy = 0.75 / max(n, 1)
    for ln in lines:
        candidates = [ln, sanitize_for_mathtext(ln)]
        shown = False
        for cand in candidates:
            if not cand:
                continue
            try:
                ax.text(
                    0.5,
                    y,
                    f"${cand}$",
                    ha="center",
                    va="top",
                    fontsize=fontsize,
                    color="#111111",
                    transform=ax.transAxes,
                    zorder=1,
                )
                shown = True
                break
            except Exception:
                continue
        if not shown:
            ax.text(
                0.5,
                y,
                ln,
                ha="center",
                va="top",
                fontsize=13,
                color="#333333",
                family="monospace",
                transform=ax.transAxes,
                zorder=1,
            )
        y -= dy

    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        out_png,
        dpi=dpi,
        facecolor="white",
        bbox_inches="tight",
        pad_inches=0.3,
    )
    plt.close(fig)
    return out_png


def render_latex_preview_png(
    latex_body: str,
    out_png: Path,
    *,
    dpi: int = 160,
    fontsize: float = 18,
) -> Path:
    """Render formula lines to PNG. Prefer TeX, else mathtext, else plain.

    dpi/fontsize control sharpness for zoomed or popup views.
    """
    lines = _prepare_lines(latex_body)
    out_png = Path(out_png)

    if _has_tex_engine():
        try:
            return _render_with_usetex(lines, out_png)
        except Exception:
            pass

    try:
        return _render_with_mathtext(lines, out_png, dpi=dpi, fontsize=fontsize)
    except Exception:
        # Ultimate fallback: plain text figure
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, max(1.5, 0.4 * len(lines) + 0.8)))
        ax.axis("off")
        ax.text(
            0.02,
            0.98,
            "\n".join(lines),
            va="top",
            family="monospace",
            fontsize=12,
            transform=ax.transAxes,
        )
        out_png.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_png, dpi=140, facecolor="white", bbox_inches="tight")
        plt.close(fig)
        return out_png
