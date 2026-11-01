"""Shared style for Contract Workbench theory curves.

Use explicit ASCII axis names so labels remain readable when mathtext fonts
are incomplete; formula goes in the title/subtitle.
"""

from __future__ import annotations

from pathlib import Path


def apply_theory_style(
    ax,
    *,
    xlabel: str,
    ylabel: str,
    title: str,
    subtitle: str | None = None,
    legend: bool = True,
) -> None:
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    full_title = title if not subtitle else f"{title}\n{subtitle}"
    ax.set_title(full_title, fontsize=11)
    ax.grid(True, which="both", alpha=0.35)
    ax.tick_params(labelsize=10)
    if legend:
        ax.legend(loc="best", fontsize=9)
    for spine in ax.spines.values():
        spine.set_visible(True)


def save_figure(fig, out: Path, dpi: int = 140) -> None:
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=dpi, bbox_inches="tight", facecolor="white")
