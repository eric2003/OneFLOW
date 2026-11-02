"""Shared style for Contract Workbench theory curves."""
from __future__ import annotations
from pathlib import Path

def apply_theory_style(ax, *, xlabel, ylabel, title, subtitle=None, legend=True):
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(f"{title}\n{subtitle}" if subtitle else title, fontsize=12)
    ax.grid(True, which="both", alpha=0.35)
    ax.tick_params(labelsize=10)
    if legend:
        ax.legend(loc="best", fontsize=10)

def save_figure(fig, out, dpi=140):
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=dpi, bbox_inches="tight", facecolor="white")
