"""Catalog of intermediate theory curves for the Contract Workbench GUI."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class TheoryCurve:
    id: str
    title: str
    domain: str
    description: str
    csv_name: str
    plot_module: str
    proves: str
    does_not_prove: str


def verification_root(repository_root: Path) -> Path:
    return (repository_root / "tools" / "verification").resolve()


def list_theory_curves() -> list[TheoryCurve]:
    return [
        TheoryCurve(
            id="sa_fw_r",
            title="SA fw(r) wall destruction",
            domain="turbulence",
            description="Destruction multiplier vs ratio r; log-layer point (1,1).",
            csv_name="sa_fw_r_theory.csv",
            plot_module="plot_sa_fw",
            proves="Closed-form fw(r) table and figure.",
            does_not_prove="Full CalcSrcSa residual or skin friction.",
        ),
        TheoryCurve(
            id="sa_fv1_chi",
            title="SA fv1(chi) eddy-viscosity damping",
            domain="turbulence",
            description="Maps working variable ratio chi to damping of mu_t.",
            csv_name="sa_fv1_chi_theory.csv",
            plot_module="plot_sa_fv1",
            proves="Closed-form fv1(chi) table and figure.",
            does_not_prove="Production/destruction balance or experiments.",
        ),
        TheoryCurve(
            id="sa_fv2_chi",
            title="SA fv2(chi) strain correction",
            domain="turbulence",
            description="Couples chi into modified strain used by production.",
            csv_name="sa_fv2_chi_theory.csv",
            plot_module="plot_sa_fv2",
            proves="Closed-form fv2(chi) table and figure.",
            does_not_prove="Full modified-strain residual under flow state.",
        ),
    ]


def csv_path(repository_root: Path, curve: TheoryCurve) -> Path:
    return verification_root(repository_root) / "data" / curve.csv_name


def render_png(repository_root: Path, curve: TheoryCurve, out_png: Path) -> Path:
    import runpy
    import sys

    ver_root = verification_root(repository_root)
    script = ver_root / f"{curve.plot_module}.py"
    if not script.is_file():
        raise FileNotFoundError(f"missing plot script: {script}")
    csv = csv_path(repository_root, curve)
    if not csv.is_file():
        raise FileNotFoundError(f"missing theory CSV: {csv}")
    old_argv = sys.argv[:]
    try:
        sys.argv = [str(script), "--csv", str(csv), "--out", str(out_png)]
        runpy.run_path(str(script), run_name="__main__")
    finally:
        sys.argv = old_argv
    if not out_png.is_file():
        raise RuntimeError(f"plot script did not write {out_png}")
    return out_png
