"""Catalog of intermediate theory knowledge entries for Contract Workbench.

Paths are always relative to the repository root (portable across OS/users).
Entries are documentation + theory tables, not production solvers.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class TheoryEntry:
    """One browsable unit: formula + role + refs + optional theory curve."""

    id: str
    # Taxonomy (coarse -> fine), e.g. turbulence / SA / closure / fw
    family: str
    genus: str
    section: str
    title: str
    # Short labels for tree and search
    role_tags: tuple[str, ...]
    formula: str
    principle: str
    background: str
    references: tuple[str, ...]
    # Repo-relative POSIX-style paths (use Path with repository_root)
    csv_relpath: str | None
    plot_module: str | None
    pure_apis: tuple[str, ...]
    doc_relpath: str | None
    proves: str
    does_not_prove: str


def list_theory_entries() -> list[TheoryEntry]:
    return [
        TheoryEntry(
            id="sa_fv1_chi",
            family="turbulence",
            genus="SA",
            section="closure / eddy viscosity",
            title="fv1(chi) eddy-viscosity damping",
            role_tags=("closure", "eddy-viscosity", "near-wall damping"),
            formula=(
                "chi = nu_tilde / nu\n"
                "fv1(chi) = chi^3 / (chi^3 + cv1^3)\n"
                "mu_t = rho * nu_tilde * fv1\n"
                "standard: cv1 = 7.1"
            ),
            principle=(
                "SA transports a working variable nu_tilde, not mu_t directly. "
                "fv1 maps chi into a damping factor so mu_t vanishes at the wall "
                "(chi->0 => fv1->0) and approaches undamped core behavior "
                "(chi->inf => fv1->1)."
            ),
            background=(
                "Without fv1, a naive mu_t = rho*nu_tilde would over-predict "
                "near-wall turbulent stress. fv1 is calibrated with the rest of "
                "the SA closure (fv2, modified strain, fw)."
            ),
            references=(
                "Spalart, P. R. and Allmaras, S. R. (1992/1994). "
                "A one-equation turbulence model for aerodynamic flows.",
                "NASA Turbulence Modeling Resource — Spalart-Allmaras "
                "https://turbmodels.larc.nasa.gov/spalart.html",
                "OneFLOW: doc/source/numerical/sa_eddy_viscosity.rst",
            ),
            csv_relpath="tools/verification/data/sa_fv1_chi_theory.csv",
            plot_module="plot_sa_fv1",
            pure_apis=("SaChi", "SaFv1", "SaEddyViscosity"),
            doc_relpath="doc/source/numerical/sa_eddy_viscosity.rst",
            proves="Closed-form fv1(chi) vs pure API SaFv1 on a chi grid.",
            does_not_prove="Full residual, wall distance, or experimental Cf.",
        ),
        TheoryEntry(
            id="sa_fv2_chi",
            family="turbulence",
            genus="SA",
            section="closure / production strain",
            title="fv2(chi) strain correction",
            role_tags=("closure", "production", "modified strain"),
            formula=(
                "fv2 = 1 - chi / (1 + chi * fv1)\n"
                "rs ~ nu_tilde / (kappa^2 * d^2 * Re)\n"
                "S_bar = rs * fv2\n"
                "S_tilde = Omega + S_bar   (or limited branch)"
            ),
            principle=(
                "Production uses a modified strain S_tilde, not raw vorticity "
                "alone. fv2 feeds chi back into that strain so near-wall "
                "damping stays consistent with fv1."
            ),
            background=(
                "S_bar vanishes with nu_tilde at the wall and scales with "
                "1/d^2. The limited branch (c2=0.7, c3=0.9) prevents "
                "non-physical S_tilde when S_bar is strongly negative."
            ),
            references=(
                "Spalart-Allmaras original model and later implementation notes.",
                "NASA TMR SA page.",
                "OneFLOW: doc/source/numerical/sa_closure.rst",
            ),
            csv_relpath="tools/verification/data/sa_fv2_chi_theory.csv",
            plot_module="plot_sa_fv2",
            pure_apis=("SaFv2", "SaSBar", "SaModifiedStrainRate"),
            doc_relpath="doc/source/numerical/sa_closure.rst",
            proves="Closed-form fv2 and modified-strain branches vs pure APIs.",
            does_not_prove="Assembled production residual on a mesh.",
        ),
        TheoryEntry(
            id="sa_fw_r",
            family="turbulence",
            genus="SA",
            section="source / destruction",
            title="fw(r) wall destruction function",
            role_tags=("source", "destruction", "wall calibration"),
            formula=(
                "r = min( rs / S_tilde , 10 )\n"
                "g = r + cw2 * (r^6 - r)\n"
                "fw = g * ((1 + cw3^6) / (g^6 + cw3^6))^(1/6)\n"
                "D ~ cw1 * fw * (nu_tilde / d)^2\n"
                "standard: cw2 = 0.3, cw3 = 2"
            ),
            principle=(
                "Destruction must recover the log layer at r=1 (fw=1) and "
                "strengthen the sink in the outer region (r>1 => fw>1). "
                "Near the wall r->0 => fw->0 with nu_tilde."
            ),
            background=(
                "A bare (nu_tilde/d)^2 sink is not enough to match both the "
                "log layer and outer wake. fw is the calibrated multiplier "
                "in the destruction term of the working-variable equation."
            ),
            references=(
                "Spalart and Allmaras (1992/1994).",
                "NASA TMR — Spalart-Allmaras.",
                "OneFLOW: doc/source/numerical/sa_closure.rst",
            ),
            csv_relpath="tools/verification/data/sa_fw_r_theory.csv",
            plot_module="plot_sa_fw",
            pure_apis=("SaRRatio", "SaG", "SaFw", "SaDestruction"),
            doc_relpath="doc/source/numerical/sa_closure.rst",
            proves="Closed-form fw(r) vs SaG/SaFw sample points.",
            does_not_prove="Skin friction or full-case CFD accuracy.",
        ),
        TheoryEntry(
            id="sa_source_residual",
            family="turbulence",
            genus="SA",
            section="source / residual assembly",
            title="SA source residual assembly",
            role_tags=("source", "production", "destruction", "diffusion"),
            formula=(
                "P = cb1 * (1 - ft2) * S_tilde * nu_tilde\n"
                "D = (cw1k * fw - cb1 * ft2) * nu_tilde_rs\n"
                "D_diff = min( rp * P , (cb2/sigma) * |grad nu|^2 / Re )\n"
                "S = P - D + D_diff"
            ),
            principle=(
                "The working-variable residual combines production, "
                "destruction, and a limited cross-diffusion term. Pure "
                "APIs isolate each algebraic piece for unit tests."
            ),
            background=(
                "Caller supplies precombined coefficients (cw1k, nuTildeRs, "
                "cb2OverSigma) from active constants and local state. See "
                "sa_source.rst for the mapping table to TurbCom::CalcSrcSa."
            ),
            references=(
                "Spalart-Allmaras source term structure.",
                "NASA TMR SA.",
                "OneFLOW: doc/source/numerical/sa_source.rst",
            ),
            csv_relpath=None,
            plot_module=None,
            pure_apis=(
                "SaProduction",
                "SaDestruction",
                "SaDiffusion",
                "SaLimitedDiffusion",
                "SaSourceResidual",
            ),
            doc_relpath="doc/source/numerical/sa_source.rst",
            proves="Algebra of residual pieces against closed-form numbers.",
            does_not_prove="Discrete conservation on a mesh or experimental match.",
        ),
    ]


def verification_root(repository_root: Path) -> Path:
    return (repository_root / "tools" / "verification").resolve()


def resolve_repo_path(repository_root: Path, relpath: str | None) -> Path | None:
    if not relpath:
        return None
    # Normalize separators so Windows/Linux both work from POSIX-style relpaths
    parts = relpath.replace("\\", "/").split("/")
    return repository_root.joinpath(*parts).resolve()


def csv_path(repository_root: Path, entry: TheoryEntry) -> Path | None:
    return resolve_repo_path(repository_root, entry.csv_relpath)


def doc_path(repository_root: Path, entry: TheoryEntry) -> Path | None:
    return resolve_repo_path(repository_root, entry.doc_relpath)


def render_png(repository_root: Path, entry: TheoryEntry, out_png: Path) -> Path:
    import runpy
    import sys

    if not entry.plot_module or not entry.csv_relpath:
        raise ValueError(f"entry {entry.id} has no theory curve plot")
    script = verification_root(repository_root) / f"{entry.plot_module}.py"
    csv = csv_path(repository_root, entry)
    if not script.is_file():
        raise FileNotFoundError(f"missing plot script: {script}")
    if csv is None or not csv.is_file():
        raise FileNotFoundError(f"missing theory CSV: {csv}")
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    old_argv = sys.argv[:]
    try:
        sys.argv = [str(script), "--csv", str(csv), "--out", str(out_png)]
        runpy.run_path(str(script), run_name="__main__")
    finally:
        sys.argv = old_argv
    if not out_png.is_file():
        raise RuntimeError(f"plot script did not write {out_png}")
    return out_png


def taxonomy_path(entry: TheoryEntry) -> tuple[str, str, str]:
    return (entry.family, entry.genus, entry.section)
