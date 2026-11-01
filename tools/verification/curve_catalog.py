"""Contract Workbench knowledge catalog (repository-relative paths)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class TheoryEntry:
    id: str
    family: str
    genus: str
    section: str
    title: str
    role_tags: tuple[str, ...]
    formula_plain: str
    formula_latex: str
    principle: str
    background: str
    references: tuple[str, ...]
    csv_relpath: str | None
    plot_module: str | None
    pure_apis: tuple[str, ...]
    source_relpaths: tuple[str, ...]
    doc_relpath: str | None
    code_sketch: str
    proves: str
    does_not_prove: str


def list_theory_entries() -> list[TheoryEntry]:
    return [
        # ========== SA ==========
        TheoryEntry(
            id="sa_fv1_chi",
            family="turbulence",
            genus="SA",
            section="closure / eddy viscosity",
            title="fv1(chi) eddy-viscosity damping",
            role_tags=("closure", "eddy-viscosity"),
            formula_plain=(
                "chi = nu_tilde / nu\n"
                "fv1 = chi^3 / (chi^3 + cv1^3)\n"
                "mu_t = rho * nu_tilde * fv1\n"
                "cv1 = 7.1"
            ),
            formula_latex=(
                r"\chi=\tilde\nu/\nu"
                "\n"
                r"f_{v1}=\frac{\chi^{3}}{\chi^{3}+c_{v1}^{3}}"
                "\n"
                r"\mu_{t}=\rho\tilde\nu f_{v1}"
            ),
            principle=(
                "SA transports nu_tilde, not mu_t. fv1 damps mu_t near the wall "
                "and approaches 1 in the core."
            ),
            background=(
                "Without fv1, mu_t=rho*nu_tilde over-predicts near-wall stress."
            ),
            references=(
                "Spalart and Allmaras (1992/1994).",
                "NASA TMR SA: https://turbmodels.larc.nasa.gov/spalart.html",
                "OneFLOW: doc/source/numerical/sa_eddy_viscosity.rst",
            ),
            csv_relpath="tools/verification/data/sa_fv1_chi_theory.csv",
            plot_module="plot_sa_fv1",
            pure_apis=("SaChi", "SaFv1", "SaEddyViscosity"),
            source_relpaths=(
                "codes/turb/include/TurbSaPolicy.h",
                "codes/turb/src/TurbSaPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sa_eddy_viscosity.rst",
            code_sketch=(
                "// GATED SKETCH — review before use.\n"
                "Real SaFv1(Real chi, Real cv1Cubed) {\n"
                "  Real chi3 = chi*chi*chi;\n"
                "  return chi3 / (chi3 + cv1Cubed);\n"
                "}\n"
            ),
            proves="Closed-form fv1 vs SaFv1.",
            does_not_prove="Full residual or experimental Cf.",
        ),
        TheoryEntry(
            id="sa_fv2_chi",
            family="turbulence",
            genus="SA",
            section="closure / production strain",
            title="fv2(chi) and modified strain",
            role_tags=("closure", "production", "modified strain"),
            formula_plain=(
                "fv2 = 1 - chi/(1 + chi*fv1)\n"
                "S_bar = rs * fv2\n"
                "S_tilde = Omega + S_bar  (or limited branch)"
            ),
            formula_latex=(
                r"f_{v2}=1-\frac{\chi}{1+\chi f_{v1}}"
                "\n"
                r"\bar{S}=r_s f_{v2}"
                "\n"
                r"\tilde{S}=\Omega+\bar{S}"
            ),
            principle=(
                "Production uses modified strain so near-wall damping stays "
                "consistent with fv1."
            ),
            background=("Limited branch uses c2=0.7, c3=0.9 when needed."),
            references=(
                "NASA TMR SA; OneFLOW sa_closure.rst",
            ),
            csv_relpath="tools/verification/data/sa_fv2_chi_theory.csv",
            plot_module="plot_sa_fv2",
            pure_apis=("SaFv2", "SaSBar", "SaModifiedStrainRate"),
            source_relpaths=(
                "codes/turb/include/TurbSaPolicy.h",
                "codes/turb/src/TurbSaPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sa_closure.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "Real SaFv2(Real chi, Real fv1) {\n"
                "  return 1.0 - chi / (1.0 + chi * fv1);\n"
                "}\n"
            ),
            proves="fv2 / strain vs pure APIs.",
            does_not_prove="Assembled production residual on a mesh.",
        ),
        TheoryEntry(
            id="sa_fw_r",
            family="turbulence",
            genus="SA",
            section="source / destruction",
            title="fw(r) wall destruction function",
            role_tags=("source", "destruction"),
            formula_plain=(
                "r = min(rs/S_tilde, 10)\n"
                "g = r + cw2*(r^6 - r)\n"
                "fw = g*((1+cw3^6)/(g^6+cw3^6))^(1/6)\n"
                "D ~ cw1*fw*(nu_tilde/d)^2"
            ),
            formula_latex=(
                r"r=\min(r_s/\tilde{S},10)"
                "\n"
                r"g=r+c_{w2}(r^{6}-r)"
                "\n"
                r"f_w=g\left(\frac{1+c_{w3}^{6}}{g^{6}+c_{w3}^{6}}\right)^{1/6}"
            ),
            principle=(
                "r=1 => fw=1 (log layer); r>1 strengthens outer destruction; "
                "r->0 => fw->0 near wall."
            ),
            background=("Calibrated multiplier on the (nu/d)^2 sink."),
            references=(
                "Spalart and Allmaras; NASA TMR SA; sa_closure.rst",
            ),
            csv_relpath="tools/verification/data/sa_fw_r_theory.csv",
            plot_module="plot_sa_fw",
            pure_apis=("SaRRatio", "SaG", "SaFw", "SaDestruction"),
            source_relpaths=(
                "codes/turb/include/TurbSaPolicy.h",
                "codes/turb/src/TurbSaPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sa_closure.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "Real SaG(Real r, Real cw2) {\n"
                "  return r + cw2 * (std::pow(r,6) - r);\n"
                "}\n"
            ),
            proves="fw(r) vs SaG/SaFw.",
            does_not_prove="Skin friction or full-case CFD.",
        ),
        TheoryEntry(
            id="sa_source_residual",
            family="turbulence",
            genus="SA",
            section="source / residual assembly",
            title="SA source residual assembly",
            role_tags=("source", "production", "destruction", "diffusion"),
            formula_plain=(
                "P = cb1*(1-ft2)*S_tilde*nu_tilde\n"
                "D = (cw1k*fw - cb1*ft2)*nu_tilde_rs\n"
                "D_diff = min(rp*P, (cb2/sigma)*|grad nu|^2/Re)\n"
                "S = P - D + D_diff"
            ),
            formula_latex=(
                r"P=c_{b1}(1-f_{t2})\tilde{S}\tilde\nu"
                "\n"
                r"D=(c_{w1k}f_w-c_{b1}f_{t2})\tilde\nu r_s"
                "\n"
                r"S=P-D+D_{\mathrm{diff}}"
            ),
            principle=(
                "Working-variable residual = production - destruction + limited diffusion."
            ),
            background=("See sa_source.rst for coefficient mapping in CalcSrcSa."),
            references=("NASA TMR SA; OneFLOW sa_source.rst",),
            csv_relpath=None,
            plot_module=None,
            pure_apis=(
                "SaProduction",
                "SaDestruction",
                "SaDiffusion",
                "SaLimitedDiffusion",
                "SaSourceResidual",
            ),
            source_relpaths=(
                "codes/turb/include/TurbSaPolicy.h",
                "codes/turb/src/TurbCom.cpp",
            ),
            doc_relpath="doc/source/numerical/sa_source.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "Real SaSourceResidual(Real P, Real D, Real Ddiff) {\n"
                "  return P - D + Ddiff;\n"
                "}\n"
            ),
            proves="Residual piece algebra.",
            does_not_prove="Mesh residual or experiment.",
        ),
        # ========== SST ==========
        TheoryEntry(
            id="sst_f1",
            family="turbulence",
            genus="SST",
            section="closure / blending",
            title="F1 blending function",
            role_tags=("closure", "blending"),
            formula_plain=(
                "CD_kw = max(2*rho*sigw2*(gradk.gradw)/w, small)\n"
                "arg1 = min( max(sqrt(k)/(beta_star*w*d), 500*nu/(w*d^2*Re)), "
                "4*rho*sigw2*k/(CD_kw*d^2) )\n"
                "F1 = tanh(arg1^4)\n"
                "phi = F1*phi1 + (1-F1)*phi2"
            ),
            formula_latex=(
                r"F_1=\tanh(\mathrm{arg}_1^{4})"
                "\n"
                r"\phi=F_1\phi_1+(1-F_1)\phi_2"
            ),
            principle=(
                "F1 blends inner (k-omega) and outer (k-epsilon) SST constants."
            ),
            background=("Continuous switch; pure APIs isolate arg1, F1, Blend."),
            references=(
                "Menter SST; NASA TMR SST: https://turbmodels.larc.nasa.gov/sst.html",
                "OneFLOW: doc/source/numerical/sst_f1_blending.rst",
            ),
            csv_relpath="tools/verification/data/sst_f1_arg1_theory.csv",
            plot_module="plot_sst_f1",
            pure_apis=(
                "SstCrossDiffScalar",
                "SstBlendingArg1",
                "SstF1FromArg1",
                "SstBlend",
            ),
            source_relpaths=(
                "codes/turb/include/TurbSstPolicy.h",
                "codes/turb/src/TurbSstPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sst_f1_blending.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "Real SstF1FromArg1(Real arg1) {\n"
                "  return std::tanh(std::pow(arg1, 4));\n"
                "}\n"
                "Real SstBlend(Real F1, Real phi1, Real phi2) {\n"
                "  return F1 * phi1 + (1.0 - F1) * phi2;\n"
                "}\n"
            ),
            proves="F1 / Blend algebra vs pure APIs.",
            does_not_prove="Full SST case accuracy.",
        ),
        TheoryEntry(
            id="sst_f2_eddy",
            family="turbulence",
            genus="SST",
            section="closure / eddy viscosity",
            title="F2(arg2) limiter + eddy viscosity",
            role_tags=("closure", "eddy-viscosity"),
            formula_plain=(
                "arg2 = max(2*sqrt(k)/(beta_star*w*d), 500*nu/(w*d^2*Re))\n"
                "F2 = tanh(arg2^2)\n"
                "mu_t = a1*rho*k / max(a1*w, Omega*F2) * Re"
            ),
            formula_latex=(
                r"F_2=\tanh(\mathrm{arg}_2^{2})"
                "\n"
                r"\mu_t=\frac{a_1\rho k}{\max(a_1\omega,\Omega F_2)}\,Re"
            ),
            principle=(
                "F2 limits eddy viscosity in boundary layers (Bradshaw)."
            ),
            background=("Optional laminar-ratio cap may still apply in the solver."),
            references=(
                "Menter SST; NASA TMR SST; sst_eddy_viscosity.rst",
            ),
            csv_relpath="tools/verification/data/sst_f2_arg2_theory.csv",
            plot_module="plot_sst_f2",
            pure_apis=("SstF2Arg", "SstF2FromArg", "SstEddyViscosity"),
            source_relpaths=(
                "codes/turb/include/TurbSstPolicy.h",
                "codes/turb/src/TurbSrcFlux.cpp",
            ),
            doc_relpath="doc/source/numerical/sst_eddy_viscosity.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "Real SstF2FromArg(Real arg2) {\n"
                "  return std::tanh(arg2 * arg2);\n"
                "}\n"
            ),
            proves="F2(arg2)=tanh(arg2^2) curve vs SstF2FromArg; eddy-viscosity algebra samples.",
            does_not_prove="Wall-resolved SST validation or full mu_t field.",
        ),
        TheoryEntry(
            id="sst_k_source",
            family="turbulence",
            genus="SST",
            section="source / k-equation",
            title="SST k production and dissipation",
            role_tags=("source", "production", "dissipation"),
            formula_plain=(
                "P_k from Boussinesq / strain / vorticity options\n"
                "eps = f_beta_star * beta_star * rho * k * omega\n"
                "P_k <- min(P_k, pklim * eps)"
            ),
            formula_latex=(
                r"\varepsilon=f_{\beta^*}\beta^*\rho k\omega"
                "\n"
                r"P_k\leftarrow\min(P_k,\,p_{k\lim}\varepsilon)"
            ),
            principle=("k residual ~ production minus dissipation with limiter."),
            background=("Pure APIs isolate production options and limiter."),
            references=("Menter SST; OneFLOW sst_k_source.rst",),
            csv_relpath=None,
            plot_module=None,
            pure_apis=("SstProdkBoussinesq", "SstDisskStandard", "SstLimitProdk"),
            source_relpaths=(
                "codes/turb/include/TurbSstPolicy.h",
                "codes/turb/src/TurbSstPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sst_k_source.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "Real SstLimitProdk(Real P, Real pklim, Real eps) {\n"
                "  return std::min(P, pklim * eps);\n"
                "}\n"
            ),
            proves="k source algebra vs pure APIs.",
            does_not_prove="Two-equation field quality.",
        ),
        TheoryEntry(
            id="sst_omega_source",
            family="turbulence",
            genus="SST",
            section="source / omega-equation",
            title="SST omega production, dissipation, cross-diffusion",
            role_tags=("source", "production", "dissipation", "cross-diffusion"),
            formula_plain=(
                "P_w ~ alpha * rho * S^2 (or related forms)\n"
                "D_w = beta * rho * omega^2\n"
                "CD ~ 2*(1-F1)*rho*sigw2*(gradk·gradw)/omega"
            ),
            formula_latex=(
                r"D_\omega=\beta\rho\omega^{2}"
                "\n"
                r"CD_{\omega}\sim 2(1-F_1)\rho\sigma_{\omega 2}"
                r"\frac{\nabla k\cdot\nabla\omega}{\omega}"
            ),
            principle=(
                "Omega equation carries production, dissipation, and cross-diffusion "
                "active mainly in the outer blend region."
            ),
            background=("See sst_omega_source.rst and TurbSstPolicy pure APIs."),
            references=("Menter SST; NASA TMR SST; sst_omega_source.rst",),
            csv_relpath=None,
            plot_module=None,
            pure_apis=(
                "SstProdwMenter",
                "SstDissw",
                "SstCrossDiffTerm",
            ),
            source_relpaths=(
                "codes/turb/include/TurbSstPolicy.h",
                "codes/turb/src/TurbSstPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sst_omega_source.rst",
            code_sketch=(
                "// GATED SKETCH — see TurbSstPolicy for exact signatures.\n"
            ),
            proves="Omega source piece APIs where extracted.",
            does_not_prove="Full SST residual on a case.",
        ),

        # ========== Wilcox ==========
        TheoryEntry(
            id="wilcox_fk_xk",
            family="turbulence",
            genus="Wilcox",
            section="closure / beta-star modifier",
            title="fk(xk) Wilcox beta* modifier",
            role_tags=("closure", "Wilcox", "beta-star"),
            formula_plain=(
                "xk = (grad k · grad omega) / omega^3\n"
                "fk = 1  if xk <= 0\n"
                "fk = (1+680 xk^2)/(1+400 xk^2)  if xk > 0"
            ),
            formula_latex=(
                r"x_k=(\nabla k\cdot\nabla\omega)/\omega^{3}"
                "\n"
                r"f_k=1\ (x_k\le 0);\quad "
                r"(1+680 x_k^{2})/(1+400 x_k^{2})\ (x_k>0)"
            ),
            principle=(
                "Wilcox beta-star is modulated by fk(xk) from the k-omega "
                "cross-diffusion sense."
            ),
            background=(
                "Used on Wilcox/EASM paths in TurbCom. SST uses F1; this is the Wilcox piece."
            ),
            references=(
                "Wilcox turbulence modeling notes.",
                "OneFLOW: doc/source/numerical/wilcox_fbeta.rst",
            ),
            csv_relpath="tools/verification/data/wilcox_fk_xk_theory.csv",
            plot_module="plot_wilcox_fk",
            pure_apis=("WilcoxXkFromCrossTerm", "WilcoxFkFromXk"),
            source_relpaths=(
                "codes/turb/include/TurbSstPolicy.h",
                "codes/turb/src/TurbSstPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/wilcox_fbeta.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "Real WilcoxFkFromXk(Real xk) {\n"
                "  if (xk <= 0) return 1;\n"
                "  Real xk2 = xk*xk;\n"
                "  return (1 + 680*xk2) / (1 + 400*xk2);\n"
                "}\n"
            ),
            proves="fk(xk) closed form vs WilcoxFkFromXk sample grid.",
            does_not_prove="Full k-omega residual or SST F1 path.",
        ),

        # ========== NS / LUSGS ==========
        TheoryEntry(
            id="ns_lusgs_viscous_spectrum",
            family="mean-flow",
            genus="NS-LUSGS",
            section="implicit / viscous spectrum",
            title="LUSGS viscous spectral residual contribution",
            role_tags=("implicit", "viscous", "spectrum"),
            formula_plain=(
                "rhs <- rhs - visrad * dq\n"
                "model1: visrad = vissr / V_cell\n"
                "face density = 0.5*(rho_L + rho_R)"
            ),
            formula_latex=(
                r"r\leftarrow r-\lambda_{\mathrm{vis}}\Delta q"
                "\n"
                r"\lambda_{\mathrm{vis}}=\lambda_{\mathrm{face}}/V"
            ),
            principle=(
                "After inviscid Jacobian increment, viscous spectral radius damps neighbor dq."
            ),
            background=(
                "Uses ViscousFaceSpectralRadiusModel1/2; pure contribution is sign-clear."
            ),
            references=(
                "Yoon-Jameson LUSGS; OneFLOW lusgs_viscous_spectrum.rst",
            ),
            csv_relpath=None,
            plot_module=None,
            pure_apis=(
                "LusgsViscousSpectrumContribution",
                "LusgsFaceArithmeticMean",
                "LusgsViscousRadiusPerVolume",
            ),
            source_relpaths=(
                "codes/ns/include/NsLusgsPolicy.h",
                "codes/uns/src/UNsLusgs.cpp",
            ),
            doc_relpath="doc/source/numerical/lusgs_viscous_spectrum.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "Real LusgsViscousSpectrumContribution(Real visrad, Real dq) {\n"
                "  return -visrad * dq;\n"
                "}\n"
            ),
            proves="Contribution algebra vs pure API.",
            does_not_prove="Full LUSGS case convergence.",
        ),
        TheoryEntry(
            id="ns_lusgs_flux_increment",
            family="mean-flow",
            genus="NS-LUSGS",
            section="implicit / inviscid Jacobian",
            title="LUSGS five-equation flux increment",
            role_tags=("implicit", "inviscid", "jacobian"),
            formula_plain=(
                "Split face eigenvalues with sign of sweep matrix\n"
                "Linearization weights x1, x2 from lambda and sound speed\n"
                "df = approximate A * dq * faceArea"
            ),
            formula_latex=(
                r"\lambda_i=\frac12(\lambda_i^{\mathrm{phys}}"
                r"+s|\lambda|_{\mathrm{face},i})"
                "\n"
                r"\Delta\mathbf{F}\approx A\Delta\mathbf{q}\,S_f"
            ),
            principle=(
                "Approximate inviscid flux Jacobian times neighbor dq for LUSGS sweeps."
            ),
            background=("Pure APIs: LusgsSplitFaceEigenvalues, LinearizationWeights, FluxIncrementFiveEq."),
            references=("Yoon-Jameson LUSGS; NsLusgsPolicy",),
            csv_relpath=None,
            plot_module=None,
            pure_apis=(
                "LusgsSplitFaceEigenvalues",
                "LusgsLinearizationWeights",
                "LusgsFluxIncrementFiveEq",
            ),
            source_relpaths=(
                "codes/ns/include/NsLusgsPolicy.h",
                "codes/ns/src/NsLusgs.cpp",
            ),
            doc_relpath="doc/source/numerical/index.rst",
            code_sketch=(
                "// GATED SKETCH — see NsLusgsPolicy.h signatures.\n"
            ),
            proves="Five-eq increment unit tests where present.",
            does_not_prove="Sweep stability on arbitrary meshes.",
        ),
        TheoryEntry(
            id="ns_roe_harten",
            family="mean-flow",
            genus="NS-inviscid",
            section="flux / Roe dissipation",
            title="Harten entropy fix (eigenvalue floor)",
            role_tags=("flux", "entropy-fix"),
            formula_plain=(
                "If enabled, floor |lambda| using max-eigenvalue * centropy factors\n"
                "(production uses MAX macro semantics)"
            ),
            formula_latex=(
                r"|\lambda|\leftarrow\max(|\lambda|,m)"
                r"\quad\text{(Harten-type floor)}"
            ),
            principle=(
                "Floors small Roe wave speeds to reduce expansion shocks."
            ),
            background=("Pure policy ApplyHartenEntropyFix; tests must match MAX semantics."),
            references=("Harten entropy fix; NsInvFluxPolicy tests",),
            csv_relpath=None,
            plot_module=None,
            pure_apis=("ApplyHartenEntropyFix",),
            source_relpaths=(
                "codes/ns/include/NsInvFluxPolicy.h",
                "codes/ns/src/NsInvFluxPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/index.rst",
            code_sketch=(
                "// GATED SKETCH — mirror production MAX in tests.\n"
            ),
            proves="Entropy-fix pure policy tests.",
            does_not_prove="Full Roe flux shock-capturing quality.",
        ),
    ]


def verification_root(repository_root: Path) -> Path:
    return (repository_root / "tools" / "verification").resolve()


def resolve_repo_path(repository_root: Path, relpath: str | None) -> Path | None:
    if not relpath:
        return None
    parts = relpath.replace("\\", "/").split("/")
    return repository_root.joinpath(*parts).resolve()


def csv_path(repository_root: Path, entry: TheoryEntry) -> Path | None:
    return resolve_repo_path(repository_root, entry.csv_relpath)


def doc_path(repository_root: Path, entry: TheoryEntry) -> Path | None:
    return resolve_repo_path(repository_root, entry.doc_relpath)


def source_paths(repository_root: Path, entry: TheoryEntry) -> list[Path]:
    return [
        p
        for rel in entry.source_relpaths
        if (p := resolve_repo_path(repository_root, rel)) is not None
    ]


def read_doc_excerpt(repository_root: Path, entry: TheoryEntry, max_chars: int = 4000) -> str:
    path = doc_path(repository_root, entry)
    if path is None or not path.is_file():
        return "(Sphinx page not found under repository_root.)"
    text = path.read_text(encoding="utf-8", errors="replace")
    if len(text) > max_chars:
        return text[:max_chars] + "\n\n... [truncated] ..."
    return text


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
    old = sys.argv[:]
    try:
        sys.argv = [str(script), "--csv", str(csv), "--out", str(out_png)]
        runpy.run_path(str(script), run_name="__main__")
    finally:
        sys.argv = old
    if not out_png.is_file():
        raise RuntimeError(f"plot script did not write {out_png}")
    return out_png
