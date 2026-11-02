"""Contract Workbench knowledge catalog (repository-relative paths)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re


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
                "CD_kw = max(2*rho*sigw2*(gradk·gradw)/w, small)\n"
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
                "mu_t = rho*k / max(w, R*F2/a1) * Re (R = selected strain/vorticity magnitude)"
            ),
            formula_latex=(
                r"F_2=\tanh(\mathrm{arg}_2^{2})"
                "\n"
                r"\mu_t=\frac{\rho k}{\max(\omega,R F_2/a_1)}\,Re"
            ),
            principle=(
                "F2 limits eddy viscosity using the selected strain-rate or vorticity magnitude (Bradshaw limiter)."
            ),
            background=("Optional laminar-ratio cap may still apply in the solver."),
            references=(
                "Menter SST; NASA TMR SST; sst_eddy_viscosity.rst",
            ),
            csv_relpath="tools/verification/data/sst_f2_arg2_theory.csv",
            plot_module="plot_sst_f2",
            pure_apis=(
                "SstF2Argument", "SstF2ArgFromTerms",
                "SstF2FromArg", "SstEddyViscosity",
                "SstComputeEddyViscosity",
            ),
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
            proves=(
                "Physical-input arg2 construction selects the larger length-scale or viscous term; "
                "F2(arg2)=tanh(arg2^2) curve, complete per-cell closure composition, "
                "production-option rate selection, and shared laminar-ratio cap."
            ),
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
                "CD ~ 2*(1-F1)*rho*sigw2*(gradk.gradw)/omega"
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
                "SstCrossDiffMenter",
            ),
            source_relpaths=(
                "codes/turb/include/TurbSstPolicy.h",
                "codes/turb/src/TurbSstPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sst_omega_source.rst",
            code_sketch=(
                "// GATED SKETCH - see TurbSstPolicy for exact signatures.\n"
            ),
            proves="Omega source piece APIs where extracted.",
            does_not_prove="Full SST residual on a case.",
        ),


        TheoryEntry(
            id="sst_residual_diagonal",
            family="turbulence",
            genus="SST",
            section="source / residual linearization",
            title="SST source residual and diagonal pieces",
            role_tags=("source", "linearization", "SST"),
            formula_plain=(
                "S_k = P_k - eps\n"
                "S_w = P_w - D_w + CD\n"
                "d_k = -2 * f_beta* * beta* * omega\n"
                "d_w = -2 * f_beta * beta * omega - |CD|/(rho*omega)\n"
                "neg = 0.5 * (S - |S|)"
            ),
            formula_latex=(
                r"S_k=P_k-\varepsilon"
                "\n"
                r"S_\omega=P_\omega-D_\omega+CD"
                "\n"
                r"d_k=-2 f_{\beta^*}\beta^*\omega"
                "\n"
                r"d_\omega=-2 f_\beta\beta\omega-|CD|/(\rho\omega)"
            ),
            principle=(
                "Menter linearization splits explicit source and diagonal "
                "contributions for implicit turbulence updates."
            ),
            background=(
                "Pure APIs SstSourceK/W, SstDiagonalK/W, SstNegativeSourcePart "
                "are locked in turb_sst_residual_test."
            ),
            references=(
                "Menter SST linearization practice.",
                "OneFLOW TurbSstPolicy residual section.",
            ),
            csv_relpath=None,
            plot_module=None,
            pure_apis=(
                "SstSourceK",
                "SstSourceW",
                "SstDiagonalK",
                "SstDiagonalW",
                "SstNegativeSourcePart",
            ),
            source_relpaths=(
                "codes/turb/include/TurbSstPolicy.h",
                "codes/turb/src/TurbSstPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sst_k_source.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "Real SstDiagonalW(Real fBeta, Real beta, Real omega,\n"
                "                  Real absCD, Real rho) {\n"
                "  return -2*fBeta*beta*omega - absCD/(rho*omega);\n"
                "}\n"
            ),
            proves="Residual piece algebra vs pure APIs.",
            does_not_prove="Full SST field residual or case Cf.",
        ),


        TheoryEntry(
            id="sst_wall_omega_menter",
            family="turbulence",
            genus="SST",
            section="boundary / wall omega",
            title="Wall omega_w (Menter)",
            role_tags=("boundary", "wall", "omega"),
            formula_plain=(
                "omega_w = 60 * mu_l / (rho * beta1 * d^2 * Re)\n"
                "with beta1 ~ 0.075 (inner SST constant)"
            ),
            formula_latex=(
                r"\omega_w=\frac{60\mu_l}{\rho\beta_1 d^{2} Re}"
            ),
            principle=(
                "Rough-wall / smooth-wall SST practice sets a large omega at "
                "the wall scaling as 1/d^2 so the viscous sublayer is resolved "
                "without integrating to omega->inf."
            ),
            background=(
                "Pure API TurbWallOmegaMenter; image conditions for k and SA "
                "nu_tilde use sign-flip ghosts (TurbWallKineticEnergyImage, "
                "TurbWallSaImage)."
            ),
            references=(
                "Menter SST wall boundary practice; NASA TMR SST.",
                "OneFLOW TurbBcPolicy.",
            ),
            csv_relpath="tools/verification/data/sst_wall_omega_menter_theory.csv",
            plot_module="plot_sst_wall_omega",
            pure_apis=("TurbWallOmegaMenter", "TurbWallKineticEnergyImage"),
            source_relpaths=(
                "codes/turb/include/TurbBcPolicy.h",
                "codes/turb/src/TurbBcPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/index.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "Real TurbWallOmegaMenter(Real mu, Real rho, Real beta1,\n"
                "                         Real d2, Real Re) {\n"
                "  return 60*mu / (rho*beta1*d2*Re);\n"
                "}\n"
            ),
            proves="omega_w(d) closed form vs TurbWallOmegaMenter grid.",
            does_not_prove="y+ resolution or full SST wall shear.",
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

        TheoryEntry(
            id="wilcox_fw1998",
            family="turbulence",
            genus="Wilcox",
            section="closure / beta modifier",
            title="fw1998(xw) Wilcox beta modifier",
            role_tags=("closure", "Wilcox", "beta"),
            formula_plain="fw1998 = (1 + 70*xw) / (1 + 80*xw)",
            formula_latex=r"f_w^{1998}=(1+70 x_w)/(1+80 x_w)",
            principle=(
                "1998 Wilcox fw softens beta based on vortex-stretching parameter xw."
            ),
            background=(
                "As xw->inf, fw->7/8. Used with Wilcox beta/beta* paths in TurbCom."
            ),
            references=(
                "Wilcox (1998) turbulence model notes.",
                "OneFLOW: doc/source/numerical/wilcox_fbeta.rst",
            ),
            csv_relpath="tools/verification/data/wilcox_fw1998_theory.csv",
            plot_module="plot_wilcox_fw",
            pure_apis=("WilcoxFw1998",),
            source_relpaths=(
                "codes/turb/include/TurbSstPolicy.h",
                "codes/turb/src/TurbSstPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/wilcox_fbeta.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "Real WilcoxFw1998(Real xw) {\n"
                "  return (1 + 70*xw) / (1 + 80*xw);\n"
                "}\n"
            ),
            proves="fw1998 closed form vs WilcoxFw1998 sample grid.",
            does_not_prove="Full Wilcox residual or SST path.",
        ),
        TheoryEntry(
            id="wilcox_fw2006",
            family="turbulence",
            genus="Wilcox",
            section="closure / beta modifier",
            title="fw2006(xw) Wilcox beta modifier",
            role_tags=("closure", "Wilcox", "beta"),
            formula_plain="fw2006 = (1 + 85*xw) / (1 + 100*xw)",
            formula_latex=r"f_w^{2006}=(1+85 x_w)/(1+100 x_w)",
            principle=(
                "2006 Wilcox fw variant; limit as xw->inf is 85/100=0.85."
            ),
            background=(
                "Companion to fw1998; selection depends on model variant flags."
            ),
            references=(
                "Wilcox (2006) notes; OneFLOW wilcox_fbeta.rst",
            ),
            csv_relpath="tools/verification/data/wilcox_fw2006_theory.csv",
            plot_module="plot_wilcox_fw",
            pure_apis=("WilcoxFw2006",),
            source_relpaths=(
                "codes/turb/include/TurbSstPolicy.h",
                "codes/turb/src/TurbSstPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/wilcox_fbeta.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "Real WilcoxFw2006(Real xw) {\n"
                "  return (1 + 85*xw) / (1 + 100*xw);\n"
                "}\n"
            ),
            proves="fw2006 closed form vs WilcoxFw2006 sample grid.",
            does_not_prove="Full Wilcox residual or SST path.",
        ),



        TheoryEntry(
            id="turb_inviscid_spectrum",
            family="turbulence",
            genus="Turb-spectrum",
            section="implicit / face spectral radius",
            title="Inviscid face spectral radius",
            role_tags=("implicit", "spectrum", "turbulence"),
            formula_plain=(
                "lambda_inv = (1/2) * A * |vn_mid|\n"
                "matrix L/R += (A/2) * (-/+ vn - |vn|) contributions"
            ),
            formula_latex=(
                r"\lambda_{\mathrm{inv}}=\frac{A}{2}|v_n|"
            ),
            principle=(
                "Turbulence LUSGS uses a scalar face spectral radius from "
                "mid-face normal speed, analogous to mean-flow inviscid radius."
            ),
            background=(
                "Pure APIs TurbInviscidFaceSpectralRadius and "
                "TurbInviscidFaceMatrixContrib; SA also uses "
                "SaEffectiveKinematicViscosity and TurbDiffusiveFaceSpectralRadius."
            ),
            references=(
                "OneFLOW TurbSpectrumPolicy; LUSGS spectral-radius practice.",
            ),
            csv_relpath="tools/verification/data/turb_inviscid_spectrum_theory.csv",
            plot_module="plot_turb_inviscid_spectrum",
            pure_apis=(
                "TurbInviscidFaceSpectralRadius",
                "TurbInviscidFaceMatrixContrib",
                "SaEffectiveKinematicViscosity",
                "TurbDiffusiveFaceSpectralRadius",
            ),
            source_relpaths=(
                "codes/turb/include/TurbSpectrumPolicy.h",
                "codes/turb/src/TurbSpectrumPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/turb_spectrum.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "Real TurbInviscidFaceSpectralRadius(Real A, Real absVn) {\n"
                "  return 0.5 * A * absVn;\n"
                "}\n"
            ),
            proves="lambda_inv algebra vs pure API on |vn| grid.",
            does_not_prove="Full turb LUSGS convergence.",
        ),

        TheoryEntry(
            id="turb_diffusive_spectrum",
            family="turbulence",
            genus="Turb-spectrum",
            section="implicit / face spectral radius",
            title="Diffusive face spectral radius (SA)",
            role_tags=("implicit", "spectrum", "diffusive", "SA"),
            formula_plain=(
                "nu_et = mu_l/rho + nu_tilde\n"
                "lambda_diff = (1/Re)*(1/sigma)*nu_et*A^2/V"
            ),
            formula_latex=(
                r"\nu_{\mathrm{et}}=\mu_l/\rho+\tilde\nu"
                "\n"
                r"\lambda_{\mathrm{diff}}=\frac{1}{Re}\frac{1}{\sigma}"
                r"\nu_{\mathrm{et}}\frac{A^{2}}{V}"
            ),
            principle=(
                "SA turbulence LUSGS adds a diffusive spectral radius from "
                "effective kinematic viscosity and face geometry."
            ),
            background=(
                "Companion to turb_inviscid_spectrum. Pure APIs "
                "SaEffectiveKinematicViscosity and TurbDiffusiveFaceSpectralRadius."
            ),
            references=(
                "OneFLOW TurbSpectrumPolicy; turb_spectrum.rst",
            ),
            csv_relpath="tools/verification/data/turb_diffusive_spectrum_theory.csv",
            plot_module="plot_turb_diffusive_spectrum",
            pure_apis=(
                "SaEffectiveKinematicViscosity",
                "TurbDiffusiveFaceSpectralRadius",
            ),
            source_relpaths=(
                "codes/turb/include/TurbSpectrumPolicy.h",
                "codes/turb/src/TurbSpectrumPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/turb_spectrum.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "Real TurbDiffusiveFaceSpectralRadius(\n"
                "    Real invRe, Real invSigma, Real nuEff,\n"
                "    Real A2, Real invV) {\n"
                "  return invRe * invV * (invSigma * nuEff * A2);\n"
                "}\n"
            ),
            proves="lambda_diff algebra vs pure API on nu_eff grid.",
            does_not_prove="Full turb LUSGS case convergence.",
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
            csv_relpath="tools/verification/data/lusgs_viscous_spectrum_contrib_theory.csv",
            plot_module="plot_lusgs_viscous_spectrum",
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
            id="lusgs_split_eigenvalues",
            family="mean-flow",
            genus="NS-LUSGS",
            section="implicit / inviscid Jacobian",
            title="Face eigenvalue split (matrix sign)",
            role_tags=("implicit", "LUSGS", "eigenvalue"),
            formula_plain=(
                "lmd1 = 1/2 * (vn + s * faceAbs1)\n"
                "lmd2 = 1/2 * (vn + c + s * faceAbs2)\n"
                "lmd3 = 1/2 * (vn - c + s * faceAbs3)"
            ),
            formula_latex=(
                r"\lambda_i=\frac12(\lambda_i^{\mathrm{phys}}+s|\lambda|_{f,i})"
            ),
            principle=(
                "Sweep direction enters through matrix sign s, blending physical "
                "eigenvalues with face absolute spectral radii."
            ),
            background=(
                "Pure API LusgsSplitFaceEigenvalues; weights follow "
                "LusgsLinearizationWeights; flux via LusgsFluxIncrementFiveEq."
            ),
            references=(
                "Yoon-Jameson LUSGS; OneFLOW NsLusgsPolicy",
            ),
            csv_relpath="tools/verification/data/lusgs_split_eigenvalues_theory.csv",
            plot_module="plot_lusgs_split_eigenvalues",
            pure_apis=(
                "LusgsSplitFaceEigenvalues",
                "LusgsLinearizationWeights",
                "LusgsFluxIncrementFiveEq",
            ),
            source_relpaths=(
                "codes/ns/include/NsLusgsPolicy.h",
                "codes/ns/src/NsLusgsPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/lusgs_split_eigenvalues.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "lmd1 = 0.5 * (vn + s * faceAbs1);\n"
                "lmd2 = 0.5 * (vn + c + s * faceAbs2);\n"
                "lmd3 = 0.5 * (vn - c + s * faceAbs3);\n"
            ),
            proves="Split eigenvalues vs closed form for s in {-1,0,+1}.",
            does_not_prove="Full LUSGS sweep stability on a case.",
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
                "// GATED SKETCH - see NsLusgsPolicy.h signatures.\n"
            ),
            proves="Five-eq increment unit tests where present.",
            does_not_prove="Sweep stability on arbitrary meshes.",
        ),

        TheoryEntry(
            id="ns_inviscid_spectrum",
            family="mean-flow",
            genus="NS-spectrum",
            section="time-step / face spectral radius",
            title="Inviscid face spectral radius and cell dt",
            role_tags=("time-step", "spectrum", "mean-flow"),
            formula_plain=(
                "lambda_inv = (1/2) * A * (|vn| + c)\n"
                "dt = CFL * Vol / lambda\n"
                "dt <- dt_inv * dt_vis / (dt_inv + dt_vis)"
            ),
            formula_latex=(
                r"\lambda_{\mathrm{inv}}=\frac{A}{2}(|v_n|+c)"
                "\n"
                r"\Delta t=\mathrm{CFL}\,V/\lambda"
            ),
            principle=(
                "Local time step from inviscid (and viscous) face spectral "
                "radii; inviscid and viscous dt combined harmonically."
            ),
            background=(
                "Pure APIs in TimeStepPolicy; separate from turbulence "
                "TurbSpectrumPolicy which uses |vn| only for scalar turb."
            ),
            references=(
                "Blazek CFD; OneFLOW TimeStepPolicy",
            ),
            csv_relpath="tools/verification/data/ns_inviscid_spectrum_theory.csv",
            plot_module="plot_ns_inviscid_spectrum",
            pure_apis=(
                "InviscidFaceSpectralRadius",
                "CellTimeStepFromSpectrum",
                "CombineInvAndVisTimeStep",
            ),
            source_relpaths=(
                "codes/ns/include/TimeStepPolicy.h",
                "codes/ns/src/TimeStepPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/ns_timestep_spectrum.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "Real InviscidFaceSpectralRadius(Real vn, Real c, Real A) {\n"
                "  return 0.5 * A * (fabs(vn) + c);\n"
                "}\n"
            ),
            proves="lambda_inv and dt algebra vs pure APIs.",
            does_not_prove="Global CFL stability on a case.",
        ),


        TheoryEntry(
            id="ausm_mach_split",
            family="mean-flow",
            genus="NS-inviscid",
            section="flux / AUSM splits",
            title="AUSM Mach and pressure split (FMSplit4 / FPSplit5)",
            role_tags=("flux", "AUSM", "Mach-split"),
                        formula_plain=(
                "M1(M,s) = 0.5 * (M + s*|M|)\n"
                "M2(M,s) = s * 0.25 * (M + s)^2\n"
                "FMSplit4(M,beta,s):\n"
                "  |M| >= 1  ->  M1(M,s)\n"
                "  |M| <  1  ->  M2(M,s) * (1 - s*16*beta*M2(M,-s))\n"
                "FPSplit5(M,alpha,s):\n"
                "  |M| >= 1  ->  0.5 * (1 + s*M/|M|)\n"
                "  |M| <  1  ->  M2(M,s) * ((2*s - M) - s*16*alpha*M*M2(M,-s))\n"
                "s = ipn = +/-1;  beta = 1/8;  alpha = 3/16 (typical AUSM)"
            ),
                        formula_latex=(
                r"M_1=\frac{1}{2}(M+s|M|),\quad M_2=s\frac{1}{4}(M+s)^{2}"
                "\n"
                r"M^{\pm}_{4}=M_1\ (|M|>=1);\ "
                r"M_2\bigl(1-s\,16\beta M_2(M,-s)\bigr)\ (|M|<1)"
                "\n"
                r"P^{\pm}_{5}=\frac{1}{2}(1+s M/|M|)\ (|M|>=1);\ "
                r"M_2\bigl((2s-M)-s\,16\alpha M\,M_2(M,-s)\bigr)\ (|M|<1)"
                "\n"
                r"s=i_{pn}=\pm 1,\ \beta=\frac{1}{8},\ \alpha=\frac{3}{16}"
            ),
principle=(
                "AUSM splits mass and pressure contributions using smooth "
                "polynomials in the subsonic range and upwind steps outside "
                "the sonic interval."
            ),
            background=(
                "Pure inline APIs in NsInvFluxPolicy.h (single definition, "
                "ODR-safe). Locked by inv_flux_policy_test Mach-grid checks."
            ),
            references=(
                "Liou AUSM family; OneFLOW NsInvFluxPolicy",
                "doc/source/numerical/ausm_mach_split.rst",
            ),
            csv_relpath="tools/verification/data/ausm_mach_split_theory.csv",
            plot_module="plot_ausm_mach_split",
            pure_apis=("FMSplit1", "FMSplit2", "FMSplit4", "FPSplit5"),
            source_relpaths=(
                "codes/ns/include/NsInvFluxPolicy.h",
            ),
            doc_relpath="doc/source/numerical/ausm_mach_split.rst",
            code_sketch=(
                "// GATED SKETCH — production is inline in NsInvFluxPolicy.h\n"
                "Real FMSplit1(M,s){ return 0.5*(M + s*fabs(M)); }\n"
                "Real FMSplit2(M,s){ return s*0.25*(M+s)*(M+s); }\n"
                "Real FMSplit4(M,beta,s){\n"
                "  if (fabs(M)>=1) return FMSplit1(M,s);\n"
                "  return FMSplit2(M,s)*(1 - s*16*beta*FMSplit2(M,-s));\n"
                "}\n"
                "Real FPSplit5(M,alpha,s){\n"
                "  if (fabs(M)>=1) return 0.5*(1 + s*M/(fabs(M)+eps));\n"
                "  return FMSplit2(M,s)*((2*s-M) - s*16*alpha*M*FMSplit2(M,-s));\n"
                "}\n"
            ),
            proves="FMSplit4/FPSplit5 closed form vs pure APIs on a Mach grid.",
            does_not_prove="Full AUSM face flux or shock-capturing quality.",
        ),



        TheoryEntry(
            id="roe_dissipation_coupling",
            family="mean-flow",
            genus="NS-inviscid",
            section="flux / Roe dissipation",
            title="Roe dissipation coupling terms (term1, term2)",
            role_tags=("flux", "Roe", "dissipation"),
            formula_plain=(
                "term1 = dh * xi1 + dc * xi2\n"
                "term2 = c^2 * dc * xi1 + dh * xi2\n"
                "dh = total-enthalpy jump; dc = density jump (Roe-scaled)"
            ),
            formula_latex=(
                r"\mathrm{term}_1=\Delta H\,\xi_1+d_c\,\xi_2"
                "\n"
                r"\mathrm{term}_2=c^{2}\,d_c\,\xi_1+\Delta H\,\xi_2"
            ),
            principle=(
                "Maps wave coefficients and jumps into the two scalars that "
                "enter the five-equation Roe dissipation update."
            ),
            background=(
                "Pure API RoeDissipationCouplingTerms. Upstream: xi1/xi2 from "
                "RoeDissipationWaveCoefficients. Downstream: ApplyRoeDissipationFiveEq."
            ),
            references=(
                "Roe (1981); OneFLOW NsInvFluxPolicy",
                "doc/source/numerical/roe_dissipation_coupling.rst",
            ),
            csv_relpath="tools/verification/data/roe_dissipation_coupling_theory.csv",
            plot_module="plot_roe_dissipation_coupling",
            pure_apis=("RoeDissipationCouplingTerms",),
            source_relpaths=(
                "codes/ns/include/NsInvFluxPolicy.h",
                "codes/ns/src/NsInvFluxPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/roe_dissipation_coupling.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "term1 = dh * xi1 + dc * xi2;\n"
                "term2 = c2 * dc * xi1 + dh * xi2;\n"
            ),
            proves="term1/term2 algebra vs pure API on a Delta-H grid.",
            does_not_prove="Full five-eq Roe flux residual on a case mesh.",
        ),

        TheoryEntry(
            id="roe_dissipation_wave_coeff",
            family="mean-flow",
            genus="NS-inviscid",
            section="flux / Roe dissipation",
            title="Roe dissipation wave coefficients xi1, xi2",
            role_tags=("flux", "Roe", "dissipation"),
            formula_plain=(
                "a1=|vn|, a2=|vn+c|, a3=|vn-c| (after absolute + optional Harten)\n"
                "xi1 = (2*a1 - a2 - a3) / (2*c^2)\n"
                "xi2 = (a2 - a3) / (2*c)\n"
                "Used to form Roe dissipation coupling terms."
            ),
            formula_latex=(
                r"a_1=|v_n|,\ a_2=|v_n+c|,\ a_3=|v_n-c|"
                "\n"
                r"\xi_1=\frac{2a_1-a_2-a_3}{2c^{2}},\quad "
                r"\xi_2=\frac{a_2-a_3}{2c}"
            ),
            principle=(
                "Projects absolute wave speeds into the two scalar coefficients "
                "that weight continuity and enthalpy jumps in Roe dissipation."
            ),
            background=(
                "Pure API RoeDissipationWaveCoefficients. Downstream: "
                "RoeDissipationCouplingTerms and ApplyRoeDissipationFiveEq."
            ),
            references=(
                "Roe (1981); OneFLOW NsInvFluxPolicy",
                "doc/source/numerical/roe_dissipation_wave_coeff.rst",
            ),
            csv_relpath="tools/verification/data/roe_dissipation_wave_coeff_theory.csv",
            plot_module="plot_roe_dissipation_wave_coeff",
            pure_apis=("RoeDissipationWaveCoefficients",),
            source_relpaths=(
                "codes/ns/include/NsInvFluxPolicy.h",
                "codes/ns/src/NsInvFluxPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/roe_dissipation_wave_coeff.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "xi1 = (2*a1 - a2 - a3) / (2*c*c);\n"
                "xi2 = (a2 - a3) / (2*c);\n"
            ),
            proves="xi1/xi2 algebra vs pure API on an |vn| grid.",
            does_not_prove="Full Roe flux residual or shock-capturing quality.",
        ),

        TheoryEntry(
            id="ns_roe_harten",
            family="mean-flow",
            genus="NS-inviscid",
            section="flux / Roe dissipation",
            title="Harten entropy fix (eigenvalue floor)",
            role_tags=("flux", "entropy-fix"),
                        formula_plain=(
                "If enabled (else identity):\n"
                "  lambda_max = max(|vn+c|, |vn-c|)\n"
                "  m1 = lambda_max * centropy1\n"
                "  m2 = lambda_max * centropy2\n"
                "  |vn|     <- max(m1, |vn|)\n"
                "  |vn+c|   <- max(m2, |vn+c|)\n"
                "  |vn-c|   <- max(m2, |vn-c|)\n"
                "Uses production MAX macro semantics."
            ),
            formula_latex=(
                r"\lambda_{\max}=\max(|v_n+c|,|v_n-c|)"
                "\n"
                r"m_1=\lambda_{\max} c_1,\quad m_2=\lambda_{\max} c_2"
                "\n"
                r"|v_n|\leftarrow\max(m_1,|v_n|),\ "
                r"|v_n\pm c|\leftarrow\max(m_2,|v_n\pm c|)"
            ),
principle=(
                "Floors small Roe wave speeds to reduce expansion shocks."
            ),
            background=(
                "Pure policy ApplyHartenEntropyFix. Floors are independent per "
                "wave family after absolute eigenvalues are formed."
            ),
            references=(
                "Harten entropy fix; OneFLOW NsInvFluxPolicy",
                "doc/source/numerical/harten_entropy_fix.rst",
            ),
            csv_relpath=None,
            plot_module=None,
            pure_apis=("ApplyHartenEntropyFix",),
            source_relpaths=(
                "codes/ns/include/NsInvFluxPolicy.h",
                "codes/ns/src/NsInvFluxPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/index.rst",
            code_sketch=(
                "// GATED SKETCH - mirror production MAX in tests.\n"
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


def validate_catalog_paths(repository_root: Path) -> list[str]:
    """Return missing repository files referenced by the theory catalog."""
    issues: list[str] = []
    for entry in list_theory_entries():
        if bool(entry.csv_relpath) != bool(entry.plot_module):
            issues.append(
                f"{entry.id}: CSV and plot module must be configured together"
            )
        for relpath in entry.source_relpaths:
            path = resolve_repo_path(repository_root, relpath)
            if path is None or not path.is_file():
                issues.append(f"{entry.id}: missing source file: {relpath}")
        if entry.doc_relpath:
            path = doc_path(repository_root, entry)
            if path is None or not path.is_file():
                issues.append(f"{entry.id}: missing documentation: {entry.doc_relpath}")
        if entry.csv_relpath:
            path = csv_path(repository_root, entry)
            if path is None or not path.is_file():
                issues.append(f"{entry.id}: missing theory data: {entry.csv_relpath}")
        if entry.plot_module:
            plot_path = verification_root(repository_root) / f"{entry.plot_module}.py"
            if not plot_path.is_file():
                issues.append(
                    f"{entry.id}: missing plot module: "
                    f"tools/verification/{entry.plot_module}.py"
                )
    return issues


def validate_catalog_api_references(repository_root: Path) -> list[str]:
    """Return pure API symbols not referenced by an entry's source files."""
    issues: list[str] = []
    for entry in list_theory_entries():
        source_text = "\n".join(
            path.read_text(encoding="utf-8", errors="replace")
            for path in source_paths(repository_root, entry)
            if path.is_file()
        )
        for api in entry.pure_apis:
            if re.search(rf"\b{re.escape(api)}\s*\(", source_text) is None:
                issues.append(
                    f"{entry.id}: API reference not found in source files: {api}"
                )
    return issues


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
