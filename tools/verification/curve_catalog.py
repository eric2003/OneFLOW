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
    examples: tuple[str, ...] = ()


def list_theory_entries() -> list[TheoryEntry]:
    return [
        # ========== SA ==========


        TheoryEntry(
            id="sa_ft2",
            family="turbulence",
            genus="SA",
            section="closure / trip function",
            title="SA trip function ft2",
            role_tags=("turbulence", "SA", "closure"),
            formula_plain=(
                "ft2 = 0 if disabled; else ct3 * exp(-ct4 * chi^2)"
            ),
            formula_latex=(
                r"f_{t2}=c_{t3}\,e^{-c_{t4}\chi^{2}}"
                r"\quad(\text{enabled})"
            ),
            principle=(
                "Optional trip term that reduces production and adjusts "
                "destruction near laminar/transition regions."
            ),
            background=(
                "Pure API SaFt2 in TurbSaPolicy. Feeds SaProduction and "
                "SaDestruction when the trip flag is on."
            ),
            references=(
                "Spalart-Allmaras model; NASA TMR SA",
                "doc/source/numerical/sa_ft2.rst",
            ),
            csv_relpath="tools/verification/data/sa_ft2_theory.csv",
            plot_module="plot_sa_ft2",
            pure_apis=("SaFt2",),
            source_relpaths=(
                "codes/turb/include/TurbSaPolicy.h",
                "codes/turb/src/TurbSaPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sa_ft2.rst",
            code_sketch=(
                "// GATED SKETCH enabled\n"
                "ft2 = ct3 * exp(-ct4 * chi * chi);\n"
            ),
            proves="ft2 algebra vs pure API on a chi grid; disabled path returns 0.",
            does_not_prove="Transition prediction quality on a case mesh.",
        ),



        TheoryEntry(
            id="sa_working_variable_upper_bound",
            family="turbulence",
            genus="SA",
            section="update / working-variable bound",
            title="SA working-variable upper bound",
            role_tags=("turbulence", "SA", "update"),
            formula_plain=(
                "nu_max = rmax * mu_laminar / (rho * fv1)"
            ),
            formula_latex=(
                r"\tilde{\nu}_{\max}"
                r"=r_{\max}\frac{\mu_\ell}{\rho f_{v1}}"
            ),
            principle=(
                "Caps the SA working variable so the implied eddy viscosity "
                "ratio stays below a prescribed maximum."
            ),
            background=(
                "Pure API SaWorkingVariableUpperBound in TurbSaPolicy. Pair "
                "with SaNonNegativeWorkingVariable and SaUnderRelaxedUpdate."
            ),
            references=(
                "SA solver robustness practice; NASA TMR SA",
                "doc/source/numerical/sa_working_variable_upper_bound.rst",
            ),
            csv_relpath="tools/verification/data/sa_working_variable_upper_bound_theory.csv",
            plot_module="plot_sa_working_variable_upper_bound",
            pure_apis=("SaWorkingVariableUpperBound",),
            source_relpaths=(
                "codes/turb/include/TurbSaPolicy.h",
                "codes/turb/src/TurbSaPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sa_working_variable_upper_bound.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "nu_max = rmax * mu / (rho * fv1);\n"
            ),
            proves="Upper-bound algebra vs pure API on an fv1 grid.",
            does_not_prove="Optimal rmax or full update pipeline on a case.",
        ),

        TheoryEntry(
            id="sa_under_relaxed_update",
            family="turbulence",
            genus="SA",
            section="update / under-relaxation",
            title="SA under-relaxed update",
            role_tags=("turbulence", "SA", "update"),
            formula_plain=(
                "q = (1-alpha)*qOld + alpha*(qOld+deltaQ) = qOld + alpha*deltaQ"
            ),
            formula_latex=(
                r"q=q_0+\alpha\Delta q"
            ),
            principle=(
                "Under-relaxation blends the explicit update with the old "
                "working variable to improve nonlinear stability."
            ),
            background=(
                "Pure API SaUnderRelaxedUpdate in TurbSaPolicy."
            ),
            references=(
                "Common CFD under-relaxation practice; SA solver update path",
                "doc/source/numerical/sa_under_relaxed_update.rst",
            ),
            csv_relpath="tools/verification/data/sa_under_relaxed_update_theory.csv",
            plot_module="plot_sa_under_relaxed_update",
            pure_apis=("SaUnderRelaxedUpdate",),
            source_relpaths=(
                "codes/turb/include/TurbSaPolicy.h",
                "codes/turb/src/TurbSaPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sa_under_relaxed_update.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "q = qOld + alpha * deltaQ;\n"
            ),
            proves="Under-relaxation algebra vs pure API on an alpha grid.",
            does_not_prove="Optimal alpha or case convergence behavior.",
        ),

        TheoryEntry(
            id="sa_eddy_viscosity",
            family="turbulence",
            genus="SA",
            section="closure / eddy viscosity",
            title="SA eddy viscosity",
            role_tags=("turbulence", "SA", "closure"),
            formula_plain=(
                "mut = density * nu_tilde * fv1"
            ),
            formula_latex=(
                r"\mu_t=\rho\,\tilde{\nu}\,f_{v1}"
            ),
            principle=(
                "Final SA turbulent viscosity from density, working variable "
                "and the fv1 wall-damping function."
            ),
            background=(
                "Pure API SaEddyViscosity in TurbSaPolicy. Upstream: SaChi, "
                "SaFv1. Full pipeline: SaComputeEddyViscosity."
            ),
            references=(
                "Spalart-Allmaras model; NASA TMR SA",
                "doc/source/numerical/sa_eddy_viscosity.rst",
            ),
            csv_relpath="tools/verification/data/sa_eddy_viscosity_theory.csv",
            plot_module="plot_sa_eddy_viscosity",
            pure_apis=("SaEddyViscosity",),
            source_relpaths=(
                "codes/turb/include/TurbSaPolicy.h",
                "codes/turb/src/TurbSaPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sa_eddy_viscosity.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "mut = density * nu_tilde * fv1;\n"
            ),
            proves="Eddy-viscosity algebra vs pure API; composition with SaFv1.",
            does_not_prove="Full capped SaComputeEddyViscosity pipeline on a case.",
        ),

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
                "// GATED SKETCH - review before use.\n"
                "Real SaFv1(Real chi, Real cv1Cubed) {\n"
                "  Real chi3 = chi*chi*chi;\n"
                "  return chi3 / (chi3 + cv1Cubed);\n"
                "}\n"
            ),
            proves="Closed-form fv1 vs SaFv1.",
            does_not_prove="Full residual or experimental Cf.",
            examples=(
                "Input: chi = 7.1 and cv1 = 7.1.\n"
                "Calculation: fv1 = chi^3 / (chi^3 + cv1^3) = 1 / 2.\n"
                "Result: fv1 = 0.5, so mu_t = 0.5 * rho * nu_tilde.\n"
                "Limit: this checks the damping function only, not a wall-resolved CFD solution.",
            ),
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
            id="sa_s_bar",
            family="turbulence",
            genus="SA",
            section="closure / modified strain intermediate",
            title="SA S-bar intermediate",
            role_tags=("turbulence", "SA", "closure"),
            formula_plain=(
                "sBar = rs * fv2"
            ),
            formula_latex=(
                r"\bar{S}=r_s f_{v2}"
            ),
            principle=(
                "Intermediate strain correction formed from the wall-scaled "
                "working-variable measure and fv2, before assembling the "
                "modified strain rate used in production."
            ),
            background=(
                "Pure API SaSBar in TurbSaPolicy. Upstream: SaFv2. Downstream: "
                "SaModifiedStrainRate."
            ),
            references=(
                "Spalart-Allmaras model; NASA TMR SA",
                "doc/source/numerical/sa_s_bar.rst",
            ),
            csv_relpath="tools/verification/data/sa_s_bar_theory.csv",
            plot_module="plot_sa_s_bar",
            pure_apis=("SaSBar",),
            source_relpaths=(
                "codes/turb/include/TurbSaPolicy.h",
                "codes/turb/src/TurbSaPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sa_s_bar.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "sBar = rs * fv2;\n"
            ),
            proves="S-bar algebra vs pure API; composition with SaFv2.",
            does_not_prove="Full modified-strain / production field on a case.",
        ),

        TheoryEntry(
            id="sa_r_ratio",
            family="turbulence",
            genus="SA",
            section="closure / wall destruction ratio",
            title="SA r-ratio",
            role_tags=("turbulence", "SA", "closure"),
            formula_plain=(
                "r = min(rs / modifiedStrain, 10)"
            ),
            formula_latex=(
                r"r=\min\!\left(\frac{r_s}{\tilde{S}},\,10\right)"
            ),
            principle=(
                "Normalized ratio that feeds the SA wall-destruction function "
                "fw through g(r). Clamped at 10 for numerical robustness."
            ),
            background=(
                "Pure API SaRRatio in TurbSaPolicy. Downstream: SaG, SaFw."
            ),
            references=(
                "Spalart-Allmaras model; NASA TMR SA",
                "doc/source/numerical/sa_r_ratio.rst",
            ),
            csv_relpath="tools/verification/data/sa_r_ratio_theory.csv",
            plot_module="plot_sa_r_ratio",
            pure_apis=("SaRRatio",),
            source_relpaths=(
                "codes/turb/include/TurbSaPolicy.h",
                "codes/turb/src/TurbSaPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sa_r_ratio.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "r = min(rs / S, 10);\n"
            ),
            proves="r-ratio algebra and clamp vs pure API.",
            does_not_prove="Full fw destruction on a wall-bounded case.",
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
            id="sa_limited_diffusion",
            family="turbulence",
            genus="SA",
            section="source / limited diffusion",
            title="SA limited diffusion",
            role_tags=("turbulence", "SA", "source"),
            formula_plain=(
                "D_lim = min(rprod * P, Diff)"
            ),
            formula_latex=(
                r"D_{\mathrm{lim}}=\min(r_p P,\,D_{\mathrm{diff}})"
            ),
            principle=(
                "Caps the non-conservative diffusion contribution relative to "
                "production to improve robustness."
            ),
            background=(
                "Pure API SaLimitedDiffusion in TurbSaPolicy. Upstream: "
                "SaDiffusion, SaProduction. Downstream: SaSourceResidual."
            ),
            references=(
                "Spalart-Allmaras implementation practice; NASA TMR SA",
                "doc/source/numerical/sa_limited_diffusion.rst",
            ),
            csv_relpath="tools/verification/data/sa_limited_diffusion_theory.csv",
            plot_module="plot_sa_limited_diffusion",
            pure_apis=("SaLimitedDiffusion",),
            source_relpaths=(
                "codes/turb/include/TurbSaPolicy.h",
                "codes/turb/src/TurbSaPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sa_limited_diffusion.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "D_lim = min(rprod * P, Diff);\n"
            ),
            proves="Limiter algebra vs pure API; composition with SaDiffusion.",
            does_not_prove="Full SA residual field on a wall-bounded case.",
        ),

        TheoryEntry(
            id="sa_diffusion",
            family="turbulence",
            genus="SA",
            section="source / diffusion",
            title="SA diffusion source term",
            role_tags=("turbulence", "SA", "source"),
            formula_plain=(
                "Diff = (cb2/sigma) * |grad nu_tilde|^2 * invRe"
            ),
            formula_latex=(
                r"D_{\mathrm{diff}}"
                r"=\frac{c_{b2}}{\sigma}|\nabla\tilde{\nu}|^{2}/\mathrm{Re}"
            ),
            principle=(
                "SA non-conservative diffusion contribution proportional to "
                "the squared gradient of the working variable."
            ),
            background=(
                "Pure API SaDiffusion in TurbSaPolicy. Downstream: "
                "SaLimitedDiffusion then SaSourceResidual."
            ),
            references=(
                "Spalart-Allmaras model; NASA TMR SA",
                "doc/source/numerical/sa_diffusion.rst",
            ),
            csv_relpath="tools/verification/data/sa_diffusion_theory.csv",
            plot_module="plot_sa_diffusion",
            pure_apis=("SaDiffusion",),
            source_relpaths=(
                "codes/turb/include/TurbSaPolicy.h",
                "codes/turb/src/TurbSaPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sa_diffusion.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "Diff = cb2OverSigma * gradNuSquared * invRe;\n"
            ),
            proves="Diffusion algebra vs pure API on a |grad nu|^2 grid.",
            does_not_prove="Full SA residual or limiter interaction on a case.",
        ),

        TheoryEntry(
            id="sa_destruction",
            family="turbulence",
            genus="SA",
            section="source / destruction",
            title="SA destruction term D",
            role_tags=("turbulence", "SA", "source"),
            formula_plain=(
                "D = (cw1k * fw - cb1 * ft2) * nu_tilde_rs"
            ),
            formula_latex=(
                r"D=\bigl(c_{w1}f_w-c_{b1}f_{t2}\bigr)\tilde{\nu}_{rs}"
            ),
            principle=(
                "SA wall destruction balances production away from the wall. "
                "fw is the wall destruction function; ft2 is the optional trip."
            ),
            background=(
                "Pure API SaDestruction in TurbSaPolicy. Feeds SaSourceResidual."
            ),
            references=(
                "Spalart-Allmaras model; NASA TMR SA",
                "doc/source/numerical/sa_destruction.rst",
            ),
            csv_relpath="tools/verification/data/sa_destruction_theory.csv",
            plot_module="plot_sa_destruction",
            pure_apis=("SaDestruction",),
            source_relpaths=(
                "codes/turb/include/TurbSaPolicy.h",
                "codes/turb/src/TurbSaPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sa_destruction.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "D = (cw1k * fw - cb1 * ft2) * nu_tilde_rs;\n"
            ),
            proves="Destruction algebra vs pure API on a nu_rs grid.",
            does_not_prove="Full SA residual or near-wall case accuracy.",
        ),

        TheoryEntry(
            id="sa_production",
            family="turbulence",
            genus="SA",
            section="source / production",
            title="SA production term P",
            role_tags=("turbulence", "SA", "source"),
            formula_plain=(
                "P = cb1 * (1 - ft2) * S_tilde * nu_tilde"
            ),
            formula_latex=(
                r"P=c_{b1}(1-f_{t2})\tilde{S}\,\tilde{\nu}"
            ),
            principle=(
                "SA turbulence production is proportional to the modified "
                "strain rate and working variable, reduced by the optional "
                "ft2 trip term."
            ),
            background=(
                "Pure API SaProduction in TurbSaPolicy. Feeds SaSourceResidual "
                "with destruction and limited diffusion."
            ),
            references=(
                "Spalart-Allmaras model; NASA TMR SA",
                "doc/source/numerical/sa_production.rst",
            ),
            csv_relpath="tools/verification/data/sa_production_theory.csv",
            plot_module="plot_sa_production",
            pure_apis=("SaProduction",),
            source_relpaths=(
                "codes/turb/include/TurbSaPolicy.h",
                "codes/turb/src/TurbSaPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sa_production.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "P = cb1 * (1 - ft2) * S_tilde * nu_tilde;\n"
            ),
            proves="Production algebra vs pure API on a nu_tilde grid.",
            does_not_prove="Full SA residual on a wall-bounded case.",
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
            csv_relpath="tools/verification/data/sa_source_residual_theory.csv",
            plot_module="plot_sa_source_residual",
            pure_apis=(
                "SaProduction",
                "SaDestruction",
                "SaDiffusion",
                "SaLimitedDiffusion",
                "SaSourceResidual",
            ),
            source_relpaths=(
                "codes/turb/include/TurbSaPolicy.h",
                "codes/turb/src/TurbSaPolicy.cpp",
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
            id="sst_floor_omega",
            family="turbulence",
            genus="SST",
            section="update / positivity floor",
            title="SST omega floor",
            role_tags=("turbulence", "SST", "update"),
            formula_plain=(
                "omega = max(omega, omega_min)"
            ),
            formula_latex=(
                r"\omega\leftarrow\max(\omega,\omega_{\min})"
            ),
            principle=(
                "Enforces a positive lower bound on specific dissipation to "
                "avoid division by zero and unphysical negative values."
            ),
            background=(
                "Pure API SstFloorOmega in TurbSstPolicy. Pair with "
                "SstFloorKineticEnergy."
            ),
            references=(
                "SST solver robustness practice; NASA TMR SST",
                "doc/source/numerical/sst_floor_omega.rst",
            ),
            csv_relpath="tools/verification/data/sst_floor_omega_theory.csv",
            plot_module="plot_sst_floor_omega",
            pure_apis=("SstFloorOmega",),
            source_relpaths=(
                "codes/turb/include/TurbSstPolicy.h",
                "codes/turb/src/TurbSstPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sst_floor_omega.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "omega = max(omega, omegaMin);\n"
            ),
            proves="Floor algebra vs pure API on an omega grid.",
            does_not_prove="Optimal omega_min or case stability margins.",
        ),

        TheoryEntry(
            id="sst_blend",
            family="turbulence",
            genus="SST",
            section="closure / coefficient blend",
            title="SST coefficient blend",
            role_tags=("turbulence", "SST", "closure"),
            formula_plain=(
                "phi = F1 * phi1 + (1 - F1) * phi2"
            ),
            formula_latex=(
                r"\phi=F_1\phi_1+(1-F_1)\phi_2"
            ),
            principle=(
                "Convex combination of inner (Wilcox/k-omega) and outer "
                "(k-epsilon) model coefficients driven by the F1 blending "
                "function."
            ),
            background=(
                "Pure API SstBlend in TurbSstPolicy. Used for beta, sigma, "
                "alpha_w and related SST constants."
            ),
            references=(
                "Menter SST (1994); NASA TMR SST",
                "doc/source/numerical/sst_blend.rst",
            ),
            csv_relpath="tools/verification/data/sst_blend_theory.csv",
            plot_module="plot_sst_blend",
            pure_apis=("SstBlend",),
            source_relpaths=(
                "codes/turb/include/TurbSstPolicy.h",
                "codes/turb/src/TurbSstPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sst_blend.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "phi = F1 * phi1 + (1 - F1) * phi2;\n"
            ),
            proves="Blend algebra and convexity on an F1 grid.",
            does_not_prove="F1 construction accuracy or full residual field.",
        ),

        TheoryEntry(
            id="sst_f1",
            family="turbulence",
            genus="SST",
            section="closure / blending",
            title="F1 blending function",
            role_tags=("closure", "blending"),
            formula_plain=(
                "CD_kw = max(2*rho*sigw2*(gradk*gradw)/w, small)\n"
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
            examples=(
                "Input: arg1 = 1, phi1 = 2, and phi2 = 1 (dimensionless illustrative values).\\n"
                "Calculation: F1 = tanh(arg1^4) = tanh(1) = 0.761594; "
                "phi = F1*phi1 + (1-F1)*phi2 = 1.761594.\\n"
                "Result: the blend is weighted toward the inner-state value phi1 for this input.\\n"
                "Limit: this evaluates the blending algebra only; it does not derive arg1 from a mesh flow state or establish full SST case accuracy.",
            ),
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
            examples=(
                "Input: arg1 = 1.\n"
                "Calculation: F1 = tanh(arg1^4) = tanh(1).\n"
                "Result: F1 is approximately 0.761594.\n"
                "Interpretation: the blend is continuous; this one point does not validate "
                "the complete SST closure or its wall behavior.",
            ),
        ),





        TheoryEntry(
            id="sst_prodk_vorticity",
            family="turbulence",
            genus="SST",
            section="source / k production",
            title="SST vorticity-based k production Pk",
            role_tags=("turbulence", "SST", "source"),
            formula_plain=(
                "Pk = mut * Omega2 * invRe - 2/3 * rho * k * div"
            ),
            formula_latex=(
                r"P_k=\mu_t\Omega^{2}/\mathrm{Re}"
                r"-\tfrac23\rho k(\nabla\cdot u)"
            ),
            principle=(
                "Production estimated from vorticity magnitude instead of "
                "strain rate, with the same divergence correction as the "
                "Boussinesq form."
            ),
            background=(
                "Pure API SstProdkVorticity in TurbSstPolicy. Third Pk option "
                "alongside Boussinesq and strain-only."
            ),
            references=(
                "SST production-option variants; NASA TMR SST",
                "doc/source/numerical/sst_prodk_vorticity.rst",
            ),
            csv_relpath="tools/verification/data/sst_prodk_vorticity_theory.csv",
            plot_module="plot_sst_prodk_vorticity",
            pure_apis=("SstProdkVorticity",),
            source_relpaths=(
                "codes/turb/include/TurbSstPolicy.h",
                "codes/turb/src/TurbSstPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sst_prodk_vorticity.rst",
            code_sketch=(
                "// GATED SKETCH div=0\n"
                "Pk = mut * Omega2 * invRe;\n"
            ),
            proves="Vorticity Pk algebra vs pure API; divergence sample.",
            does_not_prove="Which production option is best for a case.",
        ),

        TheoryEntry(
            id="sst_prodk_strain_only",
            family="turbulence",
            genus="SST",
            section="source / k production",
            title="SST strain-only k production Pk",
            role_tags=("turbulence", "SST", "source"),
            formula_plain=(
                "Pk = mut * S2 * invReynolds"
            ),
            formula_latex=(
                r"P_k=\mu_t S^{2}/\mathrm{Re}"
            ),
            principle=(
                "Simplified production using only strain-rate magnitude, "
                "without divergence correction. Alternate to Boussinesq form."
            ),
            background=(
                "Pure API SstProdkStrainOnly in TurbSstPolicy. Pair with "
                "SstProdkBoussinesq and SstProdkVorticity options."
            ),
            references=(
                "SST production-option variants; NASA TMR SST",
                "doc/source/numerical/sst_prodk_strain_only.rst",
            ),
            csv_relpath="tools/verification/data/sst_prodk_strain_only_theory.csv",
            plot_module="plot_sst_prodk_strain_only",
            pure_apis=("SstProdkStrainOnly",),
            source_relpaths=(
                "codes/turb/include/TurbSstPolicy.h",
                "codes/turb/src/TurbSstPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sst_prodk_strain_only.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "Pk = mut * S2 * invRe;\n"
            ),
            proves="Strain-only Pk algebra vs pure API on an S2 grid.",
            does_not_prove="Which production option is best for a case.",
        ),

        TheoryEntry(
            id="sst_prodk_boussinesq",
            family="turbulence",
            genus="SST",
            section="source / k production",
            title="SST Boussinesq k production Pk",
            role_tags=("turbulence", "SST", "source"),
            formula_plain=(
                "Pk = mut * (S2 - 2/3 div^2) / Re - 2/3 * rho * k * div"
            ),
            formula_latex=(
                r"P_k=\mu_t\bigl(S^{2}-\tfrac23(\nabla\cdot u)^{2}\bigr)/\mathrm{Re}"
                r"-\tfrac23\rho k(\nabla\cdot u)"
            ),
            principle=(
                "Boussinesq production from eddy viscosity and mean strain, "
                "with divergence corrections for compressible flow."
            ),
            background=(
                "Pure API SstProdkBoussinesq in TurbSstPolicy. Alternatives: "
                "SstProdkStrainOnly, SstProdkVorticity. Downstream: SstLimitProdk."
            ),
            references=(
                "Menter SST; standard Boussinesq production forms",
                "doc/source/numerical/sst_prodk_boussinesq.rst",
            ),
            csv_relpath="tools/verification/data/sst_prodk_boussinesq_theory.csv",
            plot_module="plot_sst_prodk_boussinesq",
            pure_apis=("SstProdkBoussinesq",),
            source_relpaths=(
                "codes/turb/include/TurbSstPolicy.h",
                "codes/turb/src/TurbSstPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sst_prodk_boussinesq.rst",
            code_sketch=(
                "// GATED SKETCH div=0\n"
                "Pk = mut * S2 * invRe;\n"
            ),
            proves="Boussinesq Pk algebra vs pure API; divergence sample.",
            does_not_prove="Full SST residual or limiter interaction on a case.",
        ),

        TheoryEntry(
            id="sst_limit_prodk",
            family="turbulence",
            genus="SST",
            section="source / k production limiter",
            title="SST k-production limiter",
            role_tags=("turbulence", "SST", "source"),
            formula_plain=(
                "P_k <- min(P_k, pklim * eps)"
            ),
            formula_latex=(
                r"P_k\leftarrow\min(P_k,\,p_{k\lim}\,\varepsilon)"
            ),
            principle=(
                "Menter-style production limiter prevents excessive P_k relative "
                "to dissipation, improving robustness in stagnation regions."
            ),
            background=(
                "Pure API SstLimitProdk in TurbSstPolicy. eps from "
                "SstDisskStandard (or DES variant)."
            ),
            references=(
                "Menter SST production limiter; NASA TMR SST",
                "doc/source/numerical/sst_limit_prodk.rst",
            ),
            csv_relpath="tools/verification/data/sst_limit_prodk_theory.csv",
            plot_module="plot_sst_limit_prodk",
            pure_apis=("SstLimitProdk",),
            source_relpaths=(
                "codes/turb/include/TurbSstPolicy.h",
                "codes/turb/src/TurbSstPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sst_limit_prodk.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "P = min(P, pklim * eps);\n"
            ),
            proves="Limiter algebra vs pure API; composition with SstDisskStandard.",
            does_not_prove="Full SST k residual on a case mesh.",
        ),


        TheoryEntry(
            id="sst_dissk_des",
            family="turbulence",
            genus="SST",
            section="source / k dissipation DES",
            title="SST DES k dissipation",
            role_tags=("turbulence", "SST", "DES", "source"),
            formula_plain=(
                "eps = f_beta_star * rho * k * sqrt(k) / L"
            ),
            formula_latex=(
                r"\varepsilon_{\mathrm{DES}}"
                r"=f_{\beta^*}\rho\,k^{3/2}/L"
            ),
            principle=(
                "DES replaces the standard omega-based k destruction with a "
                "length-scale form based on local filter/wall distance measure L."
            ),
            background=(
                "Pure API SstDisskDes in TurbSstPolicy. Alternate to "
                "SstDisskStandard for DES-mode branches."
            ),
            references=(
                "Spalart DES; Menter SST-DES variants; NASA TMR",
                "doc/source/numerical/sst_dissk_des.rst",
            ),
            csv_relpath="tools/verification/data/sst_dissk_des_theory.csv",
            plot_module="plot_sst_dissk_des",
            pure_apis=("SstDisskDes",),
            source_relpaths=(
                "codes/turb/include/TurbSstPolicy.h",
                "codes/turb/src/TurbSstPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sst_dissk_des.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "eps = fBetas * density * k * sqrt(k) / L;\n"
            ),
            proves="DES dissipation algebra vs pure API on a k grid.",
            does_not_prove="Full DES residual or length-scale construction.",
        ),

        TheoryEntry(
            id="sst_dissk_standard",
            family="turbulence",
            genus="SST",
            section="source / k dissipation",
            title="SST standard k dissipation eps",
            role_tags=("turbulence", "SST", "source"),
            formula_plain=(
                "eps = f_beta_star * beta_star * rho * k * omega"
            ),
            formula_latex=(
                r"\varepsilon=f_{\beta^*}\beta^*\,\rho\,k\,\omega"
            ),
            principle=(
                "Standard SST/k-omega destruction of turbulent kinetic energy. "
                "f_beta* is a Wilcox-style compressibility or low-Re modifier "
                "when enabled (often 1)."
            ),
            background=(
                "Pure API SstDisskStandard in TurbSstPolicy. Used with "
                "SstLimitProdk and SstSourceK for the k equation source."
            ),
            references=(
                "Menter SST; Wilcox k-omega; NASA TMR SST",
                "doc/source/numerical/sst_dissk_standard.rst",
            ),
            csv_relpath="tools/verification/data/sst_dissk_standard_theory.csv",
            plot_module="plot_sst_dissk_standard",
            pure_apis=("SstDisskStandard",),
            source_relpaths=(
                "codes/turb/include/TurbSstPolicy.h",
                "codes/turb/src/TurbSstPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sst_dissk_standard.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "eps = fBetas * betas * density * k * omega;\n"
            ),
            proves="k-dissipation algebra vs pure API on an omega grid.",
            does_not_prove="Full SST k residual on a case mesh.",
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
            csv_relpath="tools/verification/data/sst_source_k_theory.csv",
            plot_module="plot_sst_source_k",
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
            id="sst_dissw",
            family="turbulence",
            genus="SST",
            section="source / omega dissipation",
            title="SST omega dissipation Dw",
            role_tags=("turbulence", "SST", "source"),
            formula_plain=(
                "Dw = f_beta * beta * rho * omega^2"
            ),
            formula_latex=(
                r"D_\omega=f_\beta\,\beta\,\rho\,\omega^{2}"
            ),
            principle=(
                "SST/k-omega destruction of specific dissipation rate. "
                "f_beta is a Wilcox-style modifier when enabled (often 1)."
            ),
            background=(
                "Pure API SstDissw in TurbSstPolicy. Feeds SstSourceW with "
                "production and cross-diffusion."
            ),
            references=(
                "Menter SST; Wilcox k-omega; NASA TMR SST",
                "doc/source/numerical/sst_dissw.rst",
            ),
            csv_relpath="tools/verification/data/sst_dissw_theory.csv",
            plot_module="plot_sst_dissw",
            pure_apis=("SstDissw",),
            source_relpaths=(
                "codes/turb/include/TurbSstPolicy.h",
                "codes/turb/src/TurbSstPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sst_dissw.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "Dw = fBeta * beta * density * omega * omega;\n"
            ),
            proves="omega-dissipation algebra vs pure API on an omega grid.",
            does_not_prove="Full SST omega residual on a case mesh.",
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
            csv_relpath="tools/verification/data/sst_source_w_theory.csv",
            plot_module="plot_sst_source_w",
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
            id="sst_prodw_wilcox",
            family="turbulence",
            genus="SST",
            section="source / omega production",
            title="SST Wilcox omega production Pw",
            role_tags=("turbulence", "SST", "source"),
            formula_plain=(
                "Pw = alphaw * omega / k * Pk"
            ),
            formula_latex=(
                r"P_\omega=\alpha_\omega\frac{\omega}{k}P_k"
            ),
            principle=(
                "Wilcox-style omega production scales limited k-production by "
                "omega/k rather than density/eddy-viscosity."
            ),
            background=(
                "Pure API SstProdwWilcox in TurbSstPolicy. Alternate to "
                "SstProdwMenter. Feeds SstSourceW."
            ),
            references=(
                "Wilcox k-omega; Menter SST production options; NASA TMR",
                "doc/source/numerical/sst_prodw_wilcox.rst",
            ),
            csv_relpath="tools/verification/data/sst_prodw_wilcox_theory.csv",
            plot_module="plot_sst_prodw_wilcox",
            pure_apis=("SstProdwWilcox",),
            source_relpaths=(
                "codes/turb/include/TurbSstPolicy.h",
                "codes/turb/src/TurbSstPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sst_prodw_wilcox.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "Pw = alphaw * omega / k * Pk;\n"
            ),
            proves="Wilcox Pw algebra vs pure API on a Pk grid.",
            does_not_prove="Full SST residual or Menter branch equivalence.",
        ),

        TheoryEntry(
            id="sst_prodw_menter",
            family="turbulence",
            genus="SST",
            section="source / omega production",
            title="SST Menter omega production Pw",
            role_tags=("turbulence", "SST", "source"),
            formula_plain=(
                "Pw = alphaw * rho / mut * Pk * Re"
            ),
            formula_latex=(
                r"P_\omega=\alpha_\omega\frac{\rho}{\mu_t}P_k\,\mathrm{Re}"
            ),
            principle=(
                "Menter form converts limited k-production into omega "
                "production using eddy viscosity scaling."
            ),
            background=(
                "Pure API SstProdwMenter in TurbSstPolicy. Feeds SstSourceW. "
                "Wilcox form SstProdwWilcox is an alternate branch."
            ),
            references=(
                "Menter SST (1994); NASA TMR SST",
                "doc/source/numerical/sst_prodw_menter.rst",
            ),
            csv_relpath="tools/verification/data/sst_prodw_menter_theory.csv",
            plot_module="plot_sst_prodw_menter",
            pure_apis=("SstProdwMenter",),
            source_relpaths=(
                "codes/turb/include/TurbSstPolicy.h",
                "codes/turb/src/TurbSstPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sst_prodw_menter.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "Pw = alphaw * density / mut * Pk * reynolds;\n"
            ),
            proves="Menter Pw algebra vs pure API on a Pk grid.",
            does_not_prove="Full SST residual or Wilcox production branch.",
        ),


        TheoryEntry(
            id="sst_cross_diff_wilcox",
            family="turbulence",
            genus="SST",
            section="source / cross diffusion",
            title="SST Wilcox cross-diffusion CD",
            role_tags=("turbulence", "SST", "source"),
            formula_plain=(
                "CD = sigd * rho * max(crossTerm, 0) / omega"
            ),
            formula_latex=(
                r"\mathrm{CD}_W=\sigma_d\rho\frac{\max(\nabla k\cdot\nabla\omega,0)}{\omega}"
            ),
            principle=(
                "Wilcox cross-diffusion activates only for positive "
                "grad-k · grad-omega and does not use the F1 blend factor."
            ),
            background=(
                "Pure API SstCrossDiffWilcox in TurbSstPolicy. Alternate to "
                "SstCrossDiffMenter."
            ),
            references=(
                "Wilcox k-omega cross-diffusion; NASA TMR",
                "doc/source/numerical/sst_cross_diff_wilcox.rst",
            ),
            csv_relpath="tools/verification/data/sst_cross_diff_wilcox_theory.csv",
            plot_module="plot_sst_cross_diff_wilcox",
            pure_apis=("SstCrossDiffWilcox",),
            source_relpaths=(
                "codes/turb/include/TurbSstPolicy.h",
                "codes/turb/src/TurbSstPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sst_cross_diff_wilcox.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "CD = sigd * density * max(cross, 0) / omega;\n"
            ),
            proves="Wilcox CD algebra vs pure API; negative cross clamp.",
            does_not_prove="Full residual field or Menter/Wilcox switch logic.",
        ),

        TheoryEntry(
            id="sst_cross_diff_menter",
            family="turbulence",
            genus="SST",
            section="source / cross diffusion",
            title="SST Menter cross-diffusion CD",
            role_tags=("turbulence", "SST", "source"),
            formula_plain=(
                "CD = 2 * (1 - F1) * rho * sigw2 * (grad k · grad omega) / omega"
            ),
            formula_latex=(
                r"\mathrm{CD}=2(1-F_1)\rho\sigma_{\omega2}"
                r"\frac{\nabla k\cdot\nabla\omega}{\omega}"
            ),
            principle=(
                "Cross-diffusion blends Wilcox and k-epsilon behaviors. "
                "It is active in the outer layer (F1 near 0) and vanishes "
                "near walls when F1 approaches 1."
            ),
            background=(
                "Pure API SstCrossDiffMenter in TurbSstPolicy. Feeds "
                "SstSourceW and SstDiagonalW."
            ),
            references=(
                "Menter SST (1994); NASA TMR SST",
                "doc/source/numerical/sst_cross_diff_menter.rst",
            ),
            csv_relpath="tools/verification/data/sst_cross_diff_menter_theory.csv",
            plot_module="plot_sst_cross_diff_menter",
            pure_apis=("SstCrossDiffMenter",),
            source_relpaths=(
                "codes/turb/include/TurbSstPolicy.h",
                "codes/turb/src/TurbSstPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sst_cross_diff_menter.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "CD = 2*(1-F1)*rho*sigw2*cross/omega;\n"
            ),
            proves="Cross-diffusion algebra vs pure API; F1=1 vanishing.",
            does_not_prove="Full SST residual field or F1 construction accuracy.",
        ),


        TheoryEntry(
            id="sst_negative_source_part",
            family="turbulence",
            genus="SST",
            section="residual / negative source part",
            title="SST negative source part",
            role_tags=("turbulence", "SST", "residual"),
            formula_plain=(
                "n = 0.5 * (S - |S|) = min(S, 0)"
            ),
            formula_latex=(
                r"n=\min(S,0)"
            ),
            principle=(
                "Extracts only the non-positive part of a source residual for "
                "implicit linearization / diagonal contributions."
            ),
            background=(
                "Pure API SstNegativeSourcePart in TurbSstPolicy. Used with "
                "SstSourceK / SstSourceW residual assembly."
            ),
            references=(
                "Menter SST linearization practice",
                "doc/source/numerical/sst_negative_source_part.rst",
            ),
            csv_relpath="tools/verification/data/sst_negative_source_part_theory.csv",
            plot_module="plot_sst_negative_source_part",
            pure_apis=("SstNegativeSourcePart",),
            source_relpaths=(
                "codes/turb/include/TurbSstPolicy.h",
                "codes/turb/src/TurbSstPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sst_negative_source_part.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "n = (S < 0) ? S : 0;\n"
            ),
            proves="Negative-part algebra vs pure API; composition with SstSourceK.",
            does_not_prove="Full residual Jacobian on a case mesh.",
        ),

        TheoryEntry(
            id="sst_diagonal_w",
            family="turbulence",
            genus="SST",
            section="residual / omega diagonal",
            title="SST omega residual diagonal",
            role_tags=("turbulence", "SST", "residual"),
            formula_plain=(
                "dw = -2 * f_beta * beta * omega - |CD| / (rho * omega)"
            ),
            formula_latex=(
                r"d_\omega=-2f_\beta\beta\omega"
                r"-\frac{|\mathrm{CD}|}{\rho\omega}"
            ),
            principle=(
                "Linearized omega-equation diagonal used by implicit updates. "
                "Cross-diffusion adds a non-negative magnitude term when active."
            ),
            background=(
                "Pure API SstDiagonalW in TurbSstPolicy. Pair with SstDiagonalK."
            ),
            references=(
                "Menter SST linearization; NASA TMR SST",
                "doc/source/numerical/sst_diagonal_w.rst",
            ),
            csv_relpath="tools/verification/data/sst_diagonal_w_theory.csv",
            plot_module="plot_sst_diagonal_w",
            pure_apis=("SstDiagonalW",),
            source_relpaths=(
                "codes/turb/include/TurbSstPolicy.h",
                "codes/turb/src/TurbSstPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/sst_diagonal_w.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "dw = -2*fBeta*beta*omega - absCD/(rho*omega);\n"
            ),
            proves="omega diagonal algebra vs pure API; cross-diff sample.",
            does_not_prove="Full residual Jacobian on a case mesh.",
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
            csv_relpath="tools/verification/data/sst_diagonal_k_theory.csv",
            plot_module="plot_sst_diagonal_k",
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
                "xk = (grad k * grad omega) / omega^3\n"
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
            csv_relpath="tools/verification/data/ns_lusgs_flux_increment_theory.csv",
            plot_module="plot_ns_lusgs_flux_increment",
            pure_apis=(
                "LusgsSplitFaceEigenvalues",
                "LusgsLinearizationWeights",
                "LusgsFluxIncrementFiveEq",
            ),
            source_relpaths=(
                "codes/ns/include/NsLusgsPolicy.h",
                "codes/ns/src/NsLusgsPolicy.cpp",
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
            id="viscous_energy_flux",
            family="mean-flow",
            genus="NS-viscous",
            section="energy / viscous flux",
            title="Viscous energy flux (stress work + heat)",
            role_tags=("viscous", "energy", "flux"),
            formula_plain=(
                "E_vis = u*f_mx + v*f_my + w*f_mz + q_n\n"
                "     = u · f_m + q_n"
            ),
            formula_latex=(
                r"E_{\mathrm{vis}}=\mathbf{u}\cdot\mathbf{f}_m+q_n"
            ),
            principle=(
                "Viscous energy face flux is the power of viscous stresses "
                "plus conductive heat flux through the face."
            ),
            background=(
                "Pure API ViscousEnergyFlux in NsVisFluxPolicy."
            ),
            references=(
                "Compressible NS energy equation viscous terms",
                "doc/source/numerical/viscous_energy_flux.rst",
            ),
            csv_relpath="tools/verification/data/viscous_energy_flux_theory.csv",
            plot_module="plot_viscous_energy_flux",
            pure_apis=("ViscousEnergyFlux",),
            source_relpaths=(
                "codes/ns/include/NsVisFluxPolicy.h",
                "codes/ns/src/NsVisFluxPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/viscous_energy_flux.rst",
            code_sketch="// GATED SKETCH\nE = u*fmx + v*fmy + w*fmz + qn;\n",
            proves="Energy-flux algebra vs pure API.",
            does_not_prove="Full viscous residual on a case mesh.",
        ),


        TheoryEntry(
            id="viscous_face_flux_coefficient",
            family="mean-flow",
            genus="NS-viscous",
            section="flux scaling / Reynolds",
            title="Viscous face flux coefficient (-A/Re)",
            role_tags=("viscous", "flux", "Reynolds"),
            formula_plain=(
                "coef = - invReynolds * faceArea\n"
                "Scales Stokes/heat face fluxes into residual units"
            ),
            formula_latex=(
                r"\mathrm{coef}=-\frac{A}{\mathrm{Re}}"
            ),
            principle=(
                "Non-dimensional NS places 1/Re on viscous fluxes. The face "
                "coefficient multiplies stress/heat contributions before "
                "deposition to cells."
            ),
            background=(
                "Pure API ViscousFaceFluxCoefficient in NsVisFluxPolicy."
            ),
            references=(
                "Non-dimensional compressible NS; OneFLOW NsVisFluxPolicy",
                "doc/source/numerical/viscous_face_flux_coefficient.rst",
            ),
            csv_relpath="tools/verification/data/viscous_face_flux_coefficient_theory.csv",
            plot_module="plot_viscous_face_flux_coefficient",
            pure_apis=("ViscousFaceFluxCoefficient",),
            source_relpaths=(
                "codes/ns/include/NsVisFluxPolicy.h",
                "codes/ns/src/NsVisFluxPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/viscous_face_flux_coefficient.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "coef = - invReynolds * faceArea;\n"
            ),
            proves="coef=-A/Re algebra vs pure API on an area grid.",
            does_not_prove="Full viscous residual magnitude on a mesh.",
        ),

        TheoryEntry(
            id="conductive_heat_flux_normal",
            family="mean-flow",
            genus="NS-viscous",
            section="heat flux / conduction",
            title="Conductive heat flux normal q_n = k * dT/dn",
            role_tags=("viscous", "heat-flux"),
            formula_plain=(
                "q_n = k_eff * (dT/dn)\n"
                "Feeds ViscousEnergyFlux as the heat term"
            ),
            formula_latex=(
                r"q_n=k\,\partial_n T"
            ),
            principle=(
                "Fourier conduction through a face: heat flux is conductivity "
                "times the normal temperature gradient."
            ),
            background=(
                "Pure API ConductiveHeatFluxNormal in NsVisFluxPolicy. "
                "Upstream: EffectiveThermalConductivity. Downstream: "
                "ViscousEnergyFlux."
            ),
            references=(
                "Fourier heat conduction; compressible NS energy equation",
                "doc/source/numerical/conductive_heat_flux_normal.rst",
            ),
            csv_relpath="tools/verification/data/conductive_heat_flux_normal_theory.csv",
            plot_module="plot_conductive_heat_flux_normal",
            pure_apis=("ConductiveHeatFluxNormal",),
            source_relpaths=(
                "codes/ns/include/NsVisFluxPolicy.h",
                "codes/ns/src/NsVisFluxPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/conductive_heat_flux_normal.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "qn = k * dTdn;\n"
            ),
            proves="q_n algebra vs pure API on a gradient grid.",
            does_not_prove="Wall heat-transfer coefficients on a case mesh.",
        ),

        TheoryEntry(
            id="project_stress_momentum_flux",
            family="mean-flow",
            genus="NS-viscous",
            section="stress / face projection",
            title="Project Stokes stress onto face normal (momentum flux)",
            role_tags=("viscous", "stress", "flux"),
            formula_plain=(
                "f_mx = nx*txx + ny*txy + nz*txz\n"
                "f_my = nx*txy + ny*tyy + nz*tyz\n"
                "f_mz = nx*txz + ny*tyz + nz*tzz"
            ),
            formula_latex=(
                r"\mathbf{f}_m=\boldsymbol{\tau}\cdot\mathbf{n}"
            ),
            principle=(
                "Maps the Stokes stress tensor to the viscous momentum flux "
                "through a face by contraction with the face normal."
            ),
            background=(
                "Pure API ProjectStressMomentumFlux in NsVisFluxPolicy. "
                "Upstream: StokesStressTensor. Downstream: ViscousEnergyFlux."
            ),
            references=(
                "Standard FV viscous flux assembly; OneFLOW NsVisFluxPolicy",
                "doc/source/numerical/project_stress_momentum_flux.rst",
            ),
            csv_relpath="tools/verification/data/project_stress_momentum_flux_theory.csv",
            plot_module="plot_project_stress_momentum_flux",
            pure_apis=("ProjectStressMomentumFlux",),
            source_relpaths=(
                "codes/ns/include/NsVisFluxPolicy.h",
                "codes/ns/src/NsVisFluxPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/project_stress_momentum_flux.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "fmx = nx*txx + ny*txy + nz*txz;\n"
            ),
            proves="tau·n algebra vs pure API on an nx grid.",
            does_not_prove="Full viscous residual assembly on a mesh.",
        ),

        TheoryEntry(
            id="stokes_stress_pure_shear",
            family="mean-flow",
            genus="NS-viscous",
            section="stress / Stokes tensor",
            title="Stokes stress tensor (pure shear and bulk sample)",
            role_tags=("viscous", "stress"),
            formula_plain=(
                "div*2/3 = (2/3)*(du/dx+dv/dy+dw/dz)\n"
                "txx = mu*(2*du/dx - div*2/3), ...\n"
                "txy = mu*(du/dy + dv/dx), ...\n"
                "Pure shear du/dy=1: txy=mu"
            ),
            formula_latex=(
                r"\tau_{ij}=\mu\left(\partial_i u_j+\partial_j u_i"
                r"-\frac{2}{3}\delta_{ij}\nabla\cdot\mathbf{u}\right)"
            ),
            principle=(
                "Newtonian Stokes stress for NS viscous fluxes. Pure shear "
                "isolates off-diagonal viscosity; bulk strain exercises the "
                "2/3 divergence correction."
            ),
            background=(
                "Pure API StokesStressTensor in NsVisFluxPolicy. Downstream: "
                "ProjectStressMomentumFlux, ViscousEnergyFlux."
            ),
            references=(
                "Standard Newtonian fluid constitutive law; OneFLOW NsVisFluxPolicy",
                "doc/source/numerical/stokes_stress_pure_shear.rst",
            ),
            csv_relpath="tools/verification/data/stokes_stress_pure_shear_theory.csv",
            plot_module="plot_stokes_stress_pure_shear",
            pure_apis=("StokesStressTensor",),
            source_relpaths=(
                "codes/ns/include/NsVisFluxPolicy.h",
                "codes/ns/src/NsVisFluxPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/stokes_stress_pure_shear.rst",
            code_sketch=(
                "// GATED SKETCH pure shear\n"
                "txy = mu * (dudy + dvdx);\n"
            ),
            proves="Pure-shear txy=mu and bulk 4/3 sample vs pure API.",
            does_not_prove="Full viscous face residual on a mesh.",
        ),

        TheoryEntry(
            id="effective_thermal_conductivity",
            family="mean-flow",
            genus="NS-viscous",
            section="heat flux / conductivity",
            title="Effective thermal conductivity (laminar + turbulent)",
            role_tags=("viscous", "heat-flux", "conductivity"),
            formula_plain=(
                "k_eff = (mu/Pr_l + mut/Pr_t) * cp\n"
                "     = (mu*invPr_l + mut*invPr_t) * cp"
            ),
            formula_latex=(
                r"k_{\mathrm{eff}}=\left(\frac{\mu}{\mathrm{Pr}_l}+\frac{\mu_t}{\mathrm{Pr}_t}\right)c_p"
            ),
            principle=(
                "Eddy-viscosity heat transport: laminar and turbulent "
                "contributions add after division by Prandtl numbers, "
                "then scale by specific heat to conductivity."
            ),
            background=(
                "Pure API EffectiveThermalConductivity in NsVisFluxPolicy. "
                "Downstream: ConductiveHeatFluxNormal."
            ),
            references=(
                "Standard RANS heat-flux closure; OneFLOW NsVisFluxPolicy",
                "doc/source/numerical/effective_thermal_conductivity.rst",
            ),
            csv_relpath="tools/verification/data/effective_thermal_conductivity_theory.csv",
            plot_module="plot_effective_thermal_conductivity",
            pure_apis=("EffectiveThermalConductivity",),
            source_relpaths=(
                "codes/ns/include/NsVisFluxPolicy.h",
                "codes/ns/src/NsVisFluxPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/effective_thermal_conductivity.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "k = (mu*invPrL + mut*invPrT) * cp;\n"
            ),
            proves="k_eff algebra vs pure API on a mut grid.",
            does_not_prove="Full viscous face heat flux residual on a mesh.",
        ),

        TheoryEntry(
            id="combine_inv_vis_timestep",
            family="mean-flow",
            genus="NS-spectrum",
            section="time-step / inviscid-viscous combine",
            title="Harmonic combine of inviscid and viscous time steps",
            role_tags=("time-step", "CFL"),
            formula_plain=(
                "dt = dt_inv * dt_vis / (dt_inv + dt_vis)\n"
                "Harmonic-type bound: limited by the stiffer (smaller) side"
            ),
            formula_latex=(
                r"\Delta t=\frac{\Delta t_{\mathrm{inv}}\Delta t_{\mathrm{vis}}}"
                r"{\Delta t_{\mathrm{inv}}+\Delta t_{\mathrm{vis}}}"
            ),
            principle=(
                "After separate inviscid and viscous Delta-t estimates, combine "
                "so either stiff channel reduces the final step. Equivalent to "
                "parallel resistors / harmonic mean scaled by 1/2 when equal."
            ),
            background=(
                "Pure API CombineInvAndVisTimeStep in TimeStepPolicy. Upstream: "
                "CellTimeStepFromSpectrum on inviscid and viscous spectra."
            ),
            references=(
                "Standard CFD local time-step practice; OneFLOW TimeStepPolicy",
                "doc/source/numerical/combine_inv_vis_timestep.rst",
            ),
            csv_relpath="tools/verification/data/combine_inv_vis_timestep_theory.csv",
            plot_module="plot_combine_inv_vis_timestep",
            pure_apis=("CombineInvAndVisTimeStep",),
            source_relpaths=(
                "codes/ns/include/TimeStepPolicy.h",
                "codes/ns/src/TimeStepPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/combine_inv_vis_timestep.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "dt = dtInv * dtVis / (dtInv + dtVis);\n"
            ),
            proves="Harmonic combine algebra vs pure API on a dt_inv grid.",
            does_not_prove="Optimal CFL on a case mesh or spectrum accuracy.",
        ),

        TheoryEntry(
            id="viscous_spectral_radius_model3",
            family="mean-flow",
            genus="NS-spectrum",
            section="time-step / viscous spectral radius",
            title="Viscous face spectral radius (model 3, minimal)",
            role_tags=("time-step", "spectrum", "viscous"),
            formula_plain=(
                "lambda = A^2 / (Re * rho)\n"
                "No explicit mu; cheapest viscous spectral placeholder"
            ),
            formula_latex=(
                r"\lambda_{\mathrm{vis}}=\frac{A^{2}}{\mathrm{Re}\,\rho}"
            ),
            principle=(
                "Minimal A^2 diffusion-scale bound with viscous magnitude "
                "absorbed into Re. Fallback when full mu/Pr fields are unused."
            ),
            background=(
                "Pure API ViscousFaceSpectralRadiusModel3. Family overview: "
                "doc/source/numerical/viscous_spectral_radius.rst"
            ),
            references=(
                "Blazek time-step chapters (motivation); OneFLOW TimeStepPolicy",
                "doc/source/numerical/viscous_spectral_radius_model3.rst",
            ),
            csv_relpath="tools/verification/data/viscous_spectral_radius_model3_theory.csv",
            plot_module="plot_viscous_spectral_radius_model3",
            pure_apis=("ViscousFaceSpectralRadiusModel3",),
            source_relpaths=(
                "codes/ns/include/TimeStepPolicy.h",
                "codes/ns/src/TimeStepPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/viscous_spectral_radius_model3.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "lambda = A*A / (Re * rho + SMALL);\n"
            ),
            proves="Model-3 algebra and A^2 scaling vs pure API.",
            does_not_prove="Adequacy vs model 1 in high-mut boundary layers.",
        ),

        TheoryEntry(
            id="viscous_spectral_radius_model2",
            family="mean-flow",
            genus="NS-spectrum",
            section="time-step / viscous spectral radius",
            title="Viscous face spectral radius (model 2, length-scale)",
            role_tags=("time-step", "spectrum", "viscous"),
            formula_plain=(
                "c1 = 2 * mu_tot / (rho * d * Re)\n"
                "lambda = 0.5 * c1 * A\n"
                "d = cell/wall distance scale"
            ),
            formula_latex=(
                r"c_1=\frac{2\mu_{\mathrm{tot}}}{\rho\,d\,\mathrm{Re}},\quad "
                r"\lambda_{\mathrm{vis}}=\frac{1}{2} c_1 A"
            ),
            principle=(
                "Uses an explicit length scale d (cell or wall distance) instead "
                "of area-squared only. Diffusion scale nu/d^2 mapped to a face "
                "spectral contribution proportional to A."
            ),
            background=(
                "Pure API ViscousFaceSpectralRadiusModel2 in TimeStepPolicy. "
                "Companion to model 1 (Blazek max diffusivity) and model 3 "
                "(A^2/(Re rho) only)."
            ),
            references=(
                "Blazek CFD; OneFLOW TimeStepPolicy",
                "doc/source/numerical/viscous_spectral_radius_model2.rst",
            ),
            csv_relpath="tools/verification/data/viscous_spectral_radius_model2_theory.csv",
            plot_module="plot_viscous_spectral_radius_model2",
            pure_apis=("ViscousFaceSpectralRadiusModel2",),
            source_relpaths=(
                "codes/ns/include/TimeStepPolicy.h",
                "codes/ns/src/TimeStepPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/viscous_spectral_radius_model2.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "c1 = 2*mu / (rho*d*Re + SMALL);\n"
                "lambda = 0.5 * c1 * A;\n"
            ),
            proves="Model-2 algebra vs pure API on a viscosity grid; 1/d scaling.",
            does_not_prove="Optimal choice of d on a mesh or model 1/3 equivalence.",
        ),

        TheoryEntry(
            id="viscous_spectral_radius_model1",
            family="mean-flow",
            genus="NS-spectrum",
            section="time-step / viscous spectral radius",
            title="Viscous face spectral radius (model 1, Blazek-style)",
            role_tags=("time-step", "spectrum", "viscous"),
            formula_plain=(
                "c1 = (4/3)*(mu + mut)\n"
                "c2 = gamma*(mu/Pr_l + mut/Pr_t)\n"
                "c3 = 2*max(c1,c2)/(Re*rho)\n"
                "lambda = A^2 * c3"
            ),
            formula_latex=(
                r"c_1=\frac{4}{3}(\mu+\mu_t),\ "
                r"c_2=\gamma(\mu/\mathrm{Pr}_l+\mu_t/\mathrm{Pr}_t)"
                "\n"
                r"\lambda_{\mathrm{vis}}=A^{2}\cdot\frac{2\max(c_1,c_2)}{\mathrm{Re}\,\rho}"
            ),
                        principle=(
                "Diffusion scale nu*k^2 with k~A/V yields lambda_vis ~ A^2 * nu_eff. "
                "Model 1 takes nu_eff from max of momentum (4/3 viscosity) and energy "
                "(gamma*mu/Pr) channels; engineering time-step bound, not exact eigenradius."
            ),
            background=(
                "Pure API ViscousFaceSpectralRadiusModel1 in TimeStepPolicy. "
                "Derivation: doc/source/numerical/viscous_spectral_radius_model1.rst "
                "(Derivation section). Feeds viscous dt and related spectra."
            ),
references=(
                "Blazek, Computational Fluid Dynamics; OneFLOW TimeStepPolicy",
                "doc/source/numerical/viscous_spectral_radius_model1.rst",
            ),
            csv_relpath="tools/verification/data/viscous_spectral_radius_model1_theory.csv",
            plot_module="plot_viscous_spectral_radius_model1",
            pure_apis=("ViscousFaceSpectralRadiusModel1",),
            source_relpaths=(
                "codes/ns/include/TimeStepPolicy.h",
                "codes/ns/src/TimeStepPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/viscous_spectral_radius_model1.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "c1 = (4.0/3.0)*(mu+mut);\n"
                "c2 = gama*(mu*invPrL + mut*invPrT);\n"
                "lambda = A*A * 2*max(c1,c2)/(Re*rho);\n"
            ),
            proves="Model-1 viscous spectral algebra vs pure API on a mut grid.",
            does_not_prove="Full CFL time-step on a case mesh or models 2/3.",
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
            id="euler_convective_flux",
            family="mean-flow",
            genus="NS-inviscid",
            section="face flux / Euler convective",
            title="Euler convective face flux",
            role_tags=("inviscid", "flux", "Euler"),
            formula_plain=(
                "mass = rho * vn\n"
                "f_rho = mass\n"
                "f_m = mass * u + n * p\n"
                "f_E = mass * h + vg * p"
            ),
            formula_latex=(
                r"\mathbf{F}^c=\rho v_n\begin{bmatrix}1\\ \mathbf{u}\\ H\end{bmatrix}"
                r"+p\begin{bmatrix}0\\ \mathbf{n}\\ v_g\end{bmatrix}"
            ),
            principle=(
                "Inviscid face flux is mass flux times transported quantities "
                "plus pressure acting on the face; energy includes grid-speed "
                "pressure work for ALE consistency."
            ),
            background=(
                "Pure API EulerConvectiveFlux in NsInvFluxPolicy. Building "
                "block for Roe/AUSM residual assembly."
            ),
            references=(
                "Blazek CFD; standard FV Euler flux",
                "doc/source/numerical/euler_convective_flux.rst",
            ),
            csv_relpath="tools/verification/data/euler_convective_flux_theory.csv",
            plot_module="plot_euler_convective_flux",
            pure_apis=("EulerConvectiveFlux",),
            source_relpaths=(
                "codes/ns/include/NsInvFluxPolicy.h",
                "codes/ns/src/NsInvFluxPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/euler_convective_flux.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "mass = rho * vn;\n"
                "f_rho = mass; f_m = mass*u + n*p; f_E = mass*h + vg*p;\n"
            ),
            proves="Euler face flux algebra vs pure API on a vn grid.",
            does_not_prove="Full residual on a mesh or Riemann solver quality.",
        ),

        TheoryEntry(
            id="face_mach_number",
            family="mean-flow",
            genus="NS-inviscid",
            section="face kinematics / Mach",
            title="Face Mach number M = vn / c",
            role_tags=("inviscid", "Mach", "AUSM"),
            formula_plain=(
                "M = vn_rel / soundSpeed\n"
                "Used by AUSM-type mass/pressure splits at faces"
            ),
            formula_latex=(
                r"M=\frac{v_n}{c}"
            ),
            principle=(
                "Relative normal Mach number at a face. Sign indicates inflow "
                "or outflow relative to the face normal orientation."
            ),
            background=(
                "Pure API FaceMachNumber in NsInvFluxPolicy. Consumers include "
                "AUSM face reconstructions in NsInvFlux."
            ),
            references=(
                "AUSM family (Liou); OneFLOW NsInvFluxPolicy",
                "doc/source/numerical/face_mach_number.rst",
            ),
            csv_relpath="tools/verification/data/face_mach_number_theory.csv",
            plot_module="plot_face_mach_number",
            pure_apis=("FaceMachNumber",),
            source_relpaths=(
                "codes/ns/include/NsInvFluxPolicy.h",
                "codes/ns/src/NsInvFluxPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/face_mach_number.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "M = vnRel / soundSpeed;\n"
            ),
            proves="M=vn/c algebra vs pure API on a vn grid.",
            does_not_prove="Full AUSM flux residual or entropy-fix interaction.",
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
                "// GATED SKETCH - production is inline in NsInvFluxPolicy.h\n"
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
            id="roe_pressure_from_enthalpy",
            family="mean-flow",
            genus="NS-EOS",
            section="thermodynamics / face pressure",
            title="Roe pressure from total enthalpy (LUSGS mode 0)",
            role_tags=("EOS", "pressure", "LUSGS"),
            formula_plain=(
                "p = rho * (H - 0.5*|u|^2) * (gamma - 1) / gamma\n"
                "LusgsFacePressureFromMode(0,...) uses this closure"
            ),
            formula_latex=(
                r"p=\rho\left(H-\frac{1}{2}|u|^{2}\right)\frac{\gamma-1}{\gamma}"
            ),
            principle=(
                "Recover static pressure from Roe-averaged total enthalpy by "
                "removing kinetic energy and applying the ideal-gas factor."
            ),
            background=(
                "Pure API RoePressureFromTotalEnthalpy. Selected by "
                "LusgsFacePressureFromMode when mode==0."
            ),
            references=(
                "Roe average state; OneFLOW NsInvFluxPolicy",
                "doc/source/numerical/roe_pressure_from_enthalpy.rst",
            ),
            csv_relpath="tools/verification/data/roe_pressure_from_enthalpy_theory.csv",
            plot_module="plot_roe_pressure_from_enthalpy",
            pure_apis=("RoePressureFromTotalEnthalpy", "LusgsFacePressureFromMode"),
            source_relpaths=(
                "codes/ns/include/NsInvFluxPolicy.h",
                "codes/ns/src/NsInvFluxPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/roe_pressure_from_enthalpy.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "p = rho * (H - 0.5 * speedSquared) * (gama - 1) / gama;\n"
            ),
            proves="Pressure-from-enthalpy algebra vs pure API on an H grid.",
            does_not_prove="Full LUSGS residual or mode-1 Roe-weighted pressure path.",
        ),


        TheoryEntry(
            id="ideal_gas_sound_speed_from_pressure",
            family="mean-flow",
            genus="NS-EOS",
            section="thermodynamics / ideal gas",
            title="Ideal-gas sound speed from pressure c(p, rho, gamma)",
            role_tags=("EOS", "sound-speed"),
            formula_plain=(
                "c = sqrt( gamma * p / rho )\n"
                "Pressure-density form (vs total-enthalpy form)"
            ),
            formula_latex=(
                r"c=\sqrt{\gamma\frac{p}{\rho}}"
            ),
            principle=(
                "Isentropic sound speed from local static pressure and density. "
                "Companion to IdealGasSoundSpeed (total-enthalpy form)."
            ),
            background=(
                "Pure API IdealGasSoundSpeedFromPressure in NsInvFluxPolicy. "
                "Related: IdealGasSoundSpeed, RoeAbsoluteEigenvalues."
            ),
            references=(
                "Anderson, Modern Compressible Flow; OneFLOW NsInvFluxPolicy",
                "doc/source/numerical/ideal_gas_sound_speed_from_pressure.rst",
            ),
            csv_relpath="tools/verification/data/ideal_gas_sound_speed_from_pressure_theory.csv",
            plot_module="plot_ideal_gas_sound_speed_from_pressure",
            pure_apis=("IdealGasSoundSpeedFromPressure",),
            source_relpaths=(
                "codes/ns/include/NsInvFluxPolicy.h",
                "codes/ns/src/NsInvFluxPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/ideal_gas_sound_speed_from_pressure.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "c = sqrt(fabs(gama * pressure / density));\n"
            ),
            proves="c=sqrt(gamma*p/rho) vs pure API on a pressure grid.",
            does_not_prove="Consistency with total-enthalpy form on full Roe state.",
        ),

        TheoryEntry(
            id="ideal_gas_sound_speed",
            family="mean-flow",
            genus="NS-EOS",
            section="thermodynamics / ideal gas",
            title="Ideal-gas sound speed from total enthalpy",
            role_tags=("EOS", "sound-speed"),
            formula_plain=(
                "c = sqrt( (gamma-1) * (H - 0.5*|u|^2) )\n"
                "With |u|=0: c = sqrt( (gamma-1)*H )"
            ),
            formula_latex=(
                r"c=\sqrt{(\gamma-1)\bigl(H-\frac{1}{2}|u|^{2}\bigr)}"
            ),
            principle=(
                "Isentropic sound speed from total enthalpy after removing "
                "kinetic energy. Used by Roe eigenvalues and time-step spectra."
            ),
            background=(
                "Pure API IdealGasSoundSpeed in NsInvFluxPolicy. "
                "Upstream: IdealGasStaticEnthalpy / TotalEnthalpy."
            ),
            references=(
                "Anderson, Modern Compressible Flow; OneFLOW NsInvFluxPolicy",
                "doc/source/numerical/ideal_gas_sound_speed.rst",
            ),
            csv_relpath="tools/verification/data/ideal_gas_sound_speed_theory.csv",
            plot_module="plot_ideal_gas_sound_speed",
            pure_apis=("IdealGasSoundSpeed",),
            source_relpaths=(
                "codes/ns/include/NsInvFluxPolicy.h",
                "codes/ns/src/NsInvFluxPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/ideal_gas_sound_speed.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "c2 = (gama - 1) * (H - 0.5 * speedSquared);\n"
                "c = sqrt(fabs(c2));\n"
            ),
            proves="Sound-speed algebra vs pure API on an H grid (and kinetic sample).",
            does_not_prove="Variable-gamma gas or numerical flux case residuals.",
        ),





        TheoryEntry(
            id="ideal_gas_conservative_to_prim",
            family="mean-flow",
            genus="NS-EOS",
            section="thermodynamics / ideal gas",
            title="Ideal-gas conservative to prim map",
            role_tags=("EOS", "primitive"),
            formula_plain=(
                "rho = |q_rho|\n"
                "u = q_mom / rho\n"
                "p = (gamma-1) * (q_rhoE - 0.5*rho*|u|^2)"
            ),
            formula_latex=(
                r"p=(\gamma-1)\bigl(\rho e-\tfrac12\rho|u|^{2}\bigr)"
            ),
            principle=(
                "Inverse of IdealGasPrimToConservative: recover primitives "
                "from conserved density, momentum, and total energy."
            ),
            background=(
                "Pure API IdealGasConservativeToPrim in NsInvFluxPolicy. "
                "Pair with IdealGasPrimToConservative for round-trip tests."
            ),
            references=(
                "Standard ideal-gas FV variable conversion",
                "doc/source/numerical/ideal_gas_conservative_to_prim.rst",
            ),
            csv_relpath="tools/verification/data/ideal_gas_conservative_to_prim_theory.csv",
            plot_module="plot_ideal_gas_conservative_to_prim",
            pure_apis=("IdealGasConservativeToPrim",),
            source_relpaths=(
                "codes/ns/include/NsInvFluxPolicy.h",
                "codes/ns/src/NsInvFluxPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/ideal_gas_conservative_to_prim.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "p = (gama-1)*(qRhoE - 0.5*rho*|u|^2);\n"
            ),
            proves="Cons-to-prim algebra on energy grid; kinetic subtraction.",
            does_not_prove="Robustness at near-vacuum or negative pressure.",
        ),

        TheoryEntry(
            id="ideal_gas_prim_to_conservative",
            family="mean-flow",
            genus="NS-EOS",
            section="thermodynamics / ideal gas",
            title="Ideal-gas prim to conservative map",
            role_tags=("EOS", "conservative"),
            formula_plain=(
                "q_rho = rho\n"
                "q_mom = rho * u\n"
                "q_rhoE = rho * ( p/(rho*(gamma-1)) + |u|^2/2 )"
            ),
            formula_latex=(
                r"\mathbf{q}=\begin{bmatrix}\rho\\ \rho\mathbf{u}\\ "
                r"\rho\bigl(\frac{p}{\rho(\gamma-1)}+\frac12|u|^{2}\bigr)"
                r"\end{bmatrix}"
            ),
            principle=(
                "Maps primitive thermodynamic/kinematic state to conserved "
                "variables for the ideal-gas five-equation system."
            ),
            background=(
                "Pure API IdealGasPrimToConservative; uses "
                "IdealGasSpecificInternalEnergy. Inverse: "
                "IdealGasConservativeToPrim."
            ),
            references=(
                "Standard ideal-gas FV variable conversion",
                "doc/source/numerical/ideal_gas_prim_to_conservative.rst",
            ),
            csv_relpath="tools/verification/data/ideal_gas_prim_to_conservative_theory.csv",
            plot_module="plot_ideal_gas_prim_to_conservative",
            pure_apis=("IdealGasPrimToConservative",),
            source_relpaths=(
                "codes/ns/include/NsInvFluxPolicy.h",
                "codes/ns/src/NsInvFluxPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/ideal_gas_prim_to_conservative.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "em = p/(rho*(g-1)) + 0.5*|u|^2;\n"
                "qRhoE = rho * em;\n"
            ),
            proves="Prim-to-cons algebra on pressure and kinetic samples.",
            does_not_prove="Variable-gamma gas or multi-species mixture.",
        ),

        TheoryEntry(
            id="ideal_gas_specific_internal_energy",
            family="mean-flow",
            genus="NS-EOS",
            section="thermodynamics / ideal gas",
            title="Ideal-gas specific energy e(p, rho, |u|^2, gamma)",
            role_tags=("EOS", "energy"),
            formula_plain=(
                "e = p / (rho * (gamma - 1)) + 0.5 * |u|^2\n"
                "Feeds conservative energy and total enthalpy paths"
            ),
            formula_latex=(
                r"e=\frac{p}{\rho(\gamma-1)}+\frac{1}{2}|u|^{2}"
            ),
            principle=(
                "Calorically perfect gas specific total energy: thermal part "
                "from equation of state plus kinetic energy."
            ),
            background=(
                "Pure API IdealGasSpecificInternalEnergy in NsInvFluxPolicy. "
                "Related: IdealGasStaticEnthalpy, TotalEnthalpy, "
                "IdealGasPrimToConservative."
            ),
            references=(
                "Anderson, Modern Compressible Flow; OneFLOW NsInvFluxPolicy",
                "doc/source/numerical/ideal_gas_specific_internal_energy.rst",
            ),
            csv_relpath="tools/verification/data/ideal_gas_specific_internal_energy_theory.csv",
            plot_module="plot_ideal_gas_specific_internal_energy",
            pure_apis=("IdealGasSpecificInternalEnergy",),
            source_relpaths=(
                "codes/ns/include/NsInvFluxPolicy.h",
                "codes/ns/src/NsInvFluxPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/ideal_gas_specific_internal_energy.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "e = (p/rho)/(gama-1) + 0.5*speedSquared;\n"
            ),
            proves="Specific-energy algebra vs pure API on a kinetic grid.",
            does_not_prove="Variable-gamma gas or full conservative residual.",
        ),

        TheoryEntry(
            id="total_enthalpy",
            family="mean-flow",
            genus="NS-EOS",
            section="thermodynamics / ideal gas",
            title="Total enthalpy H = h + |u|^2/2",
            role_tags=("EOS", "enthalpy"),
            formula_plain=(
                "H = h + 0.5 * |u|^2\n"
                "Feeds Roe average, sound speed, pressure-from-H"
            ),
            formula_latex=(
                r"H=h+\frac{1}{2}|u|^{2}"
            ),
            principle=(
                "Total enthalpy adds kinetic energy to static enthalpy. "
                "Central link between IdealGasStaticEnthalpy and Roe/LUSGS "
                "face thermodynamics."
            ),
            background=(
                "Pure API TotalEnthalpy in NsInvFluxPolicy. Upstream: static "
                "enthalpy; downstream: IdealGasSoundSpeed, "
                "RoePressureFromTotalEnthalpy."
            ),
            references=(
                "Anderson, Modern Compressible Flow; OneFLOW NsInvFluxPolicy",
                "doc/source/numerical/total_enthalpy.rst",
            ),
            csv_relpath="tools/verification/data/total_enthalpy_theory.csv",
            plot_module="plot_total_enthalpy",
            pure_apis=("TotalEnthalpy",),
            source_relpaths=(
                "codes/ns/include/NsInvFluxPolicy.h",
                "codes/ns/src/NsInvFluxPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/total_enthalpy.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "H = h + 0.5 * speedSquared;\n"
            ),
            proves="Total-enthalpy algebra vs pure API on a kinetic grid.",
            does_not_prove="Full Roe average state consistency on a mesh.",
        ),

        TheoryEntry(
            id="ideal_gas_static_enthalpy",
            family="mean-flow",
            genus="NS-EOS",
            section="thermodynamics / ideal gas",
            title="Ideal-gas static enthalpy h(p, rho, gamma)",
            role_tags=("EOS", "enthalpy"),
            formula_plain=(
                "h = (gamma / (gamma - 1)) * (p / rho)\n"
                "gamma = 1.4 -> slope 3.5 on the p/rho axis"
            ),
            formula_latex=(
                r"h=\frac{\gamma}{\gamma-1}\frac{p}{\rho}"
            ),
            principle=(
                "Calorically perfect gas static enthalpy from pressure and density. "
                "Feeds Roe average total enthalpy and sound-speed closures."
            ),
            background=(
                "Pure API IdealGasStaticEnthalpy in NsInvFluxPolicy. "
                "Related: IdealGasSpecificInternalEnergy, IdealGasSoundSpeed."
            ),
            references=(
                "Anderson, Modern Compressible Flow; OneFLOW NsInvFluxPolicy",
                "doc/source/numerical/ideal_gas_static_enthalpy.rst",
            ),
            csv_relpath="tools/verification/data/ideal_gas_static_enthalpy_theory.csv",
            plot_module="plot_ideal_gas_static_enthalpy",
            pure_apis=("IdealGasStaticEnthalpy",),
            source_relpaths=(
                "codes/ns/include/NsInvFluxPolicy.h",
                "codes/ns/src/NsInvFluxPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/ideal_gas_static_enthalpy.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "h = (gama / (gama - 1.0)) * (pressure / density);\n"
            ),
            proves="Static enthalpy algebra vs pure API on a p/rho grid.",
            does_not_prove="Variable-gamma gas or real-gas thermodynamics.",
        ),




        TheoryEntry(
            id="roe_total_enthalpy_jump",
            family="mean-flow",
            genus="NS-inviscid",
            section="Roe dissipation / enthalpy jump",
            title="Roe total-enthalpy jump",
            role_tags=("inviscid", "Roe", "dissipation"),
            formula_plain=(
                "ae = gamma-1\n"
                "af = 0.5*ae*|u|^2\n"
                "dh = -ae*(u·dq_mom - dqE) + af*dqRho"
            ),
            formula_latex=(
                r"dh=-(\gamma-1)\bigl(\mathbf{u}\cdot\Delta(\rho\mathbf{u})"
                r"-\Delta(\rho e)\bigr)"
                r"+\tfrac12(\gamma-1)|u|^{2}\,\Delta\rho"
            ),
            principle=(
                "Enthalpy-related jump used in Roe dissipation coupling. "
                "Rest mid-state reduces to dh=(gamma-1)*Delta(rho e)."
            ),
            background=(
                "Pure API RoeTotalEnthalpyJump in NsInvFluxPolicy. "
                "Feeds RoeDissipationCouplingTerms."
            ),
            references=(
                "Roe approximate Riemann solver; Toro Riemann solvers",
                "doc/source/numerical/roe_total_enthalpy_jump.rst",
            ),
            csv_relpath="tools/verification/data/roe_total_enthalpy_jump_theory.csv",
            plot_module="plot_roe_total_enthalpy_jump",
            pure_apis=("RoeTotalEnthalpyJump",),
            source_relpaths=(
                "codes/ns/include/NsInvFluxPolicy.h",
                "codes/ns/src/NsInvFluxPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/roe_total_enthalpy_jump.rst",
            code_sketch=(
                "// GATED SKETCH rest mid-state\n"
                "dh = (gama-1)*dqE;\n"
            ),
            proves="Enthalpy-jump algebra vs pure API on energy and kinetic samples.",
            does_not_prove="Full Roe residual or entropy-fix behavior.",
        ),

        TheoryEntry(
            id="roe_continuity_normal_jump",
            family="mean-flow",
            genus="NS-inviscid",
            section="Roe dissipation / continuity jump",
            title="Roe continuity normal jump",
            role_tags=("inviscid", "Roe", "dissipation"),
            formula_plain=(
                "jump = vn * dqRho - nx*dqRhoU - ny*dqRhoV - nz*dqRhoW\n"
                "     = vn * Delta(rho) - n · Delta(rho u)"
            ),
            formula_latex=(
                r"J_{\rho}=v_n\Delta\rho-\mathbf{n}\cdot\Delta(\rho\mathbf{u})"
            ),
            principle=(
                "Continuity contribution to Roe dissipation uses the "
                "relative normal velocity times density jump minus the "
                "projected momentum jump."
            ),
            background=(
                "Pure API RoeContinuityNormalJump in NsInvFluxPolicy. "
                "Part of Roe wave-strength assembly."
            ),
            references=(
                "Roe approximate Riemann solver; Toro Riemann solvers",
                "doc/source/numerical/roe_continuity_normal_jump.rst",
            ),
            csv_relpath="tools/verification/data/roe_continuity_normal_jump_theory.csv",
            plot_module="plot_roe_continuity_normal_jump",
            pure_apis=("RoeContinuityNormalJump",),
            source_relpaths=(
                "codes/ns/include/NsInvFluxPolicy.h",
                "codes/ns/src/NsInvFluxPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/roe_continuity_normal_jump.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "J = vn*dqRho - nx*dqU - ny*dqV - nz*dqW;\n"
            ),
            proves="Continuity jump algebra vs pure API.",
            does_not_prove="Full Roe residual or entropy-fix interaction.",
        ),

        TheoryEntry(
            id="roe_average_component",
            family="mean-flow",
            genus="NS-inviscid",
            section="Roe average / component blend",
            title="Roe average component (weighted mean)",
            role_tags=("inviscid", "Roe", "average"),
            formula_plain=(
                "avg = (q_L + q_R * ratio) * coef\n"
                "ratio = sqrt(rho_R/rho_L), coef = 1/(1+ratio)"
            ),
            formula_latex=(
                r"\tilde{q}=(q_L+q_R\,r)\,c,\quad "
                r"r=\sqrt{\rho_R/\rho_L},\; c=\frac{1}{1+r}"
            ),
            principle=(
                "Each Roe-averaged primitive is a density-weighted blend of "
                "left and right states using the Roe square-root weights."
            ),
            background=(
                "Pure API RoeAverageComponent; weights from RoeAverageWeights. "
                "Assembled in RoeAverageState."
            ),
            references=(
                "P. L. Roe, approximate Riemann solvers; Toro Riemann solvers",
                "doc/source/numerical/roe_average_component.rst",
            ),
            csv_relpath="tools/verification/data/roe_average_component_theory.csv",
            plot_module="plot_roe_average_component",
            pure_apis=("RoeAverageComponent", "RoeAverageWeights"),
            source_relpaths=(
                "codes/ns/include/NsInvFluxPolicy.h",
                "codes/ns/src/NsInvFluxPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/roe_average_component.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "avg = (left + right * ratio) * coef;\n"
            ),
            proves="Component blend algebra vs pure API; weight consistency.",
            does_not_prove="Full Roe flux residual or entropy-fix behavior.",
        ),

        TheoryEntry(
            id="roe_average_weights",
            family="mean-flow",
            genus="NS-inviscid",
            section="flux / Roe average",
            title="Roe average weights (ratio, coef)",
            role_tags=("flux", "Roe", "average"),
            formula_plain=(
                "ratio = sqrt(rho_R / rho_L)\n"
                "coef  = 1 / (1 + ratio)\n"
                "Component average: (left + right*ratio)*coef"
            ),
            formula_latex=(
                r"\mathrm{ratio}=\sqrt{\rho_R/\rho_L}"
                "\n"
                r"\mathrm{coef}=\frac{1}{1+\mathrm{ratio}}"
            ),
            principle=(
                "Density-based weights for the Roe-averaged interface state. "
                "Equal densities give ratio=1, coef=1/2 (arithmetic mean)."
            ),
            background=(
                "Pure API RoeAverageWeights. Feeds RoeAverageComponent and "
                "RoeAverageState used by the Roe flux."
            ),
            references=(
                "Roe (1981); OneFLOW NsInvFluxPolicy",
                "doc/source/numerical/roe_average_weights.rst",
            ),
            csv_relpath="tools/verification/data/roe_average_weights_theory.csv",
            plot_module="plot_roe_average_weights",
            pure_apis=("RoeAverageWeights", "RoeAverageComponent"),
            source_relpaths=(
                "codes/ns/include/NsInvFluxPolicy.h",
                "codes/ns/src/NsInvFluxPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/roe_average_weights.rst",
            code_sketch=(
                "// GATED SKETCH\n"
                "ratio = sqrt(rhoR / rhoL);\n"
                "coef  = 1.0 / (1.0 + ratio);\n"
            ),
            proves="ratio/coef algebra vs pure API on a density-ratio grid.",
            does_not_prove="Full RoeAverageState pressure/enthalpy consistency.",
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
                        csv_relpath="tools/verification/data/harten_entropy_fix_theory.csv",
            plot_module="plot_harten_entropy_fix",
            pure_apis=("ApplyHartenEntropyFix",),
            source_relpaths=(
                "codes/ns/include/NsInvFluxPolicy.h",
                "codes/ns/src/NsInvFluxPolicy.cpp",
            ),
            doc_relpath="doc/source/numerical/harten_entropy_fix.rst",
            code_sketch=(
                "// GATED SKETCH - mirror production MAX in tests.\n"
                "lambda_max = max(|vn+c|, |vn-c|);\n"
                "m1 = lambda_max * centropy1; m2 = lambda_max * centropy2;\n"
                "|vn| = max(m1, |vn|); |vn+/-c| = max(m2, |vn+/-c|);\n"
            ),
proves="Entropy-fix pure policy tests.",
            does_not_prove="Full Roe flux shock-capturing quality.",
        ),
        TheoryEntry(
            id="time_integral_rk_stage_weight",
            family="numerical-method",
            genus="Time integration",
            section="explicit multi-stage control",
            title="OneFLOW multi-stage coefficient schedule",
            role_tags=("time-integration", "multi-stage", "solver-control"),
            formula_plain=(
                "N = rk_stage\n"
                "a_s = 1 / (s + 1),  s = 0, ..., N - 1\n"
                "ctrl.lhscoef = a_s  (on the finest grid)"
            ),
            formula_latex=(
                r"N=\texttt{rk\_stage}"
                "\n"
                r"a_s=\frac{1}{s+1},\quad s=0,\ldots,N-1"
                "\n"
                r"\texttt{lhscoef}=a_s"
            ),
            principle=(
                "The configured stage count determines the number of stored "
                "stage weights. TimeIntegral::RungeKutta assigns each weight to "
                "ctrl.lhscoef before executing that stage's residual and update tasks."
            ),
            background=(
                "This documents the coefficient schedule implemented by OneFLOW. "
                "It is not a Butcher tableau and does not establish the formal order, "
                "stability region, or accuracy of the complete time integrator."
            ),
            references=(
                "OneFLOW: codes/global/src/Ctrl.cpp and solver control Init implementations.",
                "OneFLOW: codes/solver/src/TimeIntegral.cpp.",
                "OneFLOW: codes/solver/src/EulerRungeKuttaCapability.cpp.",
                "doc/source/numerical/time_integral_rk_stage_weight.rst",
            ),
            csv_relpath="tools/verification/data/time_integral_rk_stage_weight_theory.csv",
            plot_module="plot_time_integral_rk_stage_weight",
            pure_apis=("RungeKuttaStageWeight",),
            source_relpaths=(
                "codes/global/include/TimeIntegralPolicy.h",
                "codes/global/src/Ctrl.cpp",
                "codes/ns/src/NsCtrl.cpp",
                "codes/ins/src/INsCtrl.cpp",
                "codes/turb/src/TurbCtrl.cpp",
                "codes/solver/src/TimeIntegral.cpp",
                "codes/solver/src/EulerRungeKuttaCapability.cpp",
                "codes/solver/include/EulerRungeKuttaCapability.h",
            ),
            doc_relpath="doc/source/numerical/time_integral_rk_stage_weight.rst",
            code_sketch=(
                "// GATED SKETCH - mirrors the current source, not a general RK tableau.\n"
                "for (int stage = 0; stage < rk_stage; ++stage) {\n"
                "    rk_coef.push_back(1.0 / (stage + 1));\n"
                "}\n"
                "for (int stage = 0; stage < rk_coef.size(); ++stage) {\n"
                "    ctrl.lhscoef = ctrl.rk_coef[stage];\n"
                "    execute_residual_and_update_tasks();\n"
                "}\n"
            ),
            proves=(
                "The source-level mapping from rk_stage to the stored coefficient "
                "sequence and finest-grid stage loop."
            ),
            does_not_prove=(
                "A complete Butcher tableau, formal temporal order, stability, "
                "convergence rate, or CFD accuracy on a physical case."
            ),
            examples=(
                "Input: rk_stage = 4.\n"
                "Schedule: a_s = 1 / (s + 1), for s = 0, 1, 2, 3.\n"
                "Coefficients: [1, 1/2, 1/3, 1/4].\n"
                "Execution: TimeIntegral::RungeKutta assigns each coefficient to ctrl.lhscoef "
                "before the corresponding residual/update stage.\n"
                "Limit: this is the implemented coefficient schedule, not a complete Butcher "
                "tableau and not evidence of formal order or stability.",
            ),
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
