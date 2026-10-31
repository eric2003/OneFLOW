#pragma once

#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

// --- k production / dissipation (existing local layer) ---
Real SstProdkBoussinesq(
    Real eddyViscosity, Real strainRateSquared, Real divergence,
    Real invReynolds, Real density, Real kineticEnergy );
Real SstProdkStrainOnly(
    Real eddyViscosity, Real strainRateSquared, Real invReynolds );
Real SstProdkVorticity(
    Real eddyViscosity, Real vorticitySquared,
    Real invReynolds, Real density, Real kineticEnergy, Real divergence );
Real SstDisskStandard(
    Real fBetas, Real betas, Real density, Real kineticEnergy, Real omega );
Real SstDisskDes(
    Real fBetas, Real density, Real kineticEnergy, Real lengthScale );
Real SstLimitProdk( Real production, Real pklim, Real dissipation );

// --- blending ---
// @doc-algorithm sst-omega
// @doc-stage blend
// @doc-pure
// @doc-formula phi = F1*phi1 + (1-F1)*phi2
Real SstBlend( Real f1, Real value1, Real value2 );

// --- omega production / dissipation / cross-diffusion ---
// @doc-algorithm sst-omega
// @doc-stage production-menter
// @doc-pure
// @doc-formula P_w = alpha * rho / mu_t * P_k * Re
Real SstProdwMenter(
    Real alphaw, Real density, Real eddyViscosity,
    Real prodk, Real reynolds );

// @doc-algorithm sst-omega
// @doc-stage production-wilcox
// @doc-pure
// @doc-formula P_w = alpha * omega / k * P_k
Real SstProdwWilcox(
    Real alphaw, Real omega, Real kineticEnergy, Real prodk );

// @doc-algorithm sst-omega
// @doc-stage dissipation
// @doc-pure
// @doc-formula D_w = f_beta * beta * rho * omega^2
Real SstDissw( Real fBeta, Real beta, Real density, Real omega );

// @doc-algorithm sst-omega
// @doc-stage cross-diff-menter
// @doc-pure
// @doc-formula CD = 2 (1-F1) rho sigma_w2 (grad k · grad omega) / omega
Real SstCrossDiffMenter(
    Real f1, Real density, Real sigw2, Real crossTerm, Real omega );

// @doc-algorithm sst-omega
// @doc-stage cross-diff-wilcox
// @doc-pure
// @doc-formula CD = sigd * rho * max(cross, 0) / omega
Real SstCrossDiffWilcox(
    Real sigd, Real density, Real crossTerm, Real omega );

EndNameSpace
