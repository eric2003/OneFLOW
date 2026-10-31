#pragma once

#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

// @doc-algorithm sst-k
// @doc-stage production-boussinesq
// @doc-pure
// @doc-formula P_k = mu_t (Sij:Sij - 2/3 (div v)^2)/Re - 2/3 rho k div v
Real SstProdkBoussinesq(
    Real eddyViscosity, Real strainRateSquared, Real divergence,
    Real invReynolds, Real density, Real kineticEnergy );

// @doc-algorithm sst-k
// @doc-stage production-sij
// @doc-pure
// @doc-formula P_k = mu_t * Sij:Sij / Re  (SST-2003)
Real SstProdkStrainOnly(
    Real eddyViscosity, Real strainRateSquared, Real invReynolds );

// @doc-algorithm sst-k
// @doc-stage production-vorticity
// @doc-pure
// @doc-formula P_k = mu_t * Omega:Omega / Re - 2/3 rho k div v
Real SstProdkVorticity(
    Real eddyViscosity, Real vorticitySquared,
    Real invReynolds, Real density, Real kineticEnergy, Real divergence );

// @doc-algorithm sst-k
// @doc-stage dissipation
// @doc-pure
// @doc-formula eps = f_beta* * beta* * rho * k * omega
Real SstDisskStandard(
    Real fBetas, Real betas, Real density, Real kineticEnergy, Real omega );

// @doc-algorithm sst-k
// @doc-stage dissipation-des
// @doc-pure
// @doc-formula eps = f_beta* * rho * k^{3/2} / L
Real SstDisskDes(
    Real fBetas, Real density, Real kineticEnergy, Real lengthScale );

// @doc-algorithm sst-k
// @doc-stage limit-production
// @doc-pure
// @doc-formula P_k <- min(P_k, pklim * eps)
Real SstLimitProdk( Real production, Real pklim, Real dissipation );

EndNameSpace
