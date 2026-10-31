#pragma once

#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

// --- k ---
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

// --- omega ---
Real SstBlend( Real f1, Real value1, Real value2 );
Real SstProdwMenter(
    Real alphaw, Real density, Real eddyViscosity,
    Real prodk, Real reynolds );
Real SstProdwWilcox(
    Real alphaw, Real omega, Real kineticEnergy, Real prodk );
Real SstDissw( Real fBeta, Real beta, Real density, Real omega );
Real SstCrossDiffMenter(
    Real f1, Real density, Real sigw2, Real crossTerm, Real omega );
Real SstCrossDiffWilcox(
    Real sigd, Real density, Real crossTerm, Real omega );

// --- F1 blending ---
// @doc-algorithm sst-f1
// @doc-stage cross-diff-scalar
// @doc-pure
// @doc-formula CD_raw = 2 rho sigw2 (grad k · grad omega) / omega
Real SstCrossDiffScalar(
    Real density, Real sigw2, Real crossTerm, Real omega );

// @doc-algorithm sst-f1
// @doc-stage cdkw-floor
// @doc-pure
Real SstCdkwWithFloor( Real crossDiff, Real floorValue );

// @doc-algorithm sst-f1
// @doc-stage arg1-terms
// @doc-pure
// @doc-formula arg1 = min( max(sqrt(k)/(beta* d omega), 500 nu/(omega d^2 Re)), 4 rho sigw2 k /(CD d^2) )
Real SstBlendingArg1(
    Real termSqrtKe, Real termViscous, Real termCrossDiff );

// @doc-algorithm sst-f1
// @doc-stage f1
// @doc-pure
// @doc-formula F1 = tanh(arg1^4)
Real SstF1FromArg1( Real arg1 );

// @doc-algorithm sst-f1
// @doc-stage cross-term
// @doc-pure
// @doc-formula cross = grad k · grad omega
Real SstGradientDotProduct(
    Real dkdx, Real dkdy, Real dkdz,
    Real dwdx, Real dwdy, Real dwdz );

EndNameSpace
