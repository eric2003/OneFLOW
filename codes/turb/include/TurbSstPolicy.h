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


// --- Wilcox / EASM beta* and beta modifiers ---
// @doc-algorithm wilcox-fbeta
// @doc-stage xk
// @doc-pure
// @doc-formula x_k = (grad k · grad omega) / omega^3
Real WilcoxXkFromCrossTerm( Real crossTerm, Real omega );

// @doc-algorithm wilcox-fbeta
// @doc-stage fk
// @doc-pure
// @doc-formula f_k = 1 (xk<=0); (1+680 xk^2)/(1+400 xk^2) (xk>0)
Real WilcoxFkFromXk( Real xk );

// @doc-algorithm wilcox-fbeta
// @doc-stage fw-1998
// @doc-pure
// @doc-formula f_w = (1+70 x_w)/(1+80 x_w)
Real WilcoxFw1998( Real xw );

// @doc-algorithm wilcox-fbeta
// @doc-stage fw-2006
// @doc-pure
// @doc-formula f_w = (1+85 x_w)/(1+100 x_w)
Real WilcoxFw2006( Real xw );

// @doc-algorithm wilcox-fbeta
// @doc-stage vortex-stretch
// @doc-pure
// Magnitude used in x_w before dividing by (beta* omega)^3
Real WilcoxVortexStretchingInvariant(
    Real w12, Real w13, Real w23,
    Real s11, Real s22, Real s33, Real s12, Real s13, Real s23 );

EndNameSpace
