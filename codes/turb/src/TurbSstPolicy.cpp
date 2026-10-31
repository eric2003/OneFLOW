#include "TurbSstPolicy.h"
#include "HXMath.h"
#include <cmath>

BeginNameSpace( ONEFLOW )

Real SstProdkBoussinesq(
    Real eddyViscosity, Real strainRateSquared, Real divergence,
    Real invReynolds, Real density, Real kineticEnergy )
{
    return eddyViscosity * ( strainRateSquared - two3rd * SQR( divergence ) ) * invReynolds
         - two3rd * density * kineticEnergy * divergence;
}

Real SstProdkStrainOnly(
    Real eddyViscosity, Real strainRateSquared, Real invReynolds )
{
    return eddyViscosity * strainRateSquared * invReynolds;
}

Real SstProdkVorticity(
    Real eddyViscosity, Real vorticitySquared,
    Real invReynolds, Real density, Real kineticEnergy, Real divergence )
{
    return eddyViscosity * vorticitySquared * invReynolds
         - two3rd * density * kineticEnergy * divergence;
}

Real SstDisskStandard(
    Real fBetas, Real betas, Real density, Real kineticEnergy, Real omega )
{
    return fBetas * betas * density * kineticEnergy * omega;
}

Real SstDisskDes(
    Real fBetas, Real density, Real kineticEnergy, Real lengthScale )
{
    return fBetas * density * kineticEnergy * std::sqrt( kineticEnergy )
         / lengthScale;
}

Real SstLimitProdk( Real production, Real pklim, Real dissipation )
{
    return MIN( production, pklim * dissipation );
}

EndNameSpace
