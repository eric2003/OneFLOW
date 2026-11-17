#include "TurbInvFluxPolicy.h"
#include "HXMath.h"

BeginNameSpace( ONEFLOW )

void SplitRelativeNormalVelocity(
    Real vnRel, Real & vnPlus, Real & vnMinus )
{
    vnPlus = half * ( vnRel + ABS( vnRel ) );
    vnMinus = half * ( vnRel - ABS( vnRel ) );
}

Real TurbConvectiveDensityWeight( Real rhoCoefFlag, Real density )
{
    return rhoCoefFlag * one + ( one - rhoCoefFlag ) * density;
}

Real UpwindScalarConvectiveFlux(
    Real vnPlus, Real weightLeft, Real phiLeft,
    Real vnMinus, Real weightRight, Real phiRight )
{
    return vnPlus * weightLeft * phiLeft
         + vnMinus * weightRight * phiRight;
}

EndNameSpace
