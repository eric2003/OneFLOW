#include "TurbSpectrumPolicy.h"
#include "HXMath.h"

BeginNameSpace( ONEFLOW )

Real TurbInviscidFaceSpectralRadius( Real faceArea, Real absVnMid )
{
    return faceArea * half * absVnMid;
}

void TurbInviscidFaceMatrixContrib(
    Real faceArea, Real vnLeft, Real vnRight, Real absVnMid,
    Real & matrixLeftAdd, Real & matrixRightAdd )
{
    matrixLeftAdd = faceArea * half * ( - vnLeft - absVnMid );
    matrixRightAdd = faceArea * half * ( vnRight - absVnMid );
}

Real SaEffectiveKinematicViscosity(
    Real laminarViscosity, Real density, Real nuTilde )
{
    return laminarViscosity / ( density + SMALL ) + nuTilde;
}

Real TurbDiffusiveFaceSpectralRadius(
    Real invReynolds, Real invSigma, Real nuEff,
    Real faceAreaSquared, Real invCellVolume )
{
    return invReynolds * invCellVolume * ( invSigma * nuEff * faceAreaSquared );
}

EndNameSpace
