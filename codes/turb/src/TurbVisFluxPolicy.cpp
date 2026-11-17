#include "TurbVisFluxPolicy.h"
#include "HXMath.h"

BeginNameSpace( ONEFLOW )

Real ArithmeticMean( Real a, Real b )
{
    return half * ( a + b );
}

Real SaFaceEffectiveDiffusivity(
    Real laminarKinematicViscosity,
    Real workingVariableAverage,
    Real invSigma )
{
    return ( laminarKinematicViscosity + workingVariableAverage ) * invSigma;
}

Real SstFaceEffectiveDiffusivity(
    Real laminarViscosity, Real turbulentViscosity, Real sigma )
{
    return laminarViscosity + turbulentViscosity * sigma;
}

Real TurbScalarDiffusiveFlux(
    Real invReynolds, Real effectiveDiffusivity,
    Real normalGradient, Real faceArea )
{
    return - invReynolds * effectiveDiffusivity * normalGradient * faceArea;
}

EndNameSpace
