#include "NsLusgsPolicy.h"
#include "HXMath.h"

BeginNameSpace( ONEFLOW )

void LusgsSplitFaceEigenvalues(
    Real vnRel, Real soundSpeed, Real signOfMatrix,
    Real faceAbs1, Real faceAbs2, Real faceAbs3,
    Real & lmd1, Real & lmd2, Real & lmd3 )
{
    lmd1 = half * ( vnRel + signOfMatrix * faceAbs1 );
    lmd2 = half * ( vnRel + soundSpeed + signOfMatrix * faceAbs2 );
    lmd3 = half * ( vnRel - soundSpeed + signOfMatrix * faceAbs3 );
}

void LusgsLinearizationWeights(
    Real lmd1, Real lmd2, Real lmd3,
    Real soundSpeedSquared, Real soundSpeed,
    Real & x1, Real & x2 )
{
    x1 = ( lmd1 + lmd1 - lmd2 - lmd3 ) / ( soundSpeedSquared + soundSpeedSquared );
    x2 = ( lmd2 - lmd3 ) / ( soundSpeed + soundSpeed );
}

void LusgsFluxIncrementFiveEq(
    Real lmd1, Real term1, Real term2,
    Real u, Real v, Real w, Real totalEnthalpy, Real vnFluid,
    Real xfn, Real yfn, Real zfn, Real faceArea,
    Real dqRho, Real dqRhoU, Real dqRhoV, Real dqRhoW, Real dqRhoE,
    Real & dfRho, Real & dfMx, Real & dfMy, Real & dfMz, Real & dfE )
{
    dfRho = ( lmd1 * dqRho - term1 ) * faceArea;
    dfMx  = ( lmd1 * dqRhoU - u * term1 + xfn * term2 ) * faceArea;
    dfMy  = ( lmd1 * dqRhoV - v * term1 + yfn * term2 ) * faceArea;
    dfMz  = ( lmd1 * dqRhoW - w * term1 + zfn * term2 ) * faceArea;
    dfE   = ( lmd1 * dqRhoE - totalEnthalpy * term1 + vnFluid * term2 ) * faceArea;
}

EndNameSpace
