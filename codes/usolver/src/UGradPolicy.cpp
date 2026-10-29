#include "UGradPolicy.h"

BeginNameSpace( ONEFLOW )

Real ArithmeticFaceAverage( Real qLeft, Real qRight )
{
    return half * ( qLeft + qRight );
}

void InverseDistanceFaceWeights(
    Real dxl, Real dyl, Real dzl,
    Real dxr, Real dyr, Real dzr,
    Real & weightLeft, Real & weightRight )
{
    const Real delt1 = DIST( dxl, dyl, dzl );
    const Real delt2 = DIST( dxr, dyr, dzr );
    const Real delta = 1.0 / ( delt1 + delt2 + SMALL );
    weightLeft = delt2 * delta;
    weightRight = delt1 * delta;
}

Real WeightedFaceAverage( Real qLeft, Real qRight, Real weightLeft, Real weightRight )
{
    return weightLeft * qLeft + weightRight * qRight;
}

void AreaWeightedFaceNormal(
    Real xfn, Real yfn, Real zfn, Real area,
    Real & snx, Real & sny, Real & snz )
{
    snx = xfn * area;
    sny = yfn * area;
    snz = zfn * area;
}


void ScaleGradientByVolume(
    Real volume,
    Real & dqdx, Real & dqdy, Real & dqdz )
{
    const Real ovol = one / volume;
    dqdx *= ovol;
    dqdy *= ovol;
    dqdz *= ovol;
}

void CopyGradient(
    Real srcDqdx, Real srcDqdy, Real srcDqdz,
    Real & dstDqdx, Real & dstDqdy, Real & dstDqdz )
{
    dstDqdx = srcDqdx;
    dstDqdy = srcDqdy;
    dstDqdz = srcDqdz;
}

void AccumulateFaceGradient(
    Real snx, Real sny, Real snz, Real faceValue, Real sign,
    Real & dqdx, Real & dqdy, Real & dqdz )
{
    const Real w = sign * faceValue;
    dqdx += snx * w;
    dqdy += sny * w;
    dqdz += snz * w;
}

EndNameSpace
