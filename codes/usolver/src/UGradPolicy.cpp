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


Real FaceSkewAngleDegrees(
    Real dx, Real dy, Real dz, Real projectedLength )
{
    Real result = projectedLength / ( DIST( dx, dy, dz ) + SMALL );
    if ( result > 1.0 ) result = 1.0;
    if ( result < -1.0 ) result = -1.0;
    return asin( result ) * 180.0 / PI;
}

void UnitizeDisplacement(
    Real & dx, Real & dy, Real & dz, Real & invLength )
{
    invLength = 1.0 / DIST( dx, dy, dz );
    dx *= invLength;
    dy *= invLength;
    dz *= invLength;
}


void CorrectFaceGradient(
    Real qLeft, Real qRight,
    Real & dqdx, Real & dqdy, Real & dqdz,
    Real unitDx, Real unitDy, Real unitDz, Real invLength )
{
    const Real projected = dqdx * unitDx + dqdy * unitDy + dqdz * unitDz;
    const Real target = ( qRight - qLeft ) * invLength;
    const Real df = target - projected;
    dqdx += df * unitDx;
    dqdy += df * unitDy;
    dqdz += df * unitDz;
}

EndNameSpace
