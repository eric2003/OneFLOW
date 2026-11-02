#include "NsBcPolicy.h"
#include "HXMath.h"
#include "NsInvFluxPolicy.h"
#include <cmath>

BeginNameSpace( ONEFLOW )

void ReflectVelocityRelativeToMovingFace(
    Real u, Real v, Real w,
    Real xfn, Real yfn, Real zfn, Real vfn,
    Real & uReflected, Real & vReflected, Real & wReflected )
{
    const Real vnRel = RelativeFaceNormalVelocity(
        xfn, yfn, zfn, vfn, u, v, w );
    uReflected = u - two * xfn * vnRel;
    vReflected = v - two * yfn * vnRel;
    wReflected = w - two * zfn * vnRel;
}

void ImageVelocityAboutWall(
    Real u, Real v, Real w,
    Real uWall, Real vWall, Real wWall,
    Real & uImage, Real & vImage, Real & wImage )
{
    uImage = - u + two * uWall;
    vImage = - v + two * vWall;
    wImage = - w + two * wWall;
}

Real ArithmeticAverage( Real a, Real b )
{
    return half * ( a + b );
}

void SubsonicAcousticRiemannCombine(
    Real vnInner, Real soundInner,
    Real vnRef, Real soundRef,
    Real gamaMinusOne,
    Real & vnBoundary, Real & soundBoundary )
{
    const Real rPlus = vnInner + two * soundInner / gamaMinusOne;
    const Real rMinus = vnRef - two * soundRef / gamaMinusOne;
    vnBoundary = half * ( rPlus + rMinus );
    soundBoundary = fourth * ( rPlus - rMinus ) * gamaMinusOne;
}

void TangentialVelocityFromNormalProjection(
    Real u, Real v, Real w,
    Real xfn, Real yfn, Real zfn, Real vn,
    Real & uTangential, Real & vTangential, Real & wTangential )
{
    uTangential = u - xfn * vn;
    vTangential = v - yfn * vn;
    wTangential = w - zfn * vn;
}

void ReconstructVelocityFromTangentialAndNormal(
    Real uTangential, Real vTangential, Real wTangential,
    Real xfn, Real yfn, Real zfn, Real vnBoundary,
    Real & u, Real & v, Real & w )
{
    u = uTangential + xfn * vnBoundary;
    v = vTangential + yfn * vnBoundary;
    w = wTangential + zfn * vnBoundary;
}

Real DensityFromIsentropicSoundSpeed(
    Real soundBoundary, Real entropy, Real gama, Real gamaMinusOne )
{
    return std::pow(
        soundBoundary * soundBoundary / ( entropy * gama ),
        one / gamaMinusOne );
}

Real PressureFromSoundSpeedDensity( Real soundBoundary, Real density, Real gama )
{
    return soundBoundary * soundBoundary * density / gama;
}

Real IdealGasDensityFromPressureTemperature(
    Real pressure, Real temperature, Real gasCoef, Real reciMolecularWeight )
{
    return pressure / ( gasCoef * temperature * reciMolecularWeight );
}

Real MirrorScalarWithFloor(
    Real wallValue, Real interiorValue, Real floorValue )
{
    const Real mirrored = two * wallValue - interiorValue;
    return ( mirrored < floorValue ) ? floorValue : mirrored;
}

Real MirrorDensityAboutWall(
    Real densityWall, Real densityInner )
{
    const Real mirrored = two * densityWall - densityInner;
    return MAX( half * densityWall, mirrored );
}

EndNameSpace
