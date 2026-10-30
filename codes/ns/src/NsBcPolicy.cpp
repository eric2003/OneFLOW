#include "NsBcPolicy.h"
#include "HXMath.h"
#include "NsInvFluxPolicy.h"

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

EndNameSpace
