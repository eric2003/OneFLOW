#include "NsInvFluxPolicy.h"

BeginNameSpace( ONEFLOW )

Real RelativeFaceNormalVelocity(
    Real xfn, Real yfn, Real zfn, Real vfn,
    Real u, Real v, Real w )
{
    return xfn * u + yfn * v + zfn * w - vfn;
}

void EulerConvectiveFlux(
    Real rho, Real u, Real v, Real w, Real p, Real h,
    Real vnRel,
    Real xfn, Real yfn, Real zfn, Real vfn,
    Real & fluxRho,
    Real & fluxMx,
    Real & fluxMy,
    Real & fluxMz,
    Real & fluxE )
{
    const Real massFlux = rho * vnRel;
    fluxRho = massFlux;
    fluxMx  = massFlux * u + xfn * p;
    fluxMy  = massFlux * v + yfn * p;
    fluxMz  = massFlux * w + zfn * p;
    fluxE   = massFlux * h + vfn * p;
}

EndNameSpace
