#pragma once

#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

// @doc-algorithm boundary-velocity
// @doc-stage symmetry-reflect
// @doc-pure
void ReflectVelocityRelativeToMovingFace(
    Real u, Real v, Real w,
    Real xfn, Real yfn, Real zfn, Real vfn,
    Real & uReflected, Real & vReflected, Real & wReflected );

// @doc-algorithm boundary-velocity
// @doc-stage wall-image
// @doc-pure
void ImageVelocityAboutWall(
    Real u, Real v, Real w,
    Real uWall, Real vWall, Real wWall,
    Real & uImage, Real & vImage, Real & wImage );

// @doc-algorithm boundary-velocity
// @doc-stage face-average
// @doc-pure
Real ArithmeticAverage( Real a, Real b );

// @doc-algorithm farfield-riemann
// @doc-stage acoustic-combine
// @doc-pure
// @doc-formula R+/- = vn +/- 2 c/(gamma-1); vn_b = (R++R-)/2; c_b = (R+-R-)*(gamma-1)/4
void SubsonicAcousticRiemannCombine(
    Real vnInner, Real soundInner,
    Real vnRef, Real soundRef,
    Real gamaMinusOne,
    Real & vnBoundary, Real & soundBoundary );

// @doc-algorithm farfield-riemann
// @doc-stage tangential
// @doc-pure
// @doc-formula v_t = v - n * vn
void TangentialVelocityFromNormalProjection(
    Real u, Real v, Real w,
    Real xfn, Real yfn, Real zfn, Real vn,
    Real & uTangential, Real & vTangential, Real & wTangential );

// @doc-algorithm farfield-riemann
// @doc-stage reconstruct-velocity
// @doc-pure
// @doc-formula v = v_t + n * vn_b
void ReconstructVelocityFromTangentialAndNormal(
    Real uTangential, Real vTangential, Real wTangential,
    Real xfn, Real yfn, Real zfn, Real vnBoundary,
    Real & u, Real & v, Real & w );

// @doc-algorithm farfield-riemann
// @doc-stage isentropic-density
// @doc-pure
// @doc-formula rho = (c^2 / (s * gamma))^{1/(gamma-1)}
Real DensityFromIsentropicSoundSpeed(
    Real soundBoundary, Real entropy, Real gama, Real gamaMinusOne );

// @doc-algorithm farfield-riemann
// @doc-stage isentropic-pressure
// @doc-pure
// @doc-formula p = c^2 * rho / gamma
Real PressureFromSoundSpeedDensity( Real soundBoundary, Real density, Real gama );

EndNameSpace
