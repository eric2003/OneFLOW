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
void SubsonicAcousticRiemannCombine(
    Real vnInner, Real soundInner,
    Real vnRef, Real soundRef,
    Real gamaMinusOne,
    Real & vnBoundary, Real & soundBoundary );

// @doc-algorithm farfield-riemann
// @doc-stage tangential
// @doc-pure
void TangentialVelocityFromNormalProjection(
    Real u, Real v, Real w,
    Real xfn, Real yfn, Real zfn, Real vn,
    Real & uTangential, Real & vTangential, Real & wTangential );

// @doc-algorithm farfield-riemann
// @doc-stage reconstruct-velocity
// @doc-pure
void ReconstructVelocityFromTangentialAndNormal(
    Real uTangential, Real vTangential, Real wTangential,
    Real xfn, Real yfn, Real zfn, Real vnBoundary,
    Real & u, Real & v, Real & w );

// @doc-algorithm farfield-riemann
// @doc-stage isentropic-density
// @doc-pure
Real DensityFromIsentropicSoundSpeed(
    Real soundBoundary, Real entropy, Real gama, Real gamaMinusOne );

// @doc-algorithm farfield-riemann
// @doc-stage isentropic-pressure
// @doc-pure
Real PressureFromSoundSpeedDensity( Real soundBoundary, Real density, Real gama );

// @doc-algorithm isothermal-wall
// @doc-stage ideal-gas-density
// @doc-pure
// @doc-formula rho = p / (R * T * rmw)
Real IdealGasDensityFromPressureTemperature(
    Real pressure, Real temperature, Real gasCoef, Real reciMolecularWeight );

// @doc-algorithm isothermal-wall
// @doc-stage temperature-image
// @doc-pure
// @doc-formula T' = 2 T_wall - T_inner, clamped to T_min
Real MirrorScalarWithFloor(
    Real wallValue, Real interiorValue, Real floorValue );

// @doc-algorithm isothermal-wall
// @doc-stage ghost-density
// @doc-pure
// @doc-formula rho' = max(0.5 rho_w, 2 rho_w - rho_inner)
Real MirrorDensityAboutWall(
    Real densityWall, Real densityInner );

EndNameSpace
