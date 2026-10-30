#pragma once

#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

// @doc-algorithm spectral-radius
// @doc-stage inviscid-face
// @doc-pure
// @doc-formula lambda_inv = (1/2) * A * (|vn| + c)
Real InviscidFaceSpectralRadius(
    Real vnRelMid, Real soundSpeedMid, Real faceArea );

// @doc-algorithm spectral-radius
// @doc-stage viscous-face-model1
// @doc-pure
// Blazek-style viscous spectral estimate with max of momentum and energy diffusivities.
Real ViscousFaceSpectralRadiusModel1(
    Real densityMid,
    Real laminarViscosity, Real turbulentViscosity,
    Real invPrandtlLaminar, Real invPrandtlTurbulent,
    Real gama, Real reynolds, Real faceArea );

// @doc-algorithm spectral-radius
// @doc-stage viscous-face-model2
// @doc-pure
Real ViscousFaceSpectralRadiusModel2(
    Real densityMid, Real totalViscosity,
    Real cellDistance, Real reynolds, Real faceArea );

// @doc-algorithm spectral-radius
// @doc-stage viscous-face-model3
// @doc-pure
Real ViscousFaceSpectralRadiusModel3(
    Real densityMid, Real reynolds, Real faceArea );

// @doc-algorithm spectral-radius
// @doc-stage cell-inv-dt
// @doc-pure
// @doc-formula dt = CFL * Vol / lambda
Real CellTimeStepFromSpectrum(
    Real cfl, Real cellVolume, Real spectralRadius );

// @doc-algorithm spectral-radius
// @doc-stage harmonic-combine
// @doc-pure
// @doc-formula dt <- dt * dt_vis / (dt + dt_vis)
Real CombineInvAndVisTimeStep(
    Real invTimeStep, Real visTimeStep );

EndNameSpace
