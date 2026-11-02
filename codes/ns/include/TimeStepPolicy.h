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
// @doc-formula c1=(4/3)*(mu+mut); c2=gamma*(mu/Pr_l+mut/Pr_t);
// @doc-formula lambda_vis = A^2 * 2*max(c1,c2)/(Re*rho)
// Blazek-style viscous spectral estimate: max of momentum and energy
// diffusivities, scaled by face area squared (time-step bound, not exact
// Jacobian eigenvalue). See doc/source/numerical/viscous_spectral_radius_model1.rst.
Real ViscousFaceSpectralRadiusModel1(
    Real densityMid,
    Real laminarViscosity, Real turbulentViscosity,
    Real invPrandtlLaminar, Real invPrandtlTurbulent,
    Real gama, Real reynolds, Real faceArea );

// @doc-algorithm spectral-radius
// @doc-stage viscous-face-model2
// @doc-pure
// @doc-formula c1=2*mu_tot/(rho*d*Re); lambda_vis=0.5*c1*A
// Length-scale form: d is cell/wall distance. See
// doc/source/numerical/viscous_spectral_radius_model2.rst.
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
