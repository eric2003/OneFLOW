/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\*---------------------------------------------------------------------------*/

#pragma once

#include "AccelViews.h"

BeginNameSpace( ONEFLOW )

// Execution-space independent primitive-face contract.
//
// These three symbols used to live in EulerCpuAdapter.h, which forced the HIP
// backend to include a header named after the CPU adapter purely to reach a
// struct it uses at one call site (HipFluxBackend.hip includes EulerCpuAdapter.h
// and never calls an adapter member). The adapter's own comment already says
// the conversion and the equation-major pack are execution-space independent,
// so the name was the only thing tying them to the host.
//
// They live here now: no CPU or accelerator type appears, and either side can
// include this without dragging in the other's implementation.

// For one face, primitive values are [rho, u, p] for 3 equations and
// [rho, u, v, w, p] for 5 equations. Face batches use equation-major input
// and output: [rho, rho*u, rho*E] or [rho, rho*u, rho*v, rho*w, rho*E].
void PrimitiveToConserved(
    const Real * primitive, int nEquations, Real gamma, Real * conserved );

struct PrimitiveFaceStateView
{
    int nFaces = 0;
    int nEquations = 0;
    const Real * primitiveLeft = nullptr;
    const Real * primitiveRight = nullptr;
    const Real * primitiveLeftComponents[ 5 ] = {};
    const Real * primitiveRightComponents[ 5 ] = {};
    const Real * xNormal = nullptr;
    const Real * yNormal = nullptr;
    const Real * zNormal = nullptr;
    const Real * meshVelocityNormal = nullptr;
    const Real * faceArea = nullptr;
    Real gamma = 1.4;
};

void PackPrimitiveFaceStates(
    const Real * primitiveLeft,
    const Real * primitiveRight,
    int nFaces,
    int nEquations,
    Real gamma,
    Real * conservedLeft,
    Real * conservedRight );

EndNameSpace
