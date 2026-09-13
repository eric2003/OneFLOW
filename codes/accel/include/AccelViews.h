/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
\*---------------------------------------------------------------------------*/

#pragma once

#include "HXTypeBasic.h"

BeginNameSpace( ONEFLOW )

// These views are deliberately backend-neutral. The owning solver remains
// responsible for lifetime and layout; a future HIP/CUDA/Kokkos adapter only
// receives pointers and extents, not MRField internals.
//
// Data layout convention (multi-equation):
//   Equation-major: data[eq * nFaces + face]  (or data[eq * nCells + cell])
//   This matches the main solver's MRField and the port's CI(c,i,nx) convention.
//
// For nEquations == 1 (scalar convection), qLeft[face] is the scalar at that face.
// For nEquations >= 3 (Euler/NS), qLeft stores conserved variables
//   [rho, rho*u, rho*v, rho*w, rho*E] for each face, equation-major.
struct FaceStateView
{
    int nFaces = 0;
    int nEquations = 0;
    const Real * qLeft = nullptr;
    const Real * qRight = nullptr;
    const Real * xNormal = nullptr;
    const Real * yNormal = nullptr;
    const Real * zNormal = nullptr;
    const Real * meshVelocityNormal = nullptr;
    const Real * faceArea = nullptr;
    Real gamma = 1.4;  // ratio of specific heats (used when nEquations >= 3)
};

struct FaceFluxView
{
    int nFaces = 0;
    int nEquations = 0;
    Real * values = nullptr;
};

struct FaceConnectivityView
{
    int nFaces = 0;
    int nBoundaryFaces = 0;
    const int * leftCell = nullptr;
    const int * rightCell = nullptr;
};

struct ResidualView
{
    int nCells = 0;
    int nEquations = 0;
    Real * values = nullptr;
};

EndNameSpace
