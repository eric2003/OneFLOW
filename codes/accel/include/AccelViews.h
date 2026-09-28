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

#include <stdexcept>

BeginNameSpace( ONEFLOW )

enum class FieldLayout
{
    EquationMajor,
    EntityMajor
};

enum class FieldRepresentation
{
    Primitive,
    Conserved,
    Residual
};

enum class FaceAreaPolicy
{
    BackendMultiplies,
    CallerMultiplies
};

// Solver-side capability metadata. Ownership remains with the solver; halo
// exchange remains outside the numerical kernel.
struct SolverDomainCapabilities
{
    int nEquations = 0;
    int nGhostCells = 0;
    int nHaloLayers = 0;
    bool hasFaceGeometry = false;
    bool hasFaceConnectivity = false;
    bool supportsStateUpload = false;
    bool supportsResidualAdd = false;

    bool SupportsEulerState() const
    {
        return nEquations == 3 || nEquations == 5;
    }
};

struct SolverFieldView
{
    int nEntities = 0;
    int nComponents = 0;
    Real * values = nullptr;
    FieldLayout layout = FieldLayout::EquationMajor;
    FieldRepresentation representation = FieldRepresentation::Conserved;
};

struct SolverConstFieldView
{
    int nEntities = 0;
    int nComponents = 0;
    const Real * values = nullptr;
    FieldLayout layout = FieldLayout::EquationMajor;
    FieldRepresentation representation = FieldRepresentation::Conserved;
};

struct FaceGeometryView
{
    int nFaces = 0;
    const Real * xNormal = nullptr;
    const Real * yNormal = nullptr;
    const Real * zNormal = nullptr;
    const Real * meshVelocityNormal = nullptr;
    const Real * faceArea = nullptr;
    FaceAreaPolicy areaPolicy = FaceAreaPolicy::BackendMultiplies;
};

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
    FieldLayout layout = FieldLayout::EquationMajor;
    FieldRepresentation representation = FieldRepresentation::Conserved;
    FaceAreaPolicy areaPolicy = FaceAreaPolicy::BackendMultiplies;
};

inline void ValidateFaceStateView( const FaceStateView & state )
{
    if ( state.nFaces <= 0 || state.nEquations <= 0
         || state.qLeft == nullptr || state.qRight == nullptr
         || state.gamma <= 1.0 )
    {
        throw std::invalid_argument( "invalid solver face state view" );
    }
    if ( state.layout != FieldLayout::EquationMajor
         || state.representation != FieldRepresentation::Conserved
         || state.areaPolicy != FaceAreaPolicy::BackendMultiplies )
    {
        throw std::invalid_argument(
            "unsupported solver face state contract" );
    }
}

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
    // Optional explicit boundary mask. Without it, boundary faces occupy
    // [0, nBoundaryFaces) for backward compatibility.
    const unsigned char * boundaryMask = nullptr;
};

inline void ValidateFaceConnectivityView(
    const FaceConnectivityView & connectivity )
{
    if ( connectivity.nFaces <= 0
         || connectivity.nBoundaryFaces < 0
         || connectivity.nBoundaryFaces > connectivity.nFaces
         || connectivity.leftCell == nullptr
         || connectivity.rightCell == nullptr )
    {
        throw std::invalid_argument(
            "invalid solver face connectivity view" );
    }
}

struct ResidualView
{
    int nCells = 0;
    int nEquations = 0;
    Real * values = nullptr;
};

inline void ValidateSolverFieldView( const SolverConstFieldView & field )
{
    if ( field.nEntities <= 0 || field.nComponents <= 0
         || field.values == nullptr )
    {
        throw std::invalid_argument( "invalid solver field view" );
    }
}

inline void ValidateFaceGeometryView( const FaceGeometryView & geometry )
{
    if ( geometry.nFaces <= 0 || geometry.faceArea == nullptr )
    {
        throw std::invalid_argument( "invalid solver face geometry view" );
    }
}

EndNameSpace
