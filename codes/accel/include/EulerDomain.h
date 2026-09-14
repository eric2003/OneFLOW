/*---------------------------------------------------------------------------*\\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
\\*---------------------------------------------------------------------------*/

#pragma once

#include "AccelBackend.h"
#include "HXTypeBasic.h"

#include <memory>
#include <stdexcept>

BeginNameSpace( ONEFLOW )

// This is the solver-facing contract extracted from the standalone 1D Euler
// port. It deliberately carries views and extents, not MRField or zone
// ownership. The solver remains responsible for adapting its fields to these
// views and for deciding when a state is created or destroyed.
enum class EulerDomainBoundary
{
    Periodic,
    Wall,
    Inflow,
    Outflow
};

enum class EulerDomainRunMode
{
    NoTrace,
    FullTrace
};

struct EulerDomainProblem
{
    int nCells = 0;
    int nGhostCells = 0;
    int nEquations = 3;
    Real gamma = 1.4;
    Real dt = 0.0;
    Real dx = 0.0;
    EulerDomainBoundary boundary = EulerDomainBoundary::Periodic;
};

struct EulerDomainFieldView
{
    int nCells = 0;
    int nEquations = 0;
    Real * values = nullptr;
};

struct EulerDomainConstFieldView
{
    int nCells = 0;
    int nEquations = 0;
    const Real * values = nullptr;
};

struct EulerDomainRunOptions
{
    EulerDomainRunMode mode = EulerDomainRunMode::NoTrace;
    void * trace = nullptr;
    void * stats = nullptr;
};

// One state per solver/zone/grid/backend. A solver index alone is not enough
// because a multiblock solve can revisit the same solver on different local
// zones or grid levels.
struct EulerDomainStateKey
{
    int solverIndex = -1;
    int localZoneId = -1;
    int gridLevel = -1;
    AccelBackendKind backend = AccelBackendKind::CPU;

    bool operator==( const EulerDomainStateKey & other ) const
    {
        return solverIndex == other.solverIndex
            && localZoneId == other.localZoneId
            && gridLevel == other.gridLevel
            && backend == other.backend;
    }
};

inline void ValidateEulerDomainProblem( const EulerDomainProblem & problem )
{
    if ( problem.nCells <= 0 || problem.nGhostCells < 0
         || problem.nEquations != 3 || problem.gamma <= 1.0
         || problem.dt <= 0.0 || problem.dx <= 0.0 )
    {
        throw std::invalid_argument( "invalid OneFLOW Euler domain problem" );
    }
}

inline void ValidateEulerDomainField(
    const EulerDomainProblem & problem,
    const EulerDomainConstFieldView & field )
{
    ValidateEulerDomainProblem( problem );
    if ( field.nCells != problem.nCells
         || field.nEquations != problem.nEquations
         || field.values == nullptr )
    {
        throw std::invalid_argument( "inconsistent OneFLOW Euler field view" );
    }
}

inline void ValidateEulerDomainField(
    const EulerDomainProblem & problem,
    const EulerDomainFieldView & field )
{
    ValidateEulerDomainProblem( problem );
    if ( field.nCells != problem.nCells
         || field.nEquations != problem.nEquations
         || field.values == nullptr )
    {
        throw std::invalid_argument( "inconsistent OneFLOW Euler field view" );
    }
}

class EulerDomainState
{
public:
    virtual ~EulerDomainState() = default;
};

class EulerDomainBackend
{
public:
    virtual ~EulerDomainBackend() = default;

    virtual const char * Name() const = 0;
    virtual bool IsAccelerator() const = 0;
    virtual std::unique_ptr< EulerDomainState > CreateState(
        const EulerDomainProblem & problem,
        const EulerDomainStateKey & key ) const = 0;
    virtual void Upload(
        EulerDomainState & state,
        const EulerDomainConstFieldView & field ) const = 0;
    virtual void Advance(
        EulerDomainState & state,
        int steps,
        const EulerDomainRunOptions & options ) const = 0;
    virtual void Download(
        const EulerDomainState & state,
        EulerDomainFieldView & field ) const = 0;
};

EndNameSpace
