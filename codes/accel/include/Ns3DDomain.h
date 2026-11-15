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

// Navier--Stokes 3D state types. They were carved out of EulerDomain.h so that
// the Euler/1D contract does not carry the 3D mean-flow state shape. The
// dependency stays one-directional: this header includes EulerDomain.h for the
// base class and the shared problem/key/view types, never the other way round.

#include "AccelBackend.h"
#include "EulerDomain.h"
#include "HXTypeBasic.h"

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <vector>

BeginNameSpace( ONEFLOW )

// Backend-specific allocations are hidden behind this ownership seam. A HIP
// implementation may derive from it and own DeviceBuffer instances without
// exposing HIP types in the solver-facing API.
class Ns3DBackendState
{
public:
    virtual ~Ns3DBackendState() = default;
    virtual AccelBackendKind Backend() const noexcept = 0;
    virtual int DeviceId() const noexcept = 0;
};

// Long-lived 3D main-solver state metadata. This is a separate specialization
// rather than forcing the 1D port state shape onto Navier--Stokes. The HIP
// backend-specific allocations are owned through an opaque state object here,
// aligned with solver/grid lifecycle boundaries.
class Ns3DDeviceState final : public EulerDomainState
{
public:
    Ns3DDeviceState( const EulerDomainProblem & problemValue,
                     const EulerDomainStateKey & keyValue )
        : problem( problemValue ), key( keyValue )
    {
    }

    void ReserveFaces( int faceCount )
    {
        if ( faceCount < 0 )
            throw std::invalid_argument( "negative Ns3D face count" );
        if ( nFaces != faceCount ) MarkTopologyChanged();
        nFaces = faceCount;
        faceState.resize( problem.nEquations * faceCount );
        faceFlux.resize( problem.nEquations * faceCount );
        xNormal.resize( faceCount );
        yNormal.resize( faceCount );
        zNormal.resize( faceCount );
        meshVelocityNormal.resize( faceCount );
        faceArea.resize( faceCount );
        leftCell.resize( faceCount );
        rightCell.resize( faceCount );
        boundaryMetadata.resize( faceCount );
        haloMetadata.resize( faceCount );
    }

    void Upload( const EulerDomainConstFieldView & field )
    {
        ValidateEulerDomainField( problem, field );
        const std::size_t count =
            static_cast< std::size_t >( problem.nCells ) * problem.nEquations;
        conservedState.assign( field.values, field.values + count );
        oldState = conservedState;
        residual.assign( count, 0.0 );
        rkScratch.assign( count, 0.0 );
        gradient.assign( count * 3, 0.0 );
        limiter.assign( count, 0.0 );
        reconstruction.assign( count, 0.0 );
        viscousTurbulence.assign( count, 0.0 );
        uploaded = true;
        MarkFieldChanged();
    }

    void AttachBackendState( std::unique_ptr< Ns3DBackendState > state )
    {
        if ( state == nullptr )
        {
            throw std::invalid_argument(
                "cannot attach a null Ns3D backend state" );
        }
        if ( state->Backend() != key.backend )
        {
            throw std::invalid_argument(
                "Ns3D backend state kind does not match its registry key" );
        }
        if ( state->DeviceId() != key.deviceId )
        {
            throw std::invalid_argument(
                "Ns3D backend state device does not match its registry key" );
        }
        backendState = std::move( state );
    }

    Ns3DBackendState * BackendState() noexcept
    {
        return backendState.get();
    }

    const Ns3DBackendState * BackendState() const noexcept
    {
        return backendState.get();
    }

    bool HasBackendState() const noexcept
    {
        return backendState != nullptr;
    }

    void ReleaseBackendState() noexcept
    {
        backendState.reset();
    }

    const void * OwnerToken() const noexcept
    {
        return this;
    }

    std::uint64_t TopologyGeneration() const noexcept
    {
        return topologyGeneration;
    }

    std::uint64_t FieldGeneration() const noexcept
    {
        return fieldGeneration;
    }

    void MarkTopologyChanged() noexcept
    {
        ++topologyGeneration;
    }

    void MarkFieldChanged() noexcept
    {
        ++fieldGeneration;
    }

    EulerDomainProblem problem;
    EulerDomainStateKey key;
    int nFaces = 0;
    bool uploaded = false;

    std::vector< Real > conservedState;
    std::vector< Real > oldState;
    std::vector< Real > residual;
    std::vector< Real > rkScratch;
    std::vector< Real > gradient;
    std::vector< Real > limiter;
    std::vector< Real > reconstruction;
    std::vector< Real > faceState;
    std::vector< Real > faceFlux;
    std::vector< Real > xNormal;
    std::vector< Real > yNormal;
    std::vector< Real > zNormal;
    std::vector< Real > meshVelocityNormal;
    std::vector< Real > faceArea;
    std::vector< int > leftCell;
    std::vector< int > rightCell;
    std::vector< unsigned char > boundaryMetadata;
    std::vector< unsigned char > haloMetadata;
    std::vector< Real > viscousTurbulence;

private:
    std::unique_ptr< Ns3DBackendState > backendState;
    std::uint64_t topologyGeneration = 1;
    std::uint64_t fieldGeneration = 1;
};

EndNameSpace
