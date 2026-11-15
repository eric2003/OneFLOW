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

#include "AccelViews.h"
#include <stdexcept>

BeginNameSpace( ONEFLOW )

// Domain-level seam for the first production accelerator operation. It is
// intentionally independent from HIP/CUDA/Kokkos so the main CFD code can
// select an implementation at phase boundaries rather than per face.
class FluxBackend
{
public:
    virtual ~FluxBackend() = default;

    virtual void CalcInvFlux(
        const FaceStateView & state,
        FaceFluxView & flux,
        int scheme ) = 0;

    virtual void AddFaceFlux(
        const FaceFluxView & flux,
        const FaceConnectivityView & connectivity,
        ResidualView & residual ) = 0;

    // Viscous and heat-conduction flux for one stage of the mean-flow residual.
    //
    // This belongs on the flux seam rather than on a concrete backend: the
    // production HIP path currently reaches HipFluxBackend::Shared() directly
    // for every stage kernel, and adding viscous flux in that style would leave
    // CUDA and Kokkos with nothing to implement against.  The payload is the
    // backend-neutral ViscousFluxView, so an execution space opts in by
    // overriding and needs no new plumbing.
    //
    // Non-pure on purpose. A viscous case is a capability, not a requirement:
    // making it pure would force throwaway stubs into every existing
    // implementation - including the CPU one, whose host chain computes viscous
    // flux outside this seam - and would break call sites for no gain. The
    // default fails closed, so a caller that asks anyway gets a clear error
    // instead of a residual that silently omits the viscous terms.
    virtual bool SupportsViscousFlux() const { return false; }

    virtual void CalcViscousFaceFlux( const ViscousFluxView & view )
    {
        ( void ) view;
        throw std::logic_error(
            "viscous flux is not implemented for this execution space" );
    }
};

class UnimplementedFluxBackend final : public FluxBackend
{
public:
    void CalcInvFlux(
        const FaceStateView &,
        FaceFluxView &,
        int ) override
    {
        throw std::logic_error(
            "OneFLOW flux backend is not implemented for this execution space." );
    }

    void AddFaceFlux(
        const FaceFluxView &,
        const FaceConnectivityView &,
        ResidualView & ) override
    {
        throw std::logic_error(
            "OneFLOW residual backend is not implemented for this execution space." );
    }
};

EndNameSpace
