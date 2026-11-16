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

// Domain-level seam for the Spalart-Allmaras one-equation transport solve.
//
// This is a parallel seam of the same shape as FluxBackend, not a widening of
// it.  The reasons are the ones that already keep viscous flux off the domain
// backend: the SA solve is a face operator over the turbulence solver's own
// fields, its state is a single transported scalar with a source, and the
// production HIP path reaches the concrete HIP class directly.  Adding SA
// entry points to EulerDomainBackend or to FluxBackend would either widen an
// interface that has one CPU implementation, or force every execution space
// to implement turbulence it does not have.
//
// Every entry is non-pure and fails closed.  A caller that asks an execution
// space without SA support gets a clear error rather than a residual that
// silently omits the term, which is the failure mode a partial device port
// produces: the residual would still look plausible and would not converge to
// the host answer.
class TurbulenceBackend
{
public:
    virtual ~TurbulenceBackend() = default;

    // --- Per-cell pieces -------------------------------------------------
    // Turbulent (eddy) viscosity from the working variable:
    // TurbSrcFlux::CalcCellVist1Equ (codes/turb/src/TurbSrcFlux.cpp:282)
    // through SaComputeEddyViscosity (TurbSaPolicy.cpp:113).
    virtual bool SupportsEddyViscosity() const { return false; }

    virtual void CalcEddyViscosity( const TurbSaEddyViscosityView & view )
    {
        ( void ) view;
        throw std::logic_error(
            "Spalart-Allmaras eddy viscosity is not implemented for this "
            "execution space" );
    }

    // Source term including the implicit spectrum the host folds into
    // spec_sa: TurbCom::CalcSrcSa (codes/turb/src/TurbCom.cpp:747) and
    // UTurbSrcFlux::Update1Equ (:242).
    virtual bool SupportsSource() const { return false; }

    virtual void CalcSource( const TurbSaSourceView & view )
    {
        ( void ) view;
        throw std::logic_error(
            "Spalart-Allmaras source term is not implemented for this "
            "execution space" );
    }

    // --- Per-face pieces -------------------------------------------------
    // Scalar upwind convection and scalar diffusion:
    // UTurbInvFlux::CalcInvFlux (codes/uturb/src/UTurbInvFlux.cpp:111) and
    // UTurbVisFlux::CalcVisFlux (codes/uturb/src/UTurbVisFlux.cpp:106).
    virtual bool SupportsFaceFlux() const { return false; }

    virtual void CalcFaceFlux( const TurbSaFaceView & view )
    {
        ( void ) view;
        throw std::logic_error(
            "Spalart-Allmaras face flux is not implemented for this execution "
            "space" );
    }

    // --- Cell gradient ---------------------------------------------------
    // Green-Gauss gradient of the working variable, as
    // UTurbGrad::CalcGrad reaches it through Grad::CalcGrad.
    virtual bool SupportsGradient() const { return false; }

    virtual void CalcGradient( const TurbSaGradientView & view )
    {
        ( void ) view;
        throw std::logic_error(
            "Spalart-Allmaras gradient is not implemented for this execution "
            "space" );
    }

    // --- Explicit update -------------------------------------------------
    // TurbUpdate::CalcFlowField1Equ (codes/turb/src/TurbUpdate.cpp:64).
    virtual bool SupportsUpdate() const { return false; }

    virtual void Update( const TurbSaUpdateView & view )
    {
        ( void ) view;
        throw std::logic_error(
            "Spalart-Allmaras update is not implemented for this execution "
            "space" );
    }
};

class UnimplementedTurbulenceBackend final : public TurbulenceBackend
{
};

EndNameSpace
