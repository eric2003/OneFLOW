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

#include "TurbulenceBackend.h"

#include <memory>

BeginNameSpace( ONEFLOW )

// HIP execution space for the Spalart-Allmaras seam.
//
// A separate class from HipFluxBackend on purpose: the Euler/NS backend owns
// five-equation state and its kernels index fixed five-element arrays by
// equation, while SA is a one-equation scalar transport with a source.  Fusing
// them would either overflow those arrays on a sixth equation or route SA
// through UpdatePrimitiveStateKernel, which derives pressure from total energy
// and is meaningless for the working variable.
class HipTurbulenceBackend final : public TurbulenceBackend
{
public:
    static HipTurbulenceBackend & Shared();
    static void ReleaseShared();

    HipTurbulenceBackend();
    // Defined out of line: the private members below are incomplete types here,
    // and unique_ptr's deleter is only instantiated where they are complete.
    ~HipTurbulenceBackend() override;

    bool SupportsEddyViscosity() const override;
    void CalcEddyViscosity( const TurbSaEddyViscosityView & view ) override;

    bool SupportsSource() const override;
    void CalcSource( const TurbSaSourceView & view ) override;

    bool SupportsGradient() const override;
    void CalcGradient( const TurbSaGradientView & view ) override;

    bool SupportsFaceFlux() const override;
    void CalcFaceFlux( const TurbSaFaceView & view ) override;

    bool SupportsUpdate() const override;
    void Update( const TurbSaUpdateView & view ) override;

    // Implicit spectrum (UTurbSpectrum::CalcSpectrum1Equ).
    bool SupportsSpectrum() const;
    void CalcSpectrum( const TurbSaSpectrumView & view );

    // Boundary ghosts (UTurbBcSolver::CalcBc).
    bool SupportsBoundary() const;
    void CalcBoundary( const TurbSaBoundaryView & view );

    // LU-SGS scalar sweeps (UTurbLusgs::LowerSweep / UpperSweep).
    bool SupportsSweep() const;
    void CalcSweep( const TurbSaSweepView & view );

    // --- Single-sample entries for the numeric contract test ---------------
    //
    // Each evaluates the shipped device function on one cell or one face and
    // returns the scalar, so the test can compare against the host formula on
    // identical inputs without a mesh.  The mean-flow backend does the same
    // with EvaluateDeviceRoe / EvaluateDeviceBarth.
    struct SourceSample
    {
        Real density = 0.0;
        Real workingVariable = 0.0;
        Real laminarViscosity = 0.0;
        Real velocityGradient[ 9 ] = {};
        Real workingVariableGradient[ 3 ] = {};
        Real lengthScale = 0.0;
        Real cellVolume = 0.0;
        SaModelConstants constants;
    };
    void EvaluateSourceSample(
        const SourceSample & sample,
        Real & residual, Real & implicitSpectrum ) const;

    struct EddyViscositySample
    {
        Real density = 0.0;
        Real workingVariable = 0.0;
        Real laminarViscosity = 0.0;
        SaModelConstants constants;
    };
    void EvaluateEddyViscositySample(
        const EddyViscositySample & sample, Real & eddyViscosity ) const;

    struct InviscidFaceSample
    {
        Real workingVariableLeft = 0.0;
        Real workingVariableRight = 0.0;
        Real densityLeft = 0.0;
        Real densityRight = 0.0;
        Real velocityLeft[ 3 ] = {};
        Real velocityRight[ 3 ] = {};
        Real normal[ 3 ] = {};
        Real meshVelocityNormal = 0.0;
        Real faceArea = 0.0;
        Real densityWeight = 1.0;
    };
    void EvaluateInviscidFaceSample(
        const InviscidFaceSample & sample, Real & flux ) const;

    struct NormalGradientSample
    {
        Real workingVariableLeft = 0.0;
        Real workingVariableRight = 0.0;
        Real gradientLeft[ 3 ] = {};
        Real gradientRight[ 3 ] = {};
        Real faceCenter[ 3 ] = {};
        Real cellLeft[ 3 ] = {};
        Real cellRight[ 3 ] = {};
        Real normal[ 3 ] = {};
        bool boundary = false;
        bool isInterface = false;
    };
    void EvaluateNormalGradientSample(
        const NormalGradientSample & sample, Real & normalGradient ) const;

private:
    // Cell-to-face CSR over owned cells, rebuilt only when the partition
    // changes.  The SA gradient and sweep kernels both walk a cell's incident
    // faces; without it they would scan the whole face array per cell, which is
    // O(nCells * nFaces) and unusable at 3D production size.
    //
    // Held behind a pointer because the concrete type owns DeviceBuffers and
    // roster of HIP types this header deliberately does not include: the seam
    // header stays reachable from every execution space.
    struct CellFaceAdjacency;
    std::unique_ptr< CellFaceAdjacency > cellFaceAdjacency_;

    // Resident inputs for the LU-SGS sweep.  The sweep runs once per solver per
    // step and the geometry, adjacency, level schedule, residual and spectrum
    // are all unchanged between the lower and upper halves, so re-uploading
    // ~35 MB per call would cost more than the sweep itself.
    struct SweepStorage;
    std::unique_ptr< SweepStorage > sweepStorage_;

    // Returns the device CSR, (re)building it when the partition changed.
    CellFaceAdjacency & EnsureCellFaceAdjacency(
        const int * leftCell, const int * rightCell,
        int nCells, int nFaces, int nBoundaryFaces );
};

EndNameSpace
