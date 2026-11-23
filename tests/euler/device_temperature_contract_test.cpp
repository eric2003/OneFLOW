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

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

// Device temperature field contract.
//
// The viscous port needs a temperature field the device does not have today:
// UpdatePrimitiveStateKernel derives pressure from total energy and never
// builds tempr.  This test pins the new kernel to the host relation rather than
// merely checking that a buffer became non-null.
//
// Host reference: NsCalcGamaT (codes/ns/src/NsRhs.cpp) writes, for every cell
// it covers,
//
//     T[cell] = p[cell] / ( nscom.statecoef * rho[cell] * oamw )
//
// with oamw the molecular-weight ratio (the local `one` outside the chemistry
// path) and nscom.statecoef = 1 / ( gama_ref * mach_ref^2 ).  Gamma therefore
// reaches the relation through the coefficient, not as a literal, and the test
// passes the coefficient instead of 1.4.
//
// Scope: a small deterministic tetrahedron-like cell block, one shared
// geometry and connectivity set, inner and ghost cells compared separately.
// This is a component contract, not a substitute for a viscous residual gate.

#include "AccelRuntime.h"
#include "HipFluxBackend.h"
// Ns3DDeviceState and Ns3DBackendState were split out of EulerDomain.h into
// their own header; the test constructs the state and attaches its backend.
#include "Ns3DDomain.h"

#include <gtest/gtest.h>

#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <vector>

using namespace ONEFLOW;

namespace
{

// The host routine's relation, transcribed once here so the expectation is
// independent of the kernel under test.
double HostTemperature(
    double density, double pressure,
    double stateCoefficient, double molecularWeightRatio )
{
    return pressure
        / ( stateCoefficient * density * molecularWeightRatio );
}

struct Mesh
{
    int nCells = 4;
    int nBoundaryFaces = 4;
    int nFaces = 8;
    int nTotalCells = 8;

    std::vector< Real > xCell;
    std::vector< Real > yCell;
    std::vector< Real > zCell;
    std::vector< Real > xFace;
    std::vector< Real > yFace;
    std::vector< Real > zFace;
    std::vector< Real > xNormal;
    std::vector< Real > yNormal;
    std::vector< Real > zNormal;
    std::vector< Real > faceArea;
    std::vector< Real > cellVolume;
    std::vector< int > leftCell;
    std::vector< int > rightCell;

    // Four unit cubes laid along x.  Boundary faces point at their ghost on the
    // right; internal faces point at the neighbour.  The ghost cell centres sit
    // one cell width outside, which is what makes the gradient weights
    // non-degenerate.
    Mesh()
        : xCell{ 0.5, 1.5, 2.5, 3.5, -0.5, 4.5, 5.5, 6.5 }
        , yCell{ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 }
        , zCell{ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 }
        , xFace{ 0.0, 1.0, 2.0, 3.0, 4.0, 2.0, 1.0, 0.0 }
        , yFace{ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 }
        , zFace{ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 }
        , xNormal{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0 }
        , yNormal{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 }
        , zNormal{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
        , faceArea{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 }
        , cellVolume{ 1.0, 1.0, 1.0, 1.0 }
        , leftCell{ 0, 1, 2, 3, 0, 1, 2, 3 }
        , rightCell{ 4, 5, 6, 7, 1, 2, 3, 4 }
    {
    }
};

// Deterministic primitive state.  Density and pressure vary per cell so a
// component mix-up cannot cancel, and the ghost cells carry their own values
// so the temperature kernel is exercised over the full resident extent rather
// than only where the solver writes.
std::vector< Real > PrimitiveState( const Mesh & mesh )
{
    std::vector< Real > primitive(
        static_cast< std::size_t >( 5 ) * mesh.nTotalCells, 0.0 );
    for ( int cell = 0; cell < mesh.nTotalCells; ++ cell )
    {
        primitive[ 0 * mesh.nTotalCells + cell ] = 1.0 + 0.25 * cell;
        primitive[ 1 * mesh.nTotalCells + cell ] = 0.10 * cell;
        primitive[ 2 * mesh.nTotalCells + cell ] = -0.05 * cell;
        primitive[ 3 * mesh.nTotalCells + cell ] = 0.02 * cell;
        primitive[ 4 * mesh.nTotalCells + cell ] = 1.0 + 0.40 * cell;
    }
    return primitive;
}

// gama_ref and mach_ref from the flagship configuration, not a hard-coded 1.4
// substitute: the coefficient is what the solver initialises.
constexpr double kGamaRef = 1.4;
constexpr double kMachRef = 0.84;
constexpr double kStateCoefficient = 1.0 / ( kGamaRef * kMachRef * kMachRef );
constexpr double kMolecularWeightRatio = 1.0;
constexpr double kTolerance = 1.0e-14;

class DeviceTemperatureContract : public ::testing::Test
{
protected:
    static void SetUpTestSuite()
    {
        // The production default is the CPU backend, whose Allocate() is plain
        // host memory and would never exercise the device path.  The contract
        // is about the HIP kernel, so select HIP explicitly before the runtime
        // reads its configuration.
        setenv( "ONEFLOW_ACCEL_BACKEND", "HIP", 1 );
        InitializeAccelRuntime( 0, 1 );
        // A missing HIP registration would otherwise surface as every case
        // being "skipped" by gtest's suite-level exception handling, which reads
        // as success.  Throw instead: a silently skipped device contract is the
        // exact failure mode this test exists to prevent.  gtest turns this
        // into a failed suite, not a skip.
        if ( AccelRuntime::Instance().Backend().Kind() != AccelBackendKind::HIP )
        {
            throw std::runtime_error(
                "the temperature contract requires the HIP backend; build this "
                "target with ONEFLOW_ENABLE_HIP so HipBackend.cpp registers it" );
        }
    }

    static void TearDownTestSuite()
    {
        AccelRuntime::Instance().Finalize();
        unsetenv( "ONEFLOW_ACCEL_BACKEND" );
    }

    // BindState rejects a key whose device does not match the device the
    // runtime actually selected, and the selection depends on the visible DCU
    // count, so the state key has to be built from the live context rather than
    // a hard-coded 0.
    static int SelectedDevice()
    {
        return AccelRuntime::Instance().Context().selectedDevice;
    }
};

}

TEST_F( DeviceTemperatureContract, DeviceTemperatureMatchesHostRelation )
{
    const Mesh mesh;
    std::vector< Real > primitive = PrimitiveState( mesh );

    const EulerDomainProblem problem{
        mesh.nCells, mesh.nBoundaryFaces, 5, kGamaRef, 1.0e-3, 1.0,
        EulerDomainBoundary::Periodic };
    const EulerDomainStateKey key{
        0, 0, 0, AccelBackendKind::HIP, SelectedDevice() };
    Ns3DDeviceState owner( problem, key );
    owner.ReserveFaces( mesh.nFaces );
    owner.AttachBackendState( CreateHipNs3DBackendState( SelectedDevice() ) );

    HipFluxBackend backend;
    // The fallback storage keeps `backend` usable for this single-grid test;
    // binding additionally proves the derived fields land in the state-owned
    // host mirror, which is what the viscous reader will consume.
    const void * cacheKey = &owner;
    backend.BindState( cacheKey, owner );

    std::vector< Real > temperature( mesh.nTotalCells, 0.0 );
    std::vector< Real > dTdx( mesh.nTotalCells, 0.0 );
    std::vector< Real > dTdy( mesh.nTotalCells, 0.0 );
    std::vector< Real > dTdz( mesh.nTotalCells, 0.0 );

    CellTemperatureView view;
    view.nCells = mesh.nCells;
    view.nGhostCells = mesh.nBoundaryFaces;
    view.nFaces = mesh.nFaces;
    view.nBoundaryFaces = mesh.nBoundaryFaces;
    view.nEquations = 5;
    view.xFace = mesh.xFace.data();
    view.yFace = mesh.yFace.data();
    view.zFace = mesh.zFace.data();
    view.xNormal = mesh.xNormal.data();
    view.yNormal = mesh.yNormal.data();
    view.zNormal = mesh.zNormal.data();
    view.faceArea = mesh.faceArea.data();
    view.xCell = mesh.xCell.data();
    view.yCell = mesh.yCell.data();
    view.zCell = mesh.zCell.data();
    view.cellVolume = mesh.cellVolume.data();
    view.leftCell = mesh.leftCell.data();
    view.rightCell = mesh.rightCell.data();
    view.stateCoefficient = kStateCoefficient;
    view.molecularWeightRatio = kMolecularWeightRatio;
    view.temperature = temperature.data();
    view.dTdx = dTdx.data();
    view.dTdy = dTdy.data();
    view.dTdz = dTdz.data();
    view.cacheKey = cacheKey;
    view.downloadToHost = true;
    for ( int equation = 0; equation < view.nEquations; ++ equation )
    {
        view.q[ equation ] =
            primitive.data() + equation * mesh.nTotalCells;
    }

    EXPECT_TRUE( backend.SupportsCellTemperature() );
    backend.CalcCellTemperature( view );

    // Inner cells and ghost cells both matter: the viscous flux reads the
    // ghost side and the host routine populates both ranges.
    for ( int cell = 0; cell < mesh.nTotalCells; ++ cell )
    {
        const double density = primitive[ 0 * mesh.nTotalCells + cell ];
        const double pressure = primitive[ 4 * mesh.nTotalCells + cell ];
        const double expected = HostTemperature(
            density, pressure, kStateCoefficient, kMolecularWeightRatio );
        EXPECT_NEAR( temperature[ cell ], expected, kTolerance )
            << "temperature at cell " << cell;
        EXPECT_TRUE( std::isfinite( temperature[ cell ] ) );
    }

    // The derived field must also be in the state-owned mirror, because that is
    // the copy a host reader and the viscous host fallback consume.
    for ( int cell = 0; cell < mesh.nTotalCells; ++ cell )
    {
        EXPECT_NEAR(
            owner.temperature[ cell ], temperature[ cell ], 0.0 )
            << "state mirror at cell " << cell;
    }
}

// Host transcription of the discrete temperature gradient.  This is the same
// operator GradientCellKernel applies to the mean flow, written out here so the
// expectation does not come from the kernel under test.  Only the operator
// matters, not smoothness: the device is required to reproduce this exact
// discrete sum, component by component, for the supplied mesh.
void HostTemperatureGradient(
    const Mesh & mesh,
    const std::vector< Real > & temperature,
    std::vector< Real > & dTdx,
    std::vector< Real > & dTdy,
    std::vector< Real > & dTdz )
{
    dTdx.assign( mesh.nTotalCells, 0.0 );
    dTdy.assign( mesh.nTotalCells, 0.0 );
    dTdz.assign( mesh.nTotalCells, 0.0 );

    // Per-cell incident faces, built the way the backend builds its CSR: a
    // boundary face touches only its owner, an internal face touches both.
    std::vector< std::vector< int > > incident( mesh.nCells );
    for ( int face = 0; face < mesh.nFaces; ++ face )
    {
        const int left = mesh.leftCell[ face ];
        const int right = mesh.rightCell[ face ];
        if ( left >= 0 && left < mesh.nCells ) incident[ left ].push_back( face );
        if ( right >= 0 && right < mesh.nCells ) incident[ right ].push_back( face );
    }

    for ( int cell = 0; cell < mesh.nCells; ++ cell )
    {
        Real ax = 0.0;
        Real ay = 0.0;
        Real az = 0.0;
        for ( int face : incident[ cell ] )
        {
            const int left = mesh.leftCell[ face ];
            const int right = mesh.rightCell[ face ];
            const Real dxl = mesh.xFace[ face ] - mesh.xCell[ left ];
            const Real dyl = mesh.yFace[ face ] - mesh.yCell[ left ];
            const Real dzl = mesh.zFace[ face ] - mesh.zCell[ left ];
            const Real dxr = mesh.xFace[ face ] - mesh.xCell[ right ];
            const Real dyr = mesh.yFace[ face ] - mesh.yCell[ right ];
            const Real dzr = mesh.zFace[ face ] - mesh.zCell[ right ];
            const Real delt1 = std::sqrt( dxl * dxl + dyl * dyl + dzl * dzl );
            const Real delt2 = std::sqrt( dxr * dxr + dyr * dyr + dzr * dzr );
            const Real delta = 1.0 / ( delt1 + delt2 + 1.0e-40 );
            const Real cl = delt2 * delta;
            const Real cr = delt1 * delta;
            const Real nx = mesh.xNormal[ face ] * mesh.faceArea[ face ];
            const Real ny = mesh.yNormal[ face ] * mesh.faceArea[ face ];
            const Real nz = mesh.zNormal[ face ] * mesh.faceArea[ face ];
            const Real sign = left == cell ? 1.0 : -1.0;
            const Real value = cl * temperature[ left ] + cr * temperature[ right ];
            ax += sign * nx * value;
            ay += sign * ny * value;
            az += sign * nz * value;
        }
        const Real inverseVolume = 1.0 / mesh.cellVolume[ cell ];
        dTdx[ cell ] = ax * inverseVolume;
        dTdy[ cell ] = ay * inverseVolume;
        dTdz[ cell ] = az * inverseVolume;
    }

    // Ghosts copy their owner's gradient, matching the boundary-copy kernel.
    for ( int face = 0; face < mesh.nBoundaryFaces; ++ face )
    {
        const int owner = mesh.leftCell[ face ];
        const int ghost = mesh.rightCell[ face ];
        dTdx[ ghost ] = dTdx[ owner ];
        dTdy[ ghost ] = dTdy[ owner ];
        dTdz[ ghost ] = dTdz[ owner ];
    }
}

TEST_F( DeviceTemperatureContract, DeviceGradientMatchesHostOperator )
{
    const Mesh mesh;
    std::vector< Real > primitive = PrimitiveState( mesh );
    const double stateCoefficient = kStateCoefficient;
    const double molecularWeightRatio = kMolecularWeightRatio;

    const EulerDomainProblem problem{
        mesh.nCells, mesh.nBoundaryFaces, 5, kGamaRef, 1.0e-3, 1.0,
        EulerDomainBoundary::Periodic };
    const EulerDomainStateKey key{
        0, 1, 0, AccelBackendKind::HIP, SelectedDevice() };
    Ns3DDeviceState owner( problem, key );
    owner.ReserveFaces( mesh.nFaces );
    owner.AttachBackendState( CreateHipNs3DBackendState( SelectedDevice() ) );

    HipFluxBackend backend;
    const void * cacheKey = &owner;
    backend.BindState( cacheKey, owner );

    std::vector< Real > temperature( mesh.nTotalCells, 0.0 );
    std::vector< Real > dTdx( mesh.nTotalCells, 0.0 );
    std::vector< Real > dTdy( mesh.nTotalCells, 0.0 );
    std::vector< Real > dTdz( mesh.nTotalCells, 0.0 );

    CellTemperatureView view;
    view.nCells = mesh.nCells;
    view.nGhostCells = mesh.nBoundaryFaces;
    view.nFaces = mesh.nFaces;
    view.nBoundaryFaces = mesh.nBoundaryFaces;
    view.nEquations = 5;
    view.xFace = mesh.xFace.data();
    view.yFace = mesh.yFace.data();
    view.zFace = mesh.zFace.data();
    view.xNormal = mesh.xNormal.data();
    view.yNormal = mesh.yNormal.data();
    view.zNormal = mesh.zNormal.data();
    view.faceArea = mesh.faceArea.data();
    view.xCell = mesh.xCell.data();
    view.yCell = mesh.yCell.data();
    view.zCell = mesh.zCell.data();
    view.cellVolume = mesh.cellVolume.data();
    view.leftCell = mesh.leftCell.data();
    view.rightCell = mesh.rightCell.data();
    view.stateCoefficient = stateCoefficient;
    view.molecularWeightRatio = molecularWeightRatio;
    view.temperature = temperature.data();
    view.dTdx = dTdx.data();
    view.dTdy = dTdy.data();
    view.dTdz = dTdz.data();
    view.cacheKey = cacheKey;
    view.downloadToHost = true;
    for ( int equation = 0; equation < view.nEquations; ++ equation )
    {
        view.q[ equation ] =
            primitive.data() + equation * mesh.nTotalCells;
    }

    backend.CalcCellTemperature( view );

    // The host reference uses the temperature the kernel must have produced,
    // so the two halves of the check are chained: first T == host relation,
    // then dT == host operator applied to that T.
    std::vector< Real > hostTemperature( mesh.nTotalCells, 0.0 );
    for ( int cell = 0; cell < mesh.nTotalCells; ++ cell )
    {
        hostTemperature[ cell ] = HostTemperature(
            primitive[ 0 * mesh.nTotalCells + cell ],
            primitive[ 4 * mesh.nTotalCells + cell ],
            stateCoefficient, molecularWeightRatio );
    }

    std::vector< Real > hostDx;
    std::vector< Real > hostDy;
    std::vector< Real > hostDz;
    HostTemperatureGradient( mesh, hostTemperature, hostDx, hostDy, hostDz );

    // The nonconstant p/rho ratio gives a genuinely varying temperature, so a
    // gradient implementation that always returns zero cannot pass.
    EXPECT_GT( std::abs( hostDx[ 0 ] ) + std::abs( hostDy[ 0 ] )
                   + std::abs( hostDz[ 0 ] ),
               0.0 );

    for ( int cell = 0; cell < mesh.nTotalCells; ++ cell )
    {
        EXPECT_NEAR( dTdx[ cell ], hostDx[ cell ], kTolerance )
            << "dTdx at cell " << cell;
        EXPECT_NEAR( dTdy[ cell ], hostDy[ cell ], kTolerance )
            << "dTdy at cell " << cell;
        EXPECT_NEAR( dTdz[ cell ], hostDz[ cell ], kTolerance )
            << "dTdz at cell " << cell;
    }
}

TEST_F( DeviceTemperatureContract, UnboundBackendReportsCapability )
{
    // The capability query must not throw even when the execution space cannot
    // build temperature; only the entry point fails closed.
    HipFluxBackend backend;
    EXPECT_TRUE( backend.SupportsCellTemperature() );
}

TEST( CellTemperatureViewContract, RejectsIncompleteWiring )
{
    Real q[ 8 ] = {};
    Real geometry[ 8 ] = {};
    Real temperature[ 8 ] = {};
    int cells[ 4 ] = { 0, 1, 2, 3 };

    CellTemperatureView view;
    view.nCells = 4;
    view.nGhostCells = 4;
    view.nFaces = 4;
    view.nBoundaryFaces = 4;
    view.nEquations = 5;
    for ( int equation = 0; equation < view.nEquations; ++ equation )
    {
        view.q[ equation ] = q;
    }
    view.xFace = geometry;
    view.yFace = geometry;
    view.zFace = geometry;
    view.xNormal = geometry;
    view.yNormal = geometry;
    view.zNormal = geometry;
    view.faceArea = geometry;
    view.xCell = geometry;
    view.yCell = geometry;
    view.zCell = geometry;
    view.cellVolume = geometry;
    view.leftCell = cells;
    view.rightCell = cells;
    view.stateCoefficient = 1.0;
    view.molecularWeightRatio = 1.0;
    view.temperature = temperature;
    view.dTdx = temperature;
    view.dTdy = temperature;
    view.dTdz = temperature;
    view.cacheKey = q;

    EXPECT_NO_THROW( ValidateCellTemperatureView( view ) );

    // A missing thermodynamic coefficient would divide by a compile-time
    // assumption instead of the solved value, so it is a hard rejection.
    view.stateCoefficient = 0.0;
    EXPECT_THROW( ValidateCellTemperatureView( view ), std::invalid_argument );
    view.stateCoefficient = 1.0;

    view.temperature = nullptr;
    EXPECT_THROW( ValidateCellTemperatureView( view ), std::invalid_argument );
    view.temperature = temperature;

    // A missing primitive plane is rejected here rather than surfacing later
    // from the resident-state binder.
    view.q[ 2 ] = nullptr;
    EXPECT_THROW( ValidateCellTemperatureView( view ), std::invalid_argument );
}
