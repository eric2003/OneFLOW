/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\*---------------------------------------------------------------------------*/

#include "AccelRuntime.h"
#include "CpuFluxBackend.h"
#include "EulerCpuAdapter.h"
#include "HipFluxBackend.h"
#include "HipKernel.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <limits>
#include <vector>

namespace
{

double MaxDiff( const std::vector< ONEFLOW::Real > & a,
                const std::vector< ONEFLOW::Real > & b )
{
    double maxVal = 0.0;
    for ( std::size_t i = 0; i < a.size(); ++ i )
    {
        maxVal = std::max( maxVal,
            static_cast< double >( std::abs( a[ i ] - b[ i ] ) ) );
    }
    return maxVal;
}

// --- Scalar convection test (existing) ---

bool TestScalarConvection()
{
    constexpr int nFaces = 513;
    constexpr int nCells = 1026;
    constexpr int nBoundaryFaces = 17;
    std::vector< ONEFLOW::Real > qLeft( nFaces );
    std::vector< ONEFLOW::Real > qRight( nFaces );
    std::vector< ONEFLOW::Real > xNormal( nFaces );
    std::vector< ONEFLOW::Real > area( nFaces );
    std::vector< int > leftCell( nFaces );
    std::vector< int > rightCell( nFaces );
    for ( int face = 0; face < nFaces; ++ face )
    {
        qLeft[ face ] = 0.25 + 0.01 * face;
        qRight[ face ] = 0.75 - 0.003 * face;
        xNormal[ face ] = ( face % 3 == 0 ) ? -1.0 : 1.0;
        area[ face ] = 0.5 + 0.002 * ( face % 11 );
        leftCell[ face ] = face;
        rightCell[ face ] = nFaces + face;
    }

    ONEFLOW::FaceStateView state;
    state.nFaces = nFaces;
    state.nEquations = 1;
    state.qLeft = qLeft.data();
    state.qRight = qRight.data();
    state.xNormal = xNormal.data();
    state.faceArea = area.data();

    std::vector< ONEFLOW::Real > cpuFlux( nFaces );
    std::vector< ONEFLOW::Real > hipFlux( nFaces );
    ONEFLOW::FaceFluxView cpuFluxView{ nFaces, 1, cpuFlux.data() };
    ONEFLOW::FaceFluxView hipFluxView{ nFaces, 1, hipFlux.data() };

    ONEFLOW::CpuFluxBackend cpuBackend;
    ONEFLOW::HipFluxBackend hipBackend;
    cpuBackend.CalcInvFlux( state, cpuFluxView, 0 );
    hipBackend.CalcInvFlux( state, hipFluxView, 0 );

    double fluxError = MaxDiff( cpuFlux, hipFlux );
    if ( fluxError > 1.0e-15 )
    {
        std::fprintf( stderr, "Scalar flux FAIL: error %.3e\n", fluxError );
        return false;
    }

    std::vector< ONEFLOW::Real > cpuResidual( nCells, 0.0 );
    std::vector< ONEFLOW::Real > hipResidual( nCells, 0.0 );
    ONEFLOW::FaceConnectivityView conn;
    conn.nFaces = nFaces;
    conn.nBoundaryFaces = nBoundaryFaces;
    conn.leftCell = leftCell.data();
    conn.rightCell = rightCell.data();
    ONEFLOW::ResidualView cpuResView{ nCells, 1, cpuResidual.data() };
    ONEFLOW::ResidualView hipResView{ nCells, 1, hipResidual.data() };
    cpuBackend.AddFaceFlux( cpuFluxView, conn, cpuResView );
    hipBackend.AddFaceFlux( hipFluxView, conn, hipResView );

    double resError = MaxDiff( cpuResidual, hipResidual );
    if ( resError > 1.0e-15 )
    {
        std::fprintf( stderr, "Scalar residual FAIL: error %.3e\n", resError );
        return false;
    }

    std::printf( "OneFLOW HIP scalar flux: PASS (flux %.3e, residual %.3e)\n",
        fluxError, resError );
    return true;
}

// --- Euler Rusanov test (new) ---

bool TestEulerRusanov()
{
    constexpr int nFaces = 128;
    constexpr int nCells = 256;
    constexpr int nBoundaryFaces = 8;
    constexpr int nEq = 3;  // 1D Euler

    // Create a smooth 1D Euler initial condition (same as contract test)
    constexpr double kDx = 1.0 / nCells;
    constexpr double kPi = 3.14159265358979323846;
    constexpr double kGamma = 1.4;

    // Cell-centered conserved state
    std::vector< ONEFLOW::Real > cellState( nCells * nEq );
    for ( int cell = 0; cell < nCells; ++ cell )
    {
        const double x = ( cell + 0.5 ) * kDx;
        const double density = 1.0 + 0.08 * std::sin( 2.0 * kPi * x );
        const double velocity = 0.2 + 0.04 * std::cos( 2.0 * kPi * x );
        const double pressure = 1.0 + 0.05 * std::sin( 4.0 * kPi * x );
        cellState[ 0 * nCells + cell ] = density;
        cellState[ 1 * nCells + cell ] = density * velocity;
        cellState[ 2 * nCells + cell ] = pressure / ( kGamma - 1.0 )
            + 0.5 * density * velocity * velocity;
    }

    // Build face data from cell-centered (simple average for left/right)
    std::vector< ONEFLOW::Real > qLeft( nFaces * nEq );
    std::vector< ONEFLOW::Real > qRight( nFaces * nEq );
    std::vector< ONEFLOW::Real > xNormal( nFaces, 1.0 );
    std::vector< ONEFLOW::Real > area( nFaces, 1.0 );
    std::vector< int > leftCell( nFaces );
    std::vector< int > rightCell( nFaces );

    for ( int face = 0; face < nFaces; ++ face )
    {
        const int il = face;
        const int ir = ( face + 1 ) % nCells;
        leftCell[ face ] = il;
        rightCell[ face ] = ir;
        for ( int eq = 0; eq < nEq; ++ eq )
        {
            qLeft[ eq * nFaces + face ] = cellState[ eq * nCells + il ];
            qRight[ eq * nFaces + face ] = cellState[ eq * nCells + ir ];
        }
    }

    ONEFLOW::FaceStateView state;
    state.nFaces = nFaces;
    state.nEquations = nEq;
    state.qLeft = qLeft.data();
    state.qRight = qRight.data();
    state.xNormal = xNormal.data();
    state.faceArea = area.data();
    state.gamma = kGamma;

    std::vector< ONEFLOW::Real > cpuFlux( nFaces * nEq );
    std::vector< ONEFLOW::Real > hipFlux( nFaces * nEq );
    ONEFLOW::FaceFluxView cpuFluxView{ nFaces, nEq, cpuFlux.data() };
    ONEFLOW::FaceFluxView hipFluxView{ nFaces, nEq, hipFlux.data() };

    ONEFLOW::CpuFluxBackend cpuBackend;
    ONEFLOW::HipFluxBackend hipBackend;
    cpuBackend.CalcInvFlux( state, cpuFluxView, 0 );
    hipBackend.CalcInvFlux( state, hipFluxView, 0 );

    double fluxError = MaxDiff( cpuFlux, hipFlux );
    // Euler flux uses more fp ops, tolerance relaxed slightly
    if ( fluxError > 1.0e-14 )
    {
        std::fprintf( stderr, "Euler flux FAIL: error %.3e\n", fluxError );
        return false;
    }

    // Test residual accumulation
    std::vector< ONEFLOW::Real > cpuResidual( nCells * nEq, 0.0 );
    std::vector< ONEFLOW::Real > hipResidual( nCells * nEq, 0.0 );
    ONEFLOW::FaceConnectivityView conn;
    conn.nFaces = nFaces;
    conn.nBoundaryFaces = nBoundaryFaces;
    conn.leftCell = leftCell.data();
    conn.rightCell = rightCell.data();
    ONEFLOW::ResidualView cpuResView{ nCells, nEq, cpuResidual.data() };
    ONEFLOW::ResidualView hipResView{ nCells, nEq, hipResidual.data() };
    cpuBackend.AddFaceFlux( cpuFluxView, conn, cpuResView );
    hipBackend.AddFaceFlux( hipFluxView, conn, hipResView );

    double resError = MaxDiff( cpuResidual, hipResidual );
    if ( resError > 1.0e-14 )
    {
        std::fprintf( stderr, "Euler residual FAIL: error %.3e\n", resError );
        return false;
    }

    // Physicality check: flux should be finite
    for ( int i = 0; i < nFaces * nEq; ++ i )
    {
        if ( !std::isfinite( cpuFlux[ i ] ) )
        {
            std::fprintf( stderr, "Euler flux non-finite at %d\n", i );
            return false;
        }
    }

    std::printf( "OneFLOW HIP Euler flux: PASS (flux %.3e, residual %.3e, %d faces, %d eq)\n",
        fluxError, resError, nFaces, nEq );
    return true;
}

bool TestMainSolverFiveEquationLaxFriedrichs()
{
    constexpr int nFaces = 257;
    constexpr int nCells = 300;
    constexpr int nEq = 5;
    constexpr double gamma = 1.4;
    std::vector< ONEFLOW::Real > primitiveLeft( nFaces * nEq );
    std::vector< ONEFLOW::Real > primitiveRight( nFaces * nEq );
    std::vector< ONEFLOW::Real > xNormal( nFaces );
    std::vector< ONEFLOW::Real > yNormal( nFaces );
    std::vector< ONEFLOW::Real > zNormal( nFaces );
    std::vector< ONEFLOW::Real > meshVelocityNormal( nFaces );
    std::vector< ONEFLOW::Real > area( nFaces );

    for ( int face = 0; face < nFaces; ++ face )
    {
        const double phase = 0.03125 * face;
        double nx = std::cos( phase );
        double ny = std::sin( phase );
        double nz = 0.25 * std::sin( 0.7 * phase );
        const double normalMagnitude = std::sqrt( nx * nx + ny * ny + nz * nz );
        xNormal[ face ] = nx / normalMagnitude;
        yNormal[ face ] = ny / normalMagnitude;
        zNormal[ face ] = nz / normalMagnitude;
        meshVelocityNormal[ face ] = 0.02 * std::cos( 0.3 * phase );
        area[ face ] = 0.5 + 0.01 * ( face % 19 );

        primitiveLeft[ 0 * nFaces + face ] = 1.0 + 0.05 * std::sin( phase );
        primitiveLeft[ 1 * nFaces + face ] = 0.7 + 0.08 * std::cos( phase );
        primitiveLeft[ 2 * nFaces + face ] = -0.2 + 0.04 * std::sin( 0.5 * phase );
        primitiveLeft[ 3 * nFaces + face ] = 0.1 + 0.03 * std::cos( 0.8 * phase );
        primitiveLeft[ 4 * nFaces + face ] = 1.0 + 0.06 * std::cos( 1.1 * phase );
        primitiveRight[ 0 * nFaces + face ] = 0.95 + 0.04 * std::cos( 0.9 * phase );
        primitiveRight[ 1 * nFaces + face ] = 0.6 + 0.07 * std::sin( 1.2 * phase );
        primitiveRight[ 2 * nFaces + face ] = -0.1 + 0.05 * std::cos( 0.6 * phase );
        primitiveRight[ 3 * nFaces + face ] = 0.15 + 0.02 * std::sin( 0.4 * phase );
        primitiveRight[ 4 * nFaces + face ] = 0.9 + 0.05 * std::sin( phase );
    }

    ONEFLOW::PrimitiveFaceStateView state;
    state.nFaces = nFaces;
    state.nEquations = nEq;
    state.primitiveLeft = primitiveLeft.data();
    state.primitiveRight = primitiveRight.data();
    state.xNormal = xNormal.data();
    state.yNormal = yNormal.data();
    state.zNormal = zNormal.data();
    state.meshVelocityNormal = meshVelocityNormal.data();
    state.faceArea = area.data();
    state.gamma = gamma;

    std::vector< ONEFLOW::Real > cpuFlux( nFaces * nEq );
    std::vector< ONEFLOW::Real > hipFlux( nFaces * nEq );
    ONEFLOW::FaceFluxView cpuFluxView{ nFaces, nEq, cpuFlux.data() };
    ONEFLOW::FaceFluxView hipFluxView{ nFaces, nEq, hipFlux.data() };
    ONEFLOW::CpuFluxBackend cpuBackend;
    ONEFLOW::HipFluxBackend hipBackend;
    ONEFLOW::EulerCpuAdapter adapter;
    adapter.CalcInvFlux( state, cpuFluxView, cpuBackend, 1 );
    adapter.CalcInvFlux( state, hipFluxView, hipBackend, 1 );

    const double fluxError = MaxDiff( cpuFlux, hipFlux );
    if ( fluxError > 1.0e-12 )
    {
        std::fprintf( stderr,
            "Main solver 5-equation Lax-Friedrichs flux FAIL: %.3e\n",
            fluxError );
        return false;
    }
    for ( ONEFLOW::Real value : hipFlux )
    {
        if ( ! std::isfinite( value ) )
        {
            std::fprintf( stderr, "Main solver HIP flux is non-finite\n" );
            return false;
        }
    }

    std::vector< int > leftCell( nFaces );
    std::vector< int > rightCell( nFaces );
    std::vector< unsigned char > boundaryMask( nFaces );
    int nBoundaryFaces = 0;
    for ( int face = 0; face < nFaces; ++ face )
    {
        leftCell[ face ] = face % nCells;
        rightCell[ face ] = ( face * 17 + 3 ) % nCells;
        if ( rightCell[ face ] == leftCell[ face ] )
            rightCell[ face ] = ( rightCell[ face ] + 1 ) % nCells;
        boundaryMask[ face ] = face % 13 == 7 ? 1 : 0;
        nBoundaryFaces += boundaryMask[ face ] != 0;
    }

    std::vector< ONEFLOW::Real > cpuResidual( nCells * nEq );
    std::vector< ONEFLOW::Real > hipResidual( nCells * nEq );
    for ( int i = 0; i < nCells * nEq; ++ i )
    {
        cpuResidual[ i ] = 1.0e-6 * ( i % 7 );
        hipResidual[ i ] = cpuResidual[ i ];
    }
    ONEFLOW::FaceConnectivityView connectivity{
        nFaces, nBoundaryFaces, leftCell.data(), rightCell.data(),
        boundaryMask.data() };
    ONEFLOW::ResidualView cpuResidualView{
        nCells, nEq, cpuResidual.data() };
    ONEFLOW::ResidualView hipResidualView{
        nCells, nEq, hipResidual.data() };
    cpuBackend.AddFaceFlux( cpuFluxView, connectivity, cpuResidualView );
    hipBackend.AddCurrentFaceFlux( connectivity, hipResidualView );

    std::vector< ONEFLOW::Real > fusedFlux( nFaces * nEq );
    std::vector< ONEFLOW::Real > fusedResidual =
        std::vector< ONEFLOW::Real >( nCells * nEq );
    for ( int i = 0; i < nCells * nEq; ++ i )
    {
        fusedResidual[ i ] = 1.0e-6 * ( i % 7 );
    }
    ONEFLOW::FaceFluxView fusedFluxView{
        nFaces, nEq, fusedFlux.data() };
    ONEFLOW::ResidualView fusedResidualView{
        nCells, nEq, fusedResidual.data() };
    ONEFLOW::HipFluxBackend fusedBackend;
    fusedBackend.CalcAndAddPrimitiveFaceFlux(
        state, connectivity, fusedResidualView, 1, & fusedFluxView );
    const double fusedFluxError = MaxDiff( cpuFlux, fusedFlux );
    const double fusedResidualError = MaxDiff( cpuResidual, fusedResidual );
    if ( fusedFluxError > 2.0e-12 || fusedResidualError > 2.0e-12 )
    {
        std::fprintf( stderr,
            "Main solver fused primitive path FAIL: flux %.3e residual %.3e\n",
            fusedFluxError, fusedResidualError );
        return false;
    }

    std::vector< ONEFLOW::Real > noTraceResidual =
        std::vector< ONEFLOW::Real >( nCells * nEq );
    for ( int i = 0; i < nCells * nEq; ++ i )
    {
        noTraceResidual[ i ] = 1.0e-6 * ( i % 7 );
    }
    ONEFLOW::ResidualView noTraceResidualView{
        nCells, nEq, noTraceResidual.data() };
    ONEFLOW::HipFluxBackend noTraceBackend;
    noTraceBackend.CalcAndAddPrimitiveFaceFlux(
        state, connectivity, noTraceResidualView, 1 );
    const double noTraceResidualError =
        MaxDiff( cpuResidual, noTraceResidual );
    if ( noTraceResidualError > 2.0e-12 )
    {
        std::fprintf( stderr,
            "Main solver fused NoTrace path FAIL: %.3e\n",
            noTraceResidualError );
        return false;
    }

    const double residualError = MaxDiff( cpuResidual, hipResidual );
    if ( residualError > 2.0e-12 )
    {
        std::fprintf( stderr,
            "Main solver explicit-mask residual FAIL: %.3e\n",
            residualError );
        return false;
    }

    std::vector< ONEFLOW::Real > expectedResidual( nCells * nEq );
    for ( int i = 0; i < nCells * nEq; ++ i )
    {
        expectedResidual[ i ] = 1.0e-6 * ( i % 7 );
    }
    for ( int face = 0; face < nFaces; ++ face )
    {
        for ( int eq = 0; eq < nEq; ++ eq )
        {
            const ONEFLOW::Real value = hipFlux[ eq * nFaces + face ];
            expectedResidual[ eq * nCells + leftCell[ face ] ] -= value;
            if ( boundaryMask[ face ] == 0 )
            {
                expectedResidual[ eq * nCells + rightCell[ face ] ] += value;
            }
        }
    }
    const double boundarySemanticError =
        MaxDiff( expectedResidual, hipResidual );
    if ( boundarySemanticError > 2.0e-12 )
    {
        std::fprintf( stderr,
            "Main solver boundary-mask semantics FAIL: %.3e\n",
            boundarySemanticError );
        return false;
    }

    double conservationError = 0.0;
    for ( int eq = 0; eq < nEq; ++ eq )
    {
        double residualDelta = 0.0;
        double boundaryFlux = 0.0;
        for ( int cell = 0; cell < nCells; ++ cell )
        {
            const int index = eq * nCells + cell;
            residualDelta += hipResidual[ index ]
                - 1.0e-6 * ( index % 7 );
        }
        for ( int face = 0; face < nFaces; ++ face )
        {
            if ( boundaryMask[ face ] != 0 )
            {
                boundaryFlux -= hipFlux[ eq * nFaces + face ];
            }
        }
        conservationError = std::max(
            conservationError, std::abs( residualDelta - boundaryFlux ) );
    }
    if ( conservationError > 2.0e-11 )
    {
        std::fprintf( stderr,
            "Main solver internal-face conservation FAIL: %.3e\n",
            conservationError );
        return false;
    }

    std::printf(
        "OneFLOW HIP main solver one-call: PASS "
        "(flux %.3e, residual %.3e, fused %.3e/%.3e, "
        "boundary %.3e, conservation %.3e, %d faces, %d eq)\n",
        fluxError, residualError, fusedFluxError, fusedResidualError,
        boundarySemanticError, conservationError,
        nFaces, nEq );
    return true;
}

bool TestMainSolverGradient()
{
    constexpr int nCells = 2;
    constexpr int nBoundaryFaces = 2;
    constexpr int nFaces = 3;
    constexpr int nEq = 5;
    constexpr int nTotalCells = nCells + nBoundaryFaces;
    const int leftCell[ nFaces ] = { 0, 1, 0 };
    const int rightCell[ nFaces ] = { 2, 3, 1 };
    const ONEFLOW::Real xCell[ nTotalCells ] = { 0.0, 1.0, -1.0, 2.0 };
    const ONEFLOW::Real yCell[ nTotalCells ] = { 0.0, 0.2, 0.0, 0.2 };
    const ONEFLOW::Real zCell[ nTotalCells ] = { 0.0, 0.1, 0.0, 0.1 };
    const ONEFLOW::Real xFace[ nFaces ] = { -0.5, 1.5, 0.5 };
    const ONEFLOW::Real yFace[ nFaces ] = { 0.0, 0.2, 0.1 };
    const ONEFLOW::Real zFace[ nFaces ] = { 0.0, 0.1, 0.05 };
    const ONEFLOW::Real xNormal[ nFaces ] = { -1.0, 1.0, 0.98 };
    const ONEFLOW::Real yNormal[ nFaces ] = { 0.0, 0.0, 0.2 };
    const ONEFLOW::Real zNormal[ nFaces ] = { 0.0, 0.0, 0.02 };
    const ONEFLOW::Real area[ nFaces ] = { 1.0, 1.1, 1.2 };
    const ONEFLOW::Real volume[ nCells ] = { 1.0, 1.2 };
    std::vector< ONEFLOW::Real > q( nEq * nTotalCells );
    for ( int eq = 0; eq < nEq; ++ eq )
        for ( int cell = 0; cell < nTotalCells; ++ cell )
            q[ eq * nTotalCells + cell ] =
                0.7 + 0.13 * eq + 0.07 * cell + 0.01 * eq * cell;

    std::vector< ONEFLOW::Real > expected( 3 * nEq * nTotalCells, 0.0 );
    for ( int face = 0; face < nFaces; ++ face )
    {
        const int left = leftCell[ face ];
        const int right = rightCell[ face ];
        const double dxl = xFace[ face ] - xCell[ left ];
        const double dyl = yFace[ face ] - yCell[ left ];
        const double dzl = zFace[ face ] - zCell[ left ];
        const double dxr = xFace[ face ] - xCell[ right ];
        const double dyr = yFace[ face ] - yCell[ right ];
        const double dzr = zFace[ face ] - zCell[ right ];
        const double delt1 = std::sqrt( dxl * dxl + dyl * dyl + dzl * dzl );
        const double delt2 = std::sqrt( dxr * dxr + dyr * dyr + dzr * dzr );
        const double delta = 1.0 / ( delt1 + delt2 + 1.0e-40 );
        const double cl = delt2 * delta;
        const double cr = delt1 * delta;
        const double ax = xNormal[ face ] * area[ face ];
        const double ay = yNormal[ face ] * area[ face ];
        const double az = zNormal[ face ] * area[ face ];
        for ( int eq = 0; eq < nEq; ++ eq )
        {
            const double value = cl * q[ eq * nTotalCells + left ]
                + cr * q[ eq * nTotalCells + right ];
            const int base = eq * 3 * nTotalCells;
            expected[ base + left ] += ax * value;
            expected[ base + nTotalCells + left ] += ay * value;
            expected[ base + 2 * nTotalCells + left ] += az * value;
            if ( face >= nBoundaryFaces )
            {
                expected[ base + right ] -= ax * value;
                expected[ base + nTotalCells + right ] -= ay * value;
                expected[ base + 2 * nTotalCells + right ] -= az * value;
            }
        }
    }
    for ( int eq = 0; eq < nEq; ++ eq )
    {
        const int base = eq * 3 * nTotalCells;
        for ( int cell = 0; cell < nCells; ++ cell )
        {
            expected[ base + cell ] /= volume[ cell ];
            expected[ base + nTotalCells + cell ] /= volume[ cell ];
            expected[ base + 2 * nTotalCells + cell ] /= volume[ cell ];
        }
        for ( int face = 0; face < nBoundaryFaces; ++ face )
        {
            const int left = leftCell[ face ];
            const int ghost = rightCell[ face ];
            expected[ base + ghost ] = expected[ base + left ];
            expected[ base + nTotalCells + ghost ] =
                expected[ base + nTotalCells + left ];
            expected[ base + 2 * nTotalCells + ghost ] =
                expected[ base + 2 * nTotalCells + left ];
        }
    }

    std::vector< ONEFLOW::Real > actual( 3 * nEq * nTotalCells );
    ONEFLOW::CellGradientView view;
    view.nCells = nCells;
    view.nGhostCells = nBoundaryFaces;
    view.nFaces = nFaces;
    view.nBoundaryFaces = nBoundaryFaces;
    view.nEquations = nEq;
    view.xFace = xFace; view.yFace = yFace; view.zFace = zFace;
    view.xNormal = xNormal; view.yNormal = yNormal;
    view.zNormal = zNormal; view.faceArea = area;
    view.xCell = xCell; view.yCell = yCell; view.zCell = zCell;
    view.cellVolume = volume;
    view.leftCell = leftCell; view.rightCell = rightCell;
    for ( int eq = 0; eq < nEq; ++ eq )
    {
        const int base = eq * 3 * nTotalCells;
        view.q[ eq ] = q.data() + eq * nTotalCells;
        view.dqdx[ eq ] = actual.data() + base;
        view.dqdy[ eq ] = actual.data() + base + nTotalCells;
        view.dqdz[ eq ] = actual.data() + base + 2 * nTotalCells;
    }
    ONEFLOW::HipFluxBackend backend;
    backend.CalcGradient( view );
    const double error = MaxDiff( expected, actual );
    if ( error > 1.0e-13 )
    {
        std::fprintf( stderr, "Main solver HIP gradient FAIL: %.3e\n", error );
        return false;
    }
    std::printf(
        "OneFLOW HIP main solver gradient: PASS (max error %.3e)\n",
        error );
    return true;
}

bool TestMainSolverAleAreaContract()
{
    constexpr int nFaces = 1;
    constexpr int nEq = 5;
    constexpr double gamma = 1.4;
    constexpr double density = 1.2;
    constexpr double u = 0.4;
    constexpr double v = -0.2;
    constexpr double w = 0.1;
    constexpr double pressure = 0.9;
    constexpr double nx = 0.6;
    constexpr double ny = 0.8;
    constexpr double nz = 0.0;
    constexpr double meshVelocityNormal = 0.15;
    constexpr double area = 2.5;
    const ONEFLOW::Real primitive[] = {
        density, u, v, w, pressure };
    const ONEFLOW::Real normalX[] = { nx };
    const ONEFLOW::Real normalY[] = { ny };
    const ONEFLOW::Real normalZ[] = { nz };
    const ONEFLOW::Real meshVelocity[] = { meshVelocityNormal };
    const ONEFLOW::Real faceArea[] = { area };
    ONEFLOW::Real fluxValues[ nEq ] = {};

    ONEFLOW::PrimitiveFaceStateView state;
    state.nFaces = nFaces;
    state.nEquations = nEq;
    state.primitiveLeft = primitive;
    state.primitiveRight = primitive;
    state.xNormal = normalX;
    state.yNormal = normalY;
    state.zNormal = normalZ;
    state.meshVelocityNormal = meshVelocity;
    state.faceArea = faceArea;
    state.gamma = gamma;

    ONEFLOW::FaceFluxView flux{ nFaces, nEq, fluxValues };
    ONEFLOW::HipFluxBackend hipBackend;
    ONEFLOW::EulerCpuAdapter adapter;
    adapter.CalcInvFlux( state, flux, hipBackend, 1 );

    const double relativeNormalVelocity =
        nx * u + ny * v + nz * w - meshVelocityNormal;
    const double massFlux = density * relativeNormalVelocity;
    const double totalEnergy = pressure / ( gamma - 1.0 )
        + 0.5 * density * ( u * u + v * v + w * w );
    const double totalEnthalpy = ( totalEnergy + pressure ) / density;
    const double expected[] = {
        area * massFlux,
        area * ( massFlux * u + nx * pressure ),
        area * ( massFlux * v + ny * pressure ),
        area * ( massFlux * w + nz * pressure ),
        area * ( massFlux * totalEnthalpy
            + meshVelocityNormal * pressure )
    };

    double maxError = 0.0;
    for ( int eq = 0; eq < nEq; ++ eq )
    {
        if ( ! std::isfinite( fluxValues[ eq ] ) )
        {
            std::fprintf( stderr,
                "Main solver ALE/area flux is non-finite at %d\n", eq );
            return false;
        }
        maxError = std::max(
            maxError,
            static_cast< double >(
                std::abs( fluxValues[ eq ] - expected[ eq ] ) ) );
    }
    if ( maxError > 1.0e-12 )
    {
        std::fprintf( stderr,
            "Main solver ALE/face-area ownership FAIL: %.3e\n", maxError );
        return false;
    }

    std::printf(
        "OneFLOW HIP main solver ALE/face-area contract: PASS "
        "(max error %.3e)\n",
        maxError );
    return true;
}

bool TestMainSolverBoundarySemantics()
{
    constexpr int nCells = 4;
    constexpr int nBoundaryFaces = 4;
    constexpr int nFaces = 5;
    constexpr int nEq = 5;
    constexpr int nTotalCells = nCells + nBoundaryFaces;
    constexpr double gamma = 1.4;

    const int leftCell[ nFaces ] = { 0, 1, 2, 3, 0 };
    const int rightCell[ nFaces ] = { 4, 5, 6, 7, 1 };
    const unsigned char boundaryMask[ nFaces ] = { 1, 1, 1, 1, 0 };
    const int boundaryType[ nBoundaryFaces ] = { -1, 99, 1, 2 };

    const ONEFLOW::Real xCell[ nTotalCells ] =
        { 0.0, 1.0, 2.0, 3.0, -1.0, -2.0, 4.0, 5.0 };
    const ONEFLOW::Real yCell[ nTotalCells ] =
        { 0.0, 0.1, 0.2, 0.3, 0.0, 0.1, 0.2, 0.3 };
    const ONEFLOW::Real zCell[ nTotalCells ] =
        { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    const ONEFLOW::Real xFace[ nFaces ] = { -0.5, -0.5, 1.5, 2.5, 0.5 };
    const ONEFLOW::Real yFace[ nFaces ] = { 0.0, 0.1, 0.15, 0.25, 0.05 };
    const ONEFLOW::Real zFace[ nFaces ] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
    const ONEFLOW::Real xNormal[ nFaces ] =
        { 1.0, 1.0, 1.0, 1.0, 1.0 };
    const ONEFLOW::Real yNormal[ nFaces ] = {};
    const ONEFLOW::Real zNormal[ nFaces ] = {};
    const ONEFLOW::Real faceArea[ nFaces ] =
        { 1.0, 1.0, 1.0, 1.0, 1.0 };
    const ONEFLOW::Real cellVolume[ nCells ] =
        { 1.0, 1.0, 1.0, 1.0 };
    const ONEFLOW::Real meshVelocity[ nFaces ] = {};

    std::vector< ONEFLOW::Real > primitive(
        nEq * nTotalCells, 0.0 );
    for ( int cell = 0; cell < nCells; ++ cell )
    {
        primitive[ 0 * nTotalCells + cell ] = 1.0 + 0.1 * cell;
        primitive[ 1 * nTotalCells + cell ] = 0.2 + 0.1 * cell;
        primitive[ 2 * nTotalCells + cell ] = 0.05 * cell;
        primitive[ 3 * nTotalCells + cell ] = -0.02 * cell;
        primitive[ 4 * nTotalCells + cell ] = 1.0 + 0.2 * cell;
    }
    const ONEFLOW::Real ghostSeed[ nBoundaryFaces ][ nEq ] = {
        { 7.0, 0.31, 0.12, -0.08, 3.0 },
        { 8.0, 0.41, 0.22, -0.18, 4.0 },
        { 9.0, 0.51, 0.32, -0.28, 5.0 },
        { 10.0, 0.61, 0.42, -0.38, 6.0 } };
    for ( int ghost = 0; ghost < nBoundaryFaces; ++ ghost )
    {
        for ( int equation = 0; equation < nEq; ++ equation )
        {
            primitive[ equation * nTotalCells + nCells + ghost ] =
                ghostSeed[ ghost ][ equation ];
        }
    }

    const ONEFLOW::EulerDomainProblem problem{
        nCells, nBoundaryFaces, nEq, gamma, 1.0e-3, 1.0,
        ONEFLOW::EulerDomainBoundary::Periodic };
    const ONEFLOW::EulerDomainStateKey key{
        0, 0, 0, ONEFLOW::AccelBackendKind::HIP, 0 };
    ONEFLOW::Ns3DDeviceState owner( problem, key );
    owner.ReserveFaces( nFaces );
    owner.AttachBackendState( ONEFLOW::CreateHipNs3DBackendState( 0 ) );

    ONEFLOW::HipFluxBackend backend;
    const void * cacheKey = &owner;
    backend.BindState( cacheKey, owner );

    std::vector< ONEFLOW::Real > gradient(
        nEq * 3 * nTotalCells, 0.0 );
    ONEFLOW::CellGradientView gradientView;
    gradientView.nCells = nCells;
    gradientView.nGhostCells = nBoundaryFaces;
    gradientView.nFaces = nFaces;
    gradientView.nBoundaryFaces = nBoundaryFaces;
    gradientView.nEquations = nEq;
    gradientView.xFace = xFace;
    gradientView.yFace = yFace;
    gradientView.zFace = zFace;
    gradientView.xNormal = xNormal;
    gradientView.yNormal = yNormal;
    gradientView.zNormal = zNormal;
    gradientView.faceArea = faceArea;
    gradientView.xCell = xCell;
    gradientView.yCell = yCell;
    gradientView.zCell = zCell;
    gradientView.cellVolume = cellVolume;
    gradientView.leftCell = leftCell;
    gradientView.rightCell = rightCell;
    gradientView.cacheKey = cacheKey;
    for ( int equation = 0; equation < nEq; ++ equation )
    {
        gradientView.q[ equation ] =
            primitive.data() + equation * nTotalCells;
        const int base = equation * 3 * nTotalCells;
        gradientView.dqdx[ equation ] = gradient.data() + base;
        gradientView.dqdy[ equation ] =
            gradient.data() + base + nTotalCells;
        gradientView.dqdz[ equation ] =
            gradient.data() + base + 2 * nTotalCells;
    }
    backend.CalcGradient( gradientView );

    const ONEFLOW::Real inflow[ nEq ] =
        { 1.0, 0.0, 0.0, 0.0, 1.0 };
    ONEFLOW::BoundaryStateView boundaryView;
    boundaryView.nCells = nCells;
    boundaryView.nGhostCells = nBoundaryFaces;
    boundaryView.nBoundaryFaces = nBoundaryFaces;
    boundaryView.nEquations = nEq;
    boundaryView.leftCell = leftCell;
    boundaryView.rightCell = rightCell;
    boundaryView.boundaryType = boundaryType;
    boundaryView.xNormal = xNormal;
    boundaryView.yNormal = yNormal;
    boundaryView.zNormal = zNormal;
    boundaryView.meshVelocityNormal = meshVelocity;
    boundaryView.inflow = inflow;
    boundaryView.gamma = gamma;
    boundaryView.cacheKey = cacheKey;
    boundaryView.preferResidentState = true;
    for ( int equation = 0; equation < nEq; ++ equation )
        boundaryView.q[ equation ] =
            primitive.data() + equation * nTotalCells;
    backend.UpdateBoundaryState( boundaryView );

    double ghostError = 0.0;
    double boundaryBcQ[ nEq ] = {};
    for ( int equation = 0; equation < nEq; ++ equation )
    {
        ghostError = std::max(
            ghostError,
            static_cast< double >( std::abs(
                primitive[ equation * nTotalCells + 4 ]
                - ghostSeed[ 0 ][ equation ] ) ) );
        ghostError = std::max(
            ghostError,
            static_cast< double >( std::abs(
                primitive[ equation * nTotalCells + 5 ]
                - ghostSeed[ 1 ][ equation ] ) ) );
        ghostError = std::max(
            ghostError,
            static_cast< double >( std::abs(
                primitive[ equation * nTotalCells + 6 ]
                - primitive[ equation * nTotalCells + 2 ] ) ) );
    }

    const double solidDensity =
        primitive[ 0 * nTotalCells + 3 ];
    const double solidU = -primitive[ 1 * nTotalCells + 3 ];
    const double solidV = primitive[ 2 * nTotalCells + 3 ];
    const double solidW = primitive[ 3 * nTotalCells + 3 ];
    const double solidPressure =
        primitive[ 4 * nTotalCells + 3 ];
    ghostError = std::max(
        ghostError,
        static_cast< double >( std::abs(
            primitive[ 0 * nTotalCells + 7 ] - solidDensity ) ) );
    ghostError = std::max(
        ghostError,
        static_cast< double >( std::abs(
            primitive[ 1 * nTotalCells + 7 ] - solidU ) ) );
    ghostError = std::max(
        ghostError,
        static_cast< double >( std::abs(
            primitive[ 2 * nTotalCells + 7 ] - solidV ) ) );
    ghostError = std::max(
        ghostError,
        static_cast< double >( std::abs(
            primitive[ 3 * nTotalCells + 7 ] - solidW ) ) );
    ghostError = std::max(
        ghostError,
        static_cast< double >( std::abs(
            primitive[ 4 * nTotalCells + 7 ] - solidPressure ) ) );

    for ( int equation = 0; equation < nEq; ++ equation )
    {
        boundaryBcQ[ equation ] = 0.5 * (
            primitive[ equation * nTotalCells + 3 ]
            + primitive[ equation * nTotalCells + 7 ] );
    }

    std::vector< ONEFLOW::Real > qLeft( nEq * nFaces, 0.0 );
    std::vector< ONEFLOW::Real > qRight( nEq * nFaces, 0.0 );
    std::vector< ONEFLOW::Real > reconstructionBcQ(
        nEq * nFaces, 0.0 );
    for ( int equation = 0; equation < nEq; ++ equation )
    {
        reconstructionBcQ[ equation * nFaces + 3 ] =
            boundaryBcQ[ equation ] + 3.0 + equation;
    }
    ONEFLOW::ReconstructionBoundaryOperation operations[ nFaces ] = {
        ONEFLOW::ReconstructionBoundaryOperation::Preserve,
        ONEFLOW::ReconstructionBoundaryOperation::Preserve,
        ONEFLOW::ReconstructionBoundaryOperation::Average,
        ONEFLOW::ReconstructionBoundaryOperation::SolidOverride,
        ONEFLOW::ReconstructionBoundaryOperation::Preserve };
    const ONEFLOW::BoundaryGhostUpdatePolicy policies[ nFaces ] = {
        ONEFLOW::BoundaryGhostUpdatePolicy::ExternalHalo,
        ONEFLOW::BoundaryGhostUpdatePolicy::ExternalHalo,
        ONEFLOW::BoundaryGhostUpdatePolicy::LocalBoundaryModel,
        ONEFLOW::BoundaryGhostUpdatePolicy::LocalBoundaryModel,
        ONEFLOW::BoundaryGhostUpdatePolicy::ExternalHalo };

    ONEFLOW::CellFaceReconstructionView reconstructionView;
    reconstructionView.nCells = nCells;
    reconstructionView.nGhostCells = nBoundaryFaces;
    reconstructionView.nFaces = nFaces;
    reconstructionView.nBoundaryFaces = nBoundaryFaces;
    reconstructionView.nEquations = nEq;
    reconstructionView.xFace = xFace;
    reconstructionView.yFace = yFace;
    reconstructionView.zFace = zFace;
    reconstructionView.xCell = xCell;
    reconstructionView.yCell = yCell;
    reconstructionView.zCell = zCell;
    reconstructionView.leftCell = leftCell;
    reconstructionView.rightCell = rightCell;
    reconstructionView.boundaryMask = boundaryMask;
    reconstructionView.boundaryOperation = operations;
    reconstructionView.boundaryGhostPolicy = policies;
    reconstructionView.hasSolidBoundary = true;
    reconstructionView.ownerToken = cacheKey;
    reconstructionView.cacheKey = cacheKey;
    reconstructionView.topologyGeneration = owner.TopologyGeneration();
    reconstructionView.fieldGeneration = owner.FieldGeneration();
    for ( int equation = 0; equation < nEq; ++ equation )
    {
        reconstructionView.q[ equation ] =
            primitive.data() + equation * nTotalCells;
        const int base = equation * 3 * nTotalCells;
        reconstructionView.dqdx[ equation ] = gradient.data() + base;
        reconstructionView.dqdy[ equation ] =
            gradient.data() + base + nTotalCells;
        reconstructionView.dqdz[ equation ] =
            gradient.data() + base + 2 * nTotalCells;
        reconstructionView.bcQ[ equation ] =
            reconstructionBcQ.data() + equation * nFaces;
        reconstructionView.qLeft[ equation ] =
            qLeft.data() + equation * nFaces;
        reconstructionView.qRight[ equation ] =
            qRight.data() + equation * nFaces;
    }
    backend.ReconstructFaceValues( reconstructionView, true );

    double reconstructionError = 0.0;
    for ( int equation = 0; equation < nEq; ++ equation )
    {
        const double ordinaryAverage = 0.5 * (
            primitive[ equation * nTotalCells + 2 ]
            + primitive[ equation * nTotalCells + 6 ] );
        reconstructionError = std::max(
            reconstructionError,
            static_cast< double >( std::abs(
                qLeft[ equation * nFaces + 2 ] - ordinaryAverage ) ) );
        reconstructionError = std::max(
            reconstructionError,
            static_cast< double >( std::abs(
                qRight[ equation * nFaces + 2 ] - ordinaryAverage ) ) );
        reconstructionError = std::max(
            reconstructionError,
            static_cast< double >( std::abs(
                qLeft[ equation * nFaces + 3 ]
                - reconstructionBcQ[ equation * nFaces + 3 ] ) ) );
        reconstructionError = std::max(
            reconstructionError,
            static_cast< double >( std::abs(
                qRight[ equation * nFaces + 3 ]
                - reconstructionBcQ[ equation * nFaces + 3 ] ) ) );
    }
    const double preservedGhostError = std::max(
        static_cast< double >( std::abs(
            primitive[ 0 * nTotalCells + 4 ] - ghostSeed[ 0 ][ 0 ] ) ),
        static_cast< double >( std::abs(
            primitive[ 0 * nTotalCells + 5 ] - ghostSeed[ 1 ][ 0 ] ) ) );

    operations[ 3 ] = ONEFLOW::ReconstructionBoundaryOperation::Average;
    owner.MarkTopologyChanged();
    reconstructionView.topologyGeneration = owner.TopologyGeneration();
    backend.ReconstructFaceValues( reconstructionView, true );
    const double expectedAverage = 0.5 * (
        primitive[ 0 * nTotalCells + 3 ]
        + primitive[ 0 * nTotalCells + 7 ] );
    const double generationError = std::max(
        static_cast< double >( std::abs(
            qLeft[ 0 * nFaces + 3 ] - expectedAverage ) ),
        static_cast< double >( std::abs(
            qRight[ 0 * nFaces + 3 ] - expectedAverage ) ) );
    const double staleOperationError = std::abs(
        qLeft[ 0 * nFaces + 3 ] - reconstructionBcQ[ 0 * nFaces + 3 ] );

    backend.UnbindState( cacheKey );
    if ( ghostError > 2.0e-12
         || reconstructionError > 2.0e-12
         || preservedGhostError > 2.0e-12
         || generationError > 2.0e-12
         || staleOperationError < 1.0e-6 )
    {
        std::fprintf(
            stderr,
            "Main solver HIP boundary contract FAIL: ghost %.3e, "
            "reconstruction %.3e, preserved %.3e, generation %.3e, "
            "stale-op %.3e\n",
            ghostError, reconstructionError, preservedGhostError,
            generationError, staleOperationError );
        return false;
    }
    std::printf(
        "OneFLOW HIP main solver boundary contract: PASS "
        "(ghost %.3e, reconstruction %.3e, generation %.3e)\n",
        ghostError, reconstructionError, generationError );
    return true;
}

bool TestMainSolverStateUpdate()
{
    constexpr int nCells = 2;
    constexpr int nGhostCells = 2;
    constexpr int nTotalCells = nCells + nGhostCells;
    constexpr int nFaces = 3;
    constexpr int nBoundaryFaces = 2;
    constexpr int nEq = 5;
    constexpr double gamma = 1.4;
    constexpr double rkCoefficient = 0.2;

    const int leftCell[ nFaces ] = { 0, 1, 0 };
    const int rightCell[ nFaces ] = { 2, 3, 1 };
    const unsigned char boundaryMask[ nFaces ] = { 1, 1, 0 };
    const ONEFLOW::ReconstructionBoundaryOperation boundaryOperation[ nFaces ] = {
        ONEFLOW::ReconstructionBoundaryOperation::Preserve,
        ONEFLOW::ReconstructionBoundaryOperation::Preserve,
        ONEFLOW::ReconstructionBoundaryOperation::Preserve };
    const ONEFLOW::BoundaryGhostUpdatePolicy boundaryGhostPolicy[ nFaces ] = {
        ONEFLOW::BoundaryGhostUpdatePolicy::ExternalHalo,
        ONEFLOW::BoundaryGhostUpdatePolicy::ExternalHalo,
        ONEFLOW::BoundaryGhostUpdatePolicy::ExternalHalo };
    const ONEFLOW::Real xCell[ nTotalCells ] = { 0.0, 1.0, -1.0, 2.0 };
    const ONEFLOW::Real yCell[ nTotalCells ] = { 0.0, 0.2, 0.0, 0.2 };
    const ONEFLOW::Real zCell[ nTotalCells ] = { 0.0, 0.1, 0.0, 0.1 };
    const ONEFLOW::Real xFace[ nFaces ] = { -0.5, 1.5, 0.5 };
    const ONEFLOW::Real yFace[ nFaces ] = { 0.0, 0.2, 0.1 };
    const ONEFLOW::Real zFace[ nFaces ] = { 0.0, 0.1, 0.05 };
    const ONEFLOW::Real xNormal[ nFaces ] = { -1.0, 1.0, 0.98 };
    const ONEFLOW::Real yNormal[ nFaces ] = { 0.0, 0.0, 0.2 };
    const ONEFLOW::Real zNormal[ nFaces ] = { 0.0, 0.0, 0.02 };
    const ONEFLOW::Real faceMeshVelocity[ nFaces ] = { 0.0, 0.0, 0.0 };
    const ONEFLOW::Real faceArea[ nFaces ] = { 1.0, 1.1, 1.2 };
    const ONEFLOW::Real cellVolume[ nCells ] = { 1.0, 1.2 };
    const ONEFLOW::Real timeStep[ nCells ] = { 1.0e-3, 1.3e-3 };

    std::vector< ONEFLOW::Real > primitive(
        nEq * nTotalCells, 0.0 );
    for ( int cell = 0; cell < nCells; ++ cell )
    {
        primitive[ 0 * nTotalCells + cell ] = 1.0 + 0.04 * cell;
        primitive[ 1 * nTotalCells + cell ] = 0.2 + 0.03 * cell;
        primitive[ 2 * nTotalCells + cell ] = -0.1 + 0.02 * cell;
        primitive[ 3 * nTotalCells + cell ] = 0.05 - 0.01 * cell;
        primitive[ 4 * nTotalCells + cell ] = 1.0 + 0.03 * cell;
    }
    // Boundary ghost values start as copies of their owner cells, matching
    // the restricted ghost-copy semantics used by the production seam.
    for ( int equation = 0; equation < nEq; ++ equation )
    {
        primitive[ equation * nTotalCells + 2 ] =
            primitive[ equation * nTotalCells + 0 ];
        primitive[ equation * nTotalCells + 3 ] =
            primitive[ equation * nTotalCells + 1 ];
    }

    const ONEFLOW::EulerDomainProblem problem{
        nCells, nGhostCells, nEq, gamma, 1.0e-3, 1.0,
        ONEFLOW::EulerDomainBoundary::Periodic };
    const ONEFLOW::EulerDomainStateKey key{
        0, 0, 0, ONEFLOW::AccelBackendKind::HIP, 0 };
    ONEFLOW::Ns3DDeviceState owner( problem, key );
    owner.ReserveFaces( nFaces );
    owner.AttachBackendState( ONEFLOW::CreateHipNs3DBackendState( 0 ) );

    ONEFLOW::HipFluxBackend backend;
    const void * cacheKey = &owner;
    backend.BindState( cacheKey, owner );

    std::vector< ONEFLOW::Real > gradient(
        nEq * 3 * nTotalCells, 0.0 );
    ONEFLOW::CellGradientView gradientView;
    gradientView.nCells = nCells;
    gradientView.nGhostCells = nGhostCells;
    gradientView.nFaces = nFaces;
    gradientView.nBoundaryFaces = nBoundaryFaces;
    gradientView.nEquations = nEq;
    gradientView.xFace = xFace;
    gradientView.yFace = yFace;
    gradientView.zFace = zFace;
    gradientView.xNormal = xNormal;
    gradientView.yNormal = yNormal;
    gradientView.zNormal = zNormal;
    gradientView.faceArea = faceArea;
    gradientView.xCell = xCell;
    gradientView.yCell = yCell;
    gradientView.zCell = zCell;
    gradientView.cellVolume = cellVolume;
    gradientView.leftCell = leftCell;
    gradientView.rightCell = rightCell;
    gradientView.cacheKey = cacheKey;
    for ( int equation = 0; equation < nEq; ++ equation )
    {
        const int base = equation * 3 * nTotalCells;
        gradientView.q[ equation ] =
            primitive.data() + equation * nTotalCells;
        gradientView.dqdx[ equation ] = gradient.data() + base;
        gradientView.dqdy[ equation ] = gradient.data() + base + nTotalCells;
        gradientView.dqdz[ equation ] =
            gradient.data() + base + 2 * nTotalCells;
    }
    backend.CalcGradient( gradientView );

    std::vector< ONEFLOW::Real > qLeft( nEq * nFaces, 0.0 );
    std::vector< ONEFLOW::Real > qRight( nEq * nFaces, 0.0 );
    ONEFLOW::CellFaceReconstructionView reconstructionView;
    reconstructionView.nCells = nCells;
    reconstructionView.nGhostCells = nGhostCells;
    reconstructionView.nFaces = nFaces;
    reconstructionView.nBoundaryFaces = nBoundaryFaces;
    reconstructionView.nEquations = nEq;
    reconstructionView.xFace = xFace;
    reconstructionView.yFace = yFace;
    reconstructionView.zFace = zFace;
    reconstructionView.xCell = xCell;
    reconstructionView.yCell = yCell;
    reconstructionView.zCell = zCell;
    reconstructionView.leftCell = leftCell;
    reconstructionView.rightCell = rightCell;
    reconstructionView.boundaryMask = boundaryMask;
    reconstructionView.boundaryOperation = boundaryOperation;
    reconstructionView.boundaryGhostPolicy = boundaryGhostPolicy;
    reconstructionView.memorySpace = ONEFLOW::MemorySpace::Host;
    reconstructionView.limiterMode =
        ONEFLOW::ReconstructionLimiterMode::Disabled;
    reconstructionView.physicality =
        ONEFLOW::ReconstructionPhysicalityPolicy::PositiveDensityPressure;
    reconstructionView.ownerToken = cacheKey;
    reconstructionView.cacheKey = cacheKey;
    reconstructionView.topologyGeneration = owner.TopologyGeneration();
    reconstructionView.fieldGeneration = owner.FieldGeneration();
    for ( int equation = 0; equation < nEq; ++ equation )
    {
        reconstructionView.q[ equation ] =
            primitive.data() + equation * nTotalCells;
        const int base = equation * 3 * nTotalCells;
        reconstructionView.dqdx[ equation ] = gradient.data() + base;
        reconstructionView.dqdy[ equation ] =
            gradient.data() + base + nTotalCells;
        reconstructionView.dqdz[ equation ] =
            gradient.data() + base + 2 * nTotalCells;
        reconstructionView.qLeft[ equation ] =
            qLeft.data() + equation * nFaces;
        reconstructionView.qRight[ equation ] =
            qRight.data() + equation * nFaces;
    }
    backend.ReconstructFaceValues( reconstructionView, false );

    ONEFLOW::FaceConnectivityView connectivity;
    connectivity.nFaces = nFaces;
    connectivity.nBoundaryFaces = nBoundaryFaces;
    connectivity.leftCell = leftCell;
    connectivity.rightCell = rightCell;
    connectivity.boundaryMask = boundaryMask;
    std::vector< ONEFLOW::Real > residual( nCells * nEq, 0.0 );
    ONEFLOW::ResidualView residualView{
        nCells, nEq, residual.data() };
    backend.CalcAndAddReconstructedFaceFlux(
        connectivity, residualView, 1, gamma, nullptr, cacheKey, true );

    std::vector< ONEFLOW::Real > expected( nEq * nCells, 0.0 );
    for ( int cell = 0; cell < nCells; ++ cell )
    {
        const double density = primitive[ 0 * nTotalCells + cell ];
        const double u = primitive[ 1 * nTotalCells + cell ];
        const double v = primitive[ 2 * nTotalCells + cell ];
        const double w = primitive[ 3 * nTotalCells + cell ];
        const double pressure = primitive[ 4 * nTotalCells + cell ];
        double conserved[ nEq ] = {
            density,
            density * u,
            density * v,
            density * w,
            pressure / ( gamma - 1.0 )
                + 0.5 * density * ( u * u + v * v + w * w ) };
        const double scale =
            timeStep[ cell ] / cellVolume[ cell ] * rkCoefficient;
        for ( int equation = 0; equation < nEq; ++ equation )
        {
            expected[ equation * nCells + cell ] =
                residual[ equation * nCells + cell ] * scale;
            conserved[ equation ] +=
                expected[ equation * nCells + cell ];
        }
        const double updatedDensity = conserved[ 0 ];
        const double updatedU = conserved[ 1 ] / updatedDensity;
        const double updatedV = conserved[ 2 ] / updatedDensity;
        const double updatedW = conserved[ 3 ] / updatedDensity;
        const double updatedPressure = ( gamma - 1.0 )
            * ( conserved[ 4 ]
                - 0.5 * updatedDensity
                    * ( updatedU * updatedU
                        + updatedV * updatedV
                        + updatedW * updatedW ) );
        expected[ 0 * nCells + cell ] = updatedDensity;
        expected[ 1 * nCells + cell ] = updatedU;
        expected[ 2 * nCells + cell ] = updatedV;
        expected[ 3 * nCells + cell ] = updatedW;
        expected[ 4 * nCells + cell ] = updatedPressure;
    }

    ONEFLOW::CellStateUpdateView updateView;
    updateView.nCells = nCells;
    updateView.nGhostCells = nGhostCells;
    updateView.nEquations = nEq;
    updateView.timeStep = timeStep;
    updateView.cellVolume = cellVolume;
    for ( int equation = 0; equation < nEq; ++ equation )
    {
        updateView.primitive[ equation ] =
            primitive.data() + equation * nTotalCells;
    }
    updateView.gamma = gamma;
    updateView.rkCoefficient = rkCoefficient;
    updateView.cacheKey = cacheKey;
    updateView.retainResidentState = true;

    backend.ScaleCurrentResidual( updateView );
    backend.UpdatePrimitiveState( updateView, true );

    double error = 0.0;
    for ( int equation = 0; equation < nEq; ++ equation )
        for ( int cell = 0; cell < nCells; ++ cell )
            error = std::max(
                error,
                static_cast< double >(
                    std::abs(
                        primitive[ equation * nTotalCells + cell ]
                        - expected[ equation * nCells + cell ] ) ) );
    if ( error > 2.0e-12 )
    {
        std::fprintf(
            stderr, "Main solver HIP state update FAIL: %.3e\n", error );
        return false;
    }

    // Prove that the boundary update consumes the state retained by the
    // preceding device update.  Corrupting the host internal values would
    // contaminate the ghosts if UpdateBoundaryState uploaded q again.
    for ( int equation = 0; equation < nEq; ++ equation )
    {
        primitive[ equation * nTotalCells + 0 ] = 10.0 + equation;
        primitive[ equation * nTotalCells + 1 ] = 20.0 + equation;
    }
    const int boundaryType[ nBoundaryFaces ] = { 1, 1 };
    const ONEFLOW::Real boundaryMeshVelocity[ nBoundaryFaces ] = {
        0.0, 0.0 };
    const ONEFLOW::Real inflow[ nEq ] = {
        1.0, 0.0, 0.0, 0.0, 1.0 };
    ONEFLOW::BoundaryStateView boundaryView;
    boundaryView.nCells = nCells;
    boundaryView.nGhostCells = nGhostCells;
    boundaryView.nBoundaryFaces = nBoundaryFaces;
    boundaryView.nEquations = nEq;
    boundaryView.leftCell = leftCell;
    boundaryView.rightCell = rightCell;
    boundaryView.boundaryType = boundaryType;
    boundaryView.xNormal = xNormal;
    boundaryView.yNormal = yNormal;
    boundaryView.zNormal = zNormal;
    boundaryView.meshVelocityNormal = boundaryMeshVelocity;
    boundaryView.inflow = inflow;
    boundaryView.gamma = gamma;
    boundaryView.cacheKey = cacheKey;
    boundaryView.preferResidentState = true;
    for ( int equation = 0; equation < nEq; ++ equation )
    {
        boundaryView.q[ equation ] =
            primitive.data() + equation * nTotalCells;
    }
    backend.UpdateBoundaryState( boundaryView );

    double residentBoundaryError = 0.0;
    for ( int equation = 0; equation < nEq; ++ equation )
    {
        for ( int face = 0; face < nBoundaryFaces; ++ face )
        {
            residentBoundaryError = std::max(
                residentBoundaryError,
                static_cast< double >( std::abs(
                    primitive[
                        equation * nTotalCells + rightCell[ face ] ]
                    - expected[
                        equation * nCells + leftCell[ face ] ] ) ) );
        }
    }
    if ( residentBoundaryError > 2.0e-12 )
    {
        std::fprintf(
            stderr,
            "Main solver HIP resident boundary handoff FAIL: %.3e\n",
            residentBoundaryError );
        return false;
    }

    constexpr double cfl = 0.25;
    constexpr double maxTimeRatio = 1.1;
    ONEFLOW::Real deviceSpectrum[ nCells ] = {};
    ONEFLOW::Real deviceTimeStep[ nCells ] = {};
    ONEFLOW::Real expectedSpectrum[ nCells ] = {};
    ONEFLOW::Real expectedTimeStep[ nCells ] = {};
    const auto expectedPrimitive = [&]( int equation, int cell ) {
        const int internal = cell < nCells ? cell : cell - nCells;
        return expected[ equation * nCells + internal ];
    };
    for ( int face = 0; face < nFaces; ++ face )
    {
        const int left = leftCell[ face ];
        const int right = rightCell[ face ];
        const double rhoLeft = expectedPrimitive( 0, left );
        const double rhoRight = expectedPrimitive( 0, right );
        const double pressureLeft = expectedPrimitive( 4, left );
        const double pressureRight = expectedPrimitive( 4, right );
        const double vnLeft =
            xNormal[ face ] * expectedPrimitive( 1, left )
            + yNormal[ face ] * expectedPrimitive( 2, left )
            + zNormal[ face ] * expectedPrimitive( 3, left )
            - faceMeshVelocity[ face ];
        const double vnRight =
            xNormal[ face ] * expectedPrimitive( 1, right )
            + yNormal[ face ] * expectedPrimitive( 2, right )
            + zNormal[ face ] * expectedPrimitive( 3, right )
            - faceMeshVelocity[ face ];
        const double soundSpeed = std::sqrt(
            gamma * 0.5 * ( pressureLeft + pressureRight )
            / ( 0.5 * ( rhoLeft + rhoRight ) ) );
        const double faceSpectrum = 0.5 * faceArea[ face ]
            * ( std::abs( 0.5 * ( vnLeft + vnRight ) )
                + soundSpeed );
        expectedSpectrum[ left ] += faceSpectrum;
        if ( right < nCells ) expectedSpectrum[ right ] += faceSpectrum;
    }
    double minimumTimeStep = std::numeric_limits< double >::max();
    for ( int cell = 0; cell < nCells; ++ cell )
    {
        expectedTimeStep[ cell ] =
            cfl * cellVolume[ cell ] / expectedSpectrum[ cell ];
        minimumTimeStep =
            std::min( minimumTimeStep, expectedTimeStep[ cell ] );
    }
    for ( int cell = 0; cell < nCells; ++ cell )
    {
        expectedTimeStep[ cell ] = std::min(
            expectedTimeStep[ cell ],
            maxTimeRatio * minimumTimeStep );
    }

    ONEFLOW::CellTimeStepView timeStepView;
    timeStepView.nCells = nCells;
    timeStepView.nGhostCells = nGhostCells;
    timeStepView.nFaces = nFaces;
    timeStepView.nBoundaryFaces = nBoundaryFaces;
    timeStepView.nEquations = nEq;
    timeStepView.xNormal = xNormal;
    timeStepView.yNormal = yNormal;
    timeStepView.zNormal = zNormal;
    timeStepView.meshVelocityNormal = faceMeshVelocity;
    timeStepView.faceArea = faceArea;
    timeStepView.cellVolume = cellVolume;
    timeStepView.leftCell = leftCell;
    timeStepView.rightCell = rightCell;
    timeStepView.inviscidSpectrum = deviceSpectrum;
    timeStepView.timeStep = deviceTimeStep;
    timeStepView.gamma = gamma;
    timeStepView.cfl = cfl;
    timeStepView.maxTimeRatio = maxTimeRatio;
    timeStepView.cacheKey = cacheKey;
    timeStepView.preferResidentState = true;
    for ( int equation = 0; equation < nEq; ++ equation )
    {
        timeStepView.q[ equation ] =
            primitive.data() + equation * nTotalCells;
    }
    backend.CalcLocalInviscidTimeStep( timeStepView );
    backend.UnbindState( cacheKey );

    double residentTimeStepError = 0.0;
    for ( int cell = 0; cell < nCells; ++ cell )
    {
        residentTimeStepError = std::max(
            residentTimeStepError,
            static_cast< double >( std::abs(
                deviceSpectrum[ cell ] - expectedSpectrum[ cell ] ) ) );
        residentTimeStepError = std::max(
            residentTimeStepError,
            static_cast< double >( std::abs(
                deviceTimeStep[ cell ] - expectedTimeStep[ cell ] ) ) );
    }
    if ( residentTimeStepError > 2.0e-12 )
    {
        std::fprintf(
            stderr,
            "Main solver HIP resident timestep FAIL: %.3e\n",
            residentTimeStepError );
        return false;
    }
    std::printf(
        "OneFLOW HIP main solver state update: PASS "
        "(max error %.3e, %d cells, %d equations)\n",
        error, nCells, nEq );
    std::printf(
        "OneFLOW HIP main solver state/boundary residency: PASS "
        "(state %.3e, boundary %.3e, %d cells, %d equations)\n",
        error, residentBoundaryError, nCells, nEq );
    std::printf(
        "OneFLOW HIP main solver resident timestep: PASS "
        "(max error %.3e, %d cells)\n",
        residentTimeStepError, nCells );
    return true;
}

bool TestMainSolverBoundaryExchange()
{
    constexpr int nCells = 2;
    constexpr int nGhostCells = 2;
    constexpr int totalCells = nCells + nGhostCells;
    constexpr int nFaces = 2;
    constexpr int nEq = 5;

    std::vector< ONEFLOW::Real > primitive( nEq * totalCells, -1.0 );
    for ( int equation = 0; equation < nEq; ++ equation )
    {
        primitive[ equation * totalCells + 0 ] = 10.0 * equation + 1.0;
        primitive[ equation * totalCells + 1 ] = 10.0 * equation + 2.0;
    }
    std::vector< ONEFLOW::Real > gradient(
        nEq * 3 * totalCells, 0.0 );
    const ONEFLOW::Real xFace[ nFaces ] = { -0.5, 1.5 };
    const ONEFLOW::Real yFace[ nFaces ] = { 0.0, 0.0 };
    const ONEFLOW::Real zFace[ nFaces ] = { 0.0, 0.0 };
    const ONEFLOW::Real xNormal[ nFaces ] = { -1.0, 1.0 };
    const ONEFLOW::Real yNormal[ nFaces ] = { 0.0, 0.0 };
    const ONEFLOW::Real zNormal[ nFaces ] = { 0.0, 0.0 };
    const ONEFLOW::Real faceArea[ nFaces ] = { 1.0, 1.0 };
    const ONEFLOW::Real xCell[ totalCells ] = { 0.0, 1.0, -1.0, 2.0 };
    const ONEFLOW::Real yCell[ totalCells ] = {};
    const ONEFLOW::Real zCell[ totalCells ] = {};
    const ONEFLOW::Real cellVolume[ nCells ] = { 1.0, 1.0 };
    const int leftCell[ nFaces ] = { 0, 1 };
    const int rightCell[ nFaces ] = { 2, 3 };

    const ONEFLOW::EulerDomainProblem problem{
        nCells, nGhostCells, nEq, 1.4, 1.0e-3, 1.0,
        ONEFLOW::EulerDomainBoundary::Periodic };
    const ONEFLOW::EulerDomainStateKey key{
        0, 0, 0, ONEFLOW::AccelBackendKind::HIP, 0 };
    ONEFLOW::Ns3DDeviceState owner( problem, key );
    owner.ReserveFaces( nFaces );
    owner.AttachBackendState( ONEFLOW::CreateHipNs3DBackendState( 0 ) );

    ONEFLOW::HipFluxBackend backend;
    const void * cacheKey = &owner;
    backend.BindState( cacheKey, owner );

    ONEFLOW::CellGradientView view;
    view.nCells = nCells;
    view.nGhostCells = nGhostCells;
    view.nFaces = nFaces;
    view.nBoundaryFaces = nFaces;
    view.nEquations = nEq;
    view.xFace = xFace;
    view.yFace = yFace;
    view.zFace = zFace;
    view.xNormal = xNormal;
    view.yNormal = yNormal;
    view.zNormal = zNormal;
    view.faceArea = faceArea;
    view.xCell = xCell;
    view.yCell = yCell;
    view.zCell = zCell;
    view.cellVolume = cellVolume;
    view.leftCell = leftCell;
    view.rightCell = rightCell;
    view.cacheKey = cacheKey;
    view.downloadToHost = false;
    for ( int equation = 0; equation < nEq; ++ equation )
    {
        view.q[ equation ] =
            primitive.data() + equation * totalCells;
        const int base = equation * 3 * totalCells;
        view.dqdx[ equation ] = gradient.data() + base;
        view.dqdy[ equation ] = gradient.data() + base + totalCells;
        view.dqdz[ equation ] = gradient.data() + base + 2 * totalCells;
    }
    backend.CalcGradient( view );

    ONEFLOW::BoundaryFieldExchangePlan plan;
    plan.sourceZone = 0;
    plan.targetZone = 1;
    plan.sourceRank = 0;
    plan.targetRank = 1;
    plan.nEquations = nEq;
    plan.sourceCellCount = totalCells;
    plan.targetCellCount = totalCells;
    plan.sendCellIds = { 1, 0 };
    plan.receiveCellIds = { 2, 3 };

    const ONEFLOW::BoundaryBufferView send = backend.PackBoundary(
        cacheKey, plan, ONEFLOW::BoundaryTransferMode::HostStaged );
    ONEFLOW::BoundaryBufferView receive = backend.PrepareBoundaryReceive(
        cacheKey, plan, ONEFLOW::BoundaryTransferMode::HostStaged );
    std::copy(
        send.values,
        send.values + send.valueCount,
        receive.values );
    backend.UnpackBoundary(
        cacheKey, plan, ONEFLOW::BoundaryTransferMode::HostStaged );

    ONEFLOW::BoundaryFieldExchangePlan verifyPlan = plan;
    verifyPlan.sendCellIds = plan.receiveCellIds;
    verifyPlan.receiveCellIds = plan.receiveCellIds;
    const ONEFLOW::BoundaryBufferView verified = backend.PackBoundary(
        cacheKey, verifyPlan, ONEFLOW::BoundaryTransferMode::HostStaged );

    double maxError = 0.0;
    for ( int equation = 0; equation < nEq; ++ equation )
    {
        for ( int interfaceId = 0; interfaceId < 2; ++ interfaceId )
        {
            const ONEFLOW::Real expected =
                primitive[
                    equation * totalCells
                    + plan.sendCellIds[ interfaceId ] ];
            const ONEFLOW::Real actual =
                verified.values[ equation * 2 + interfaceId ];
            maxError = std::max(
                maxError,
                static_cast< double >( std::abs( actual - expected ) ) );
        }
    }

    ONEFLOW::BoundaryFieldExchangePlan directPlan = plan;
    directPlan.sendCellIds = { 0, 1 };
    directPlan.receiveCellIds = { 3, 2 };
    const ONEFLOW::BoundaryBufferView directSend = backend.PackBoundary(
        cacheKey, directPlan, ONEFLOW::BoundaryTransferMode::DeviceDirect );
    ONEFLOW::BoundaryBufferView directReceive =
        backend.PrepareBoundaryReceive(
            cacheKey, directPlan,
            ONEFLOW::BoundaryTransferMode::DeviceDirect );
    ONEFLOW::AccelRuntime::Instance().Backend().Copy(
        directReceive.values,
        directSend.values,
        directSend.valueCount * sizeof( ONEFLOW::Real ),
        ONEFLOW::AccelCopyKind::DeviceToDevice );
    backend.UnpackBoundary(
        cacheKey, directPlan,
        ONEFLOW::BoundaryTransferMode::DeviceDirect );

    verifyPlan.sendCellIds = directPlan.receiveCellIds;
    verifyPlan.receiveCellIds = directPlan.receiveCellIds;
    const ONEFLOW::BoundaryBufferView directVerified =
        backend.PackBoundary(
            cacheKey, verifyPlan,
            ONEFLOW::BoundaryTransferMode::HostStaged );
    double directError = 0.0;
    for ( int equation = 0; equation < nEq; ++ equation )
    {
        for ( int interfaceId = 0; interfaceId < 2; ++ interfaceId )
        {
            const ONEFLOW::Real expected =
                primitive[
                    equation * totalCells
                    + directPlan.sendCellIds[ interfaceId ] ];
            const ONEFLOW::Real actual =
                directVerified.values[ equation * 2 + interfaceId ];
            directError = std::max(
                directError,
                static_cast< double >( std::abs( actual - expected ) ) );
        }
    }
    backend.UnbindState( cacheKey );
    if ( maxError > 1.0e-15 || directError > 1.0e-15 )
    {
        std::fprintf(
            stderr,
            "Main solver HIP boundary exchange FAIL: %.3e %.3e\n",
            maxError, directError );
        return false;
    }
    std::printf(
        "OneFLOW HIP main solver boundary exchange: PASS "
        "(host-staged %.3e, device-direct %.3e)\n",
        maxError, directError );
    return true;
}

} // namespace

int main()
{
    ONEFLOW::InitializeAccelRuntime( 0, 1 );
    try
    {
        ONEFLOW::RunHipBackendSelfTest();

        bool ok = true;
        ok = TestScalarConvection() && ok;
        ok = TestEulerRusanov() && ok;
        ok = TestMainSolverFiveEquationLaxFriedrichs() && ok;
        ok = TestMainSolverGradient() && ok;
        ok = TestMainSolverBoundarySemantics() && ok;
        ok = TestMainSolverAleAreaContract() && ok;
        ok = TestMainSolverStateUpdate() && ok;
        ok = TestMainSolverBoundaryExchange() && ok;

        ONEFLOW::FinalizeAccelRuntime();
        return ok ? 0 : 1;
    }
    catch ( const std::exception & error )
    {
        std::fprintf( stderr, "OneFLOW HIP smoke: FAIL: %s\n", error.what() );
        ONEFLOW::FinalizeAccelRuntime();
        return 1;
    }
}
