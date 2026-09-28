/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\*---------------------------------------------------------------------------*/

#include "AccelRuntime.h"
#include "CpuFluxBackend.h"
#include "HipFluxBackend.h"
#include "HipKernel.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
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
