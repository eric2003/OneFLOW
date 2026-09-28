#include "OneDEulerBackend.h"
#include "CpuFluxBackend.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace
{

using oneflow_1d::CpuEulerBackend;
using oneflow_1d::EulerBoundary;
using oneflow_1d::EulerComponents;
using oneflow_1d::EulerTrace;

constexpr double kGamma = 1.4;
constexpr double kPi = 3.14159265358979323846;

// -----------------------------------------------------------------------
// Port's face indexing (OneDEulerPersistent.hip Face kernel):
//   face i has left = cell (i-1+nx)%nx, right = cell i%nx
//   There are nx+1 faces (0..nx), face 0 ≡ face nx for periodic.
//
// Bridge test's face f:
//   left = cell f, right = cell (f+1)%nx
//   There are nx faces (0..nx-1).
//
// Mapping: bridge face f → port face (f+1) for f=0..nx-2
//          bridge face nx-1 → port face 0
// -----------------------------------------------------------------------

int MapBridgeFaceToPortFace( int f, int nx )
{
    return ( f + 1 ) % nx;  // f=nx-1 → 0, f=0 → 1, etc.
}

void BuildFaceViews(
    const double * cellState, int nx,
    std::vector<ONEFLOW::Real> & qLeft,
    std::vector<ONEFLOW::Real> & qRight )
{
    const int nFaces = nx;
    qLeft.resize( EulerComponents * nFaces );
    qRight.resize( EulerComponents * nFaces );

    for ( int face = 0; face < nFaces; ++ face )
    {
        const int il = face;
        const int ir = ( face + 1 ) % nx;
        for ( int c = 0; c < EulerComponents; ++ c )
        {
            qLeft[ c * nFaces + face ] = cellState[ c * nx + il ];
            qRight[ c * nFaces + face ] = cellState[ c * nx + ir ];
        }
    }
}

void ExtractTraceFlux(
    const EulerTrace & trace, int nx,
    std::vector<ONEFLOW::Real> & traceFlux )
{
    const int nFaces = nx;
    const int nFacesPort = nx + 1;  // port stores nx+1 faces
    traceFlux.resize( EulerComponents * nFaces );

    for ( int face = 0; face < nFaces; ++ face )
    {
        const int portFace = MapBridgeFaceToPortFace( face, nx );
        for ( int eq = 0; eq < EulerComponents; ++ eq )
        {
            // Stage 0: numericalFlux[eq * nFacesPort + portFace]
            traceFlux[ eq * nFaces + face ] =
                trace.numericalFlux[ eq * nFacesPort + portFace ];
        }
    }
}

bool CompareFlux( int nx )
{
    const double dx = 1.0 / nx;
    const double dt = 0.0005;

    // Smooth initial condition
    std::vector<double> initialState( EulerComponents * nx );
    for ( int cell = 0; cell < nx; ++ cell )
    {
        const double x = ( cell + 0.5 ) * dx;
        const double density = 1.0 + 0.08 * std::sin( 2.0 * kPi * x );
        const double velocity = 0.2 + 0.04 * std::cos( 2.0 * kPi * x );
        const double pressure = 1.0 + 0.05 * std::sin( 4.0 * kPi * x );
        initialState[ cell ] = density;
        initialState[ nx + cell ] = density * velocity;
        initialState[ 2 * nx + cell ] = pressure / ( kGamma - 1.0 )
            + 0.5 * density * velocity * velocity;
    }

    // Path A: EulerBackend::Step (port)
    CpuEulerBackend portBackend;
    EulerTrace trace;
    portBackend.Step(
        initialState.data(), nx, kGamma, dt, dx,
        EulerBoundary::Periodic, trace );

    std::vector<ONEFLOW::Real> traceFlux;
    ExtractTraceFlux( trace, nx, traceFlux );

    // Path B: CpuFluxBackend::CalcInvFlux (accel)
    std::vector<ONEFLOW::Real> qLeft, qRight;
    BuildFaceViews( initialState.data(), nx, qLeft, qRight );

    std::vector<ONEFLOW::Real> xNormal( nx, 1.0 );
    std::vector<ONEFLOW::Real> area( nx, 1.0 );

    ONEFLOW::FaceStateView state;
    state.nFaces = nx;
    state.nEquations = EulerComponents;
    state.qLeft = qLeft.data();
    state.qRight = qRight.data();
    state.xNormal = xNormal.data();
    state.faceArea = area.data();
    state.gamma = kGamma;

    std::vector<ONEFLOW::Real> accelFlux( EulerComponents * nx );
    ONEFLOW::FaceFluxView fluxView{ nx, EulerComponents, accelFlux.data() };

    ONEFLOW::CpuFluxBackend accelBackend;
    accelBackend.CalcInvFlux( state, fluxView, 0 );

    // Compare
    double maxAbsDiff = 0.0;
    double maxRelDiff = 0.0;
    int diffCount = 0;
    int firstBadEq = -1, firstBadFace = -1;

    for ( int face = 0; face < nx; ++ face )
    {
        for ( int eq = 0; eq < EulerComponents; ++ eq )
        {
            const int idx = eq * nx + face;
            const double ref = traceFlux[ idx ];
            const double val = accelFlux[ idx ];
            const double absDiff = std::abs( val - ref );
            const double relDiff = absDiff
                / std::max( 1.0, std::abs( ref ) );

            if ( absDiff > maxAbsDiff ) maxAbsDiff = absDiff;
            if ( relDiff > maxRelDiff ) maxRelDiff = relDiff;
            if ( absDiff > 1.0e-15 )
            {
                ++ diffCount;
                if ( firstBadEq < 0 ) { firstBadEq = eq; firstBadFace = face; }
            }
        }
    }

    std::printf( "  nx=%d: max abs diff=%.3e, max rel diff=%.3e, diffs>1e-15=%d",
        nx, maxAbsDiff, maxRelDiff, diffCount );
    if ( diffCount > 0 )
        std::printf( " (first: eq=%d face=%d ref=%.6f val=%.6f)",
            firstBadEq, firstBadFace,
            traceFlux[ firstBadEq * nx + firstBadFace ],
            accelFlux[ firstBadEq * nx + firstBadFace ] );

    const bool passed = ( maxAbsDiff <= 1.0e-14 && diffCount == 0 );
    std::printf( "  [%s]\n", passed ? "PASS" : "FAIL" );
    return passed;
}

} // namespace

int main()
{
    std::printf( "FluxBackend <-> EulerBackend bridge test\n" );
    std::printf( "========================================\n" );

    bool allPassed = true;
    allPassed = CompareFlux( 32 ) && allPassed;
    allPassed = CompareFlux( 64 ) && allPassed;
    allPassed = CompareFlux( 128 ) && allPassed;
    allPassed = CompareFlux( 256 ) && allPassed;

    std::printf( "\nBridge test: %s\n", allPassed ? "ALL PASSED" : "FAILED" );
    return allPassed ? 0 : 1;
}
