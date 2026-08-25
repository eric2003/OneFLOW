/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\---------------------------------------------------------------------------*/

#include "AccelRuntime.h"
#include "CpuFluxBackend.h"
#include "HipFluxBackend.h"
#include "HipKernel.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>

int main()
{
    ONEFLOW::InitializeAccelRuntime( 0, 1 );
    try
    {
        ONEFLOW::RunHipBackendSelfTest();

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

        std::vector< ONEFLOW::Real > cpuFluxValues( nFaces );
        std::vector< ONEFLOW::Real > hipFluxValues( nFaces );
        ONEFLOW::FaceFluxView cpuFlux{ nFaces, 1, cpuFluxValues.data() };
        ONEFLOW::FaceFluxView hipFlux{ nFaces, 1, hipFluxValues.data() };

        ONEFLOW::CpuFluxBackend cpuFluxBackend;
        ONEFLOW::HipFluxBackend hipFluxBackend;
        cpuFluxBackend.CalcInvFlux( state, cpuFlux, 0 );
        hipFluxBackend.CalcInvFlux( state, hipFlux, 0 );

        std::vector< ONEFLOW::Real > cpuResidualValues( nCells, 0.0 );
        std::vector< ONEFLOW::Real > hipResidualValues( nCells, 0.0 );
        ONEFLOW::FaceConnectivityView connectivity;
        connectivity.nFaces = nFaces;
        connectivity.nBoundaryFaces = nBoundaryFaces;
        connectivity.leftCell = leftCell.data();
        connectivity.rightCell = rightCell.data();
        ONEFLOW::ResidualView cpuResidual{ nCells, 1, cpuResidualValues.data() };
        ONEFLOW::ResidualView hipResidual{ nCells, 1, hipResidualValues.data() };
        cpuFluxBackend.AddFaceFlux( cpuFlux, connectivity, cpuResidual );
        hipFluxBackend.AddFaceFlux( hipFlux, connectivity, hipResidual );

        double maxFluxError = 0.0;
        double maxResidualError = 0.0;
        for ( int face = 0; face < nFaces; ++ face )
        {
            maxFluxError = std::max(
                maxFluxError,
                std::abs( cpuFluxValues[ face ] - hipFluxValues[ face ] ) );
        }
        for ( int cell = 0; cell < nCells; ++ cell )
        {
            maxResidualError = std::max(
                maxResidualError,
                std::abs( cpuResidualValues[ cell ] - hipResidualValues[ cell ] ) );
        }
        if ( maxFluxError > 1.0e-15 || maxResidualError > 1.0e-15 )
        {
            std::fprintf(
                stderr,
                "OneFLOW HIP flux smoke: FAIL (flux %.3e, residual %.3e)\n",
                maxFluxError,
                maxResidualError );
            ONEFLOW::FinalizeAccelRuntime();
            return 1;
        }
        std::printf(
            "OneFLOW HIP flux smoke: PASS (flux %.3e, residual %.3e)\n",
            maxFluxError,
            maxResidualError );

        ONEFLOW::FinalizeAccelRuntime();
        return 0;
    }
    catch ( const std::exception & error )
    {
        std::fprintf( stderr, "OneFLOW HIP smoke: FAIL: %s\n", error.what() );
        ONEFLOW::FinalizeAccelRuntime();
        return 1;
    }
}
