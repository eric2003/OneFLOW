/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\*---------------------------------------------------------------------------*/

#include "HipKernel.h"

#include <hip/hip_runtime.h>

#include <cmath>
#include <algorithm>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

namespace ONEFLOW
{
namespace
{

__global__ void AddScaledKernel(
    const double * a,
    const double * b,
    double * c,
    int count )
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if ( i < count ) c[ i ] = a[ i ] + 2.0 * b[ i ];
}

void CheckHip( hipError_t status, const char * operation )
{
    if ( status == hipSuccess ) return;
    throw std::runtime_error(
        std::string( operation ) + ": " + hipGetErrorString( status ) );
}

}

void RunHipBackendSelfTest()
{
    constexpr int count = 4096;
    std::vector< double > a( count );
    std::vector< double > b( count );
    std::vector< double > expected( count );
    std::vector< double > result( count, 0.0 );

    for ( int i = 0; i < count; ++ i )
    {
        a[ i ] = 0.25 * static_cast< double >( i ) - 3.0;
        b[ i ] = std::sin( 0.01 * static_cast< double >( i ) );
        expected[ i ] = a[ i ] + 2.0 * b[ i ];
    }

    double * deviceA = nullptr;
    double * deviceB = nullptr;
    double * deviceC = nullptr;

    CheckHip( hipMalloc( &deviceA, count * sizeof( double ) ), "hipMalloc(a)" );
    CheckHip( hipMalloc( &deviceB, count * sizeof( double ) ), "hipMalloc(b)" );
    CheckHip( hipMalloc( &deviceC, count * sizeof( double ) ), "hipMalloc(c)" );

    try
    {
        CheckHip(
            hipMemcpy(
                deviceA, a.data(), count * sizeof( double ),
                hipMemcpyHostToDevice ),
            "hipMemcpy(a)" );
        CheckHip(
            hipMemcpy(
                deviceB, b.data(), count * sizeof( double ),
                hipMemcpyHostToDevice ),
            "hipMemcpy(b)" );

        const int blockSize = 256;
        const int gridSize = ( count + blockSize - 1 ) / blockSize;
        hipLaunchKernelGGL(
            AddScaledKernel,
            dim3( gridSize ),
            dim3( blockSize ),
            0,
            0,
            deviceA,
            deviceB,
            deviceC,
            count );
        CheckHip( hipGetLastError(), "AddScaledKernel launch" );
        CheckHip( hipDeviceSynchronize(), "AddScaledKernel synchronize" );
        CheckHip(
            hipMemcpy(
                result.data(), deviceC, count * sizeof( double ),
                hipMemcpyDeviceToHost ),
            "hipMemcpy(c)" );
    }
    catch ( ... )
    {
        hipFree( deviceC );
        hipFree( deviceB );
        hipFree( deviceA );
        throw;
    }

    CheckHip( hipFree( deviceC ), "hipFree(c)" );
    CheckHip( hipFree( deviceB ), "hipFree(b)" );
    CheckHip( hipFree( deviceA ), "hipFree(a)" );

    double maxError = 0.0;
    for ( int i = 0; i < count; ++ i )
    {
        maxError = std::max( maxError, std::abs( result[ i ] - expected[ i ] ) );
    }
    if ( maxError > 1.0e-14 )
    {
        throw std::runtime_error(
            "HIP self-test failed: max absolute error = "
            + std::to_string( maxError ) );
    }

    std::printf(
        "OneFLOW HIP self-test: PASS (double, %d values, max error %.3e)\n",
        count,
        maxError );
}

}
