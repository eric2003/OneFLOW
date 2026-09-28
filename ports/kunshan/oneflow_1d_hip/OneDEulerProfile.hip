#include "OneDEulerProfile.h"

#include <hip/hip_runtime.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <string>
#include <utility>

namespace oneflow_1d
{
namespace
{

using Clock = std::chrono::steady_clock;

double Milliseconds( Clock::time_point begin, Clock::time_point end )
{
    return std::chrono::duration< double, std::milli >( end - begin ).count();
}

void Check( hipError_t error, const char * operation )
{
    if ( error != hipSuccess )
    {
        throw std::runtime_error(
            std::string( operation ) + ": " + hipGetErrorString( error ) );
    }
}

__device__ int CellIndex( int component, int cell, int nx )
{
    return component * nx + cell;
}

__device__ int FaceIndex( int component, int face, int nx )
{
    return component * ( nx + 1 ) + face;
}

__device__ void PhysicalFlux(
    const double * state, int cell, int nx, double gamma, double flux[ 3 ],
    double & waveSpeed )
{
    const double rho = state[ CellIndex( 0, cell, nx ) ];
    const double momentum = state[ CellIndex( 1, cell, nx ) ];
    const double energy = state[ CellIndex( 2, cell, nx ) ];
    const double velocity = momentum / rho;
    const double pressure =
        ( gamma - 1.0 ) * ( energy - 0.5 * momentum * momentum / rho );
    const double soundSpeed = sqrt( gamma * pressure / rho );
    flux[ 0 ] = momentum;
    flux[ 1 ] = momentum * velocity + pressure;
    flux[ 2 ] = ( energy + pressure ) * velocity;
    waveSpeed = fabs( velocity ) + soundSpeed;
}

__global__ void FaceKernel(
    const double * state, int nx, double gamma, int boundary, double * left,
    double * right, double * flux )
{
    const int face = blockIdx.x * blockDim.x + threadIdx.x;
    if ( face > nx ) return;

    int leftCell = face - 1;
    int rightCell = face;
    if ( boundary == static_cast< int >( EulerBoundary::Periodic ) )
    {
        leftCell = ( face + nx - 1 ) % nx;
        rightCell = face % nx;
    }
    else
    {
        leftCell = leftCell < 0 ? 0 : leftCell;
        rightCell = rightCell >= nx ? nx - 1 : rightCell;
    }

    double leftFlux[ 3 ];
    double rightFlux[ 3 ];
    double leftSpeed = 0.0;
    double rightSpeed = 0.0;
    PhysicalFlux( state, leftCell, nx, gamma, leftFlux, leftSpeed );
    PhysicalFlux( state, rightCell, nx, gamma, rightFlux, rightSpeed );
    const double speed = leftSpeed > rightSpeed ? leftSpeed : rightSpeed;
    for ( int component = 0; component < EulerComponents; ++ component )
    {
        const double leftState = state[ CellIndex( component, leftCell, nx ) ];
        const double rightState = state[ CellIndex( component, rightCell, nx ) ];
        const int index = FaceIndex( component, face, nx );
        left[ index ] = leftState;
        right[ index ] = rightState;
        flux[ index ] = 0.5 * ( leftFlux[ component ] + rightFlux[ component] )
            - 0.5 * speed * ( rightState - leftState );
    }
}

__global__ void ResidualKernel(
    const double * flux, int nx, double inverseDx, double * residual )
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    const int count = EulerComponents * nx;
    if ( index >= count ) return;
    const int component = index / nx;
    const int cell = index - component * nx;
    residual[ index ] = -(
        flux[ FaceIndex( component, cell + 1, nx ) ]
        - flux[ FaceIndex( component, cell, nx ) ] ) * inverseDx;
}

__global__ void UpdateKernel(
    const double * base, const double * current, const double * residual,
    int count, double dt, double baseWeight, double updateWeight, double * next )
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    if ( index < count )
    {
        next[ index ] = baseWeight * base[ index ]
            + updateWeight * ( current[ index ] + dt * residual[ index ] );
    }
}

void Free( double * pointer ) noexcept
{
    if ( pointer != nullptr ) static_cast< void >( hipFree( pointer ) );
}

} // namespace

void ProfiledHipEulerStep(
    const double * input, int nx, double gamma, double dt, double dx,
    EulerBoundary boundary, EulerTrace & trace, EulerProfileTiming & timing )
{
    timing = EulerProfileTiming {};
    const Clock::time_point hostBegin = Clock::now();
    if ( input == nullptr || nx <= 0 || gamma <= 1.0 || dt <= 0.0 || dx <= 0.0 )
    {
        throw std::invalid_argument( "invalid profiled HIP Euler input" );
    }

    int device = 0;
    Check( hipGetDevice( &device ), "hipGetDevice" );
    hipDeviceProp_t properties {};
    Check( hipGetDeviceProperties( &properties, device ), "hipGetDeviceProperties" );
    std::strncpy(
        timing.deviceName, properties.name, sizeof( timing.deviceName ) - 1 );
    std::strncpy(
        timing.architecture,
        properties.gcnArchName,
        sizeof( timing.architecture ) - 1 );

    ResizeEulerTrace( trace, nx );
    const int faceValues = EulerComponents * ( nx + 1 );
    const int cellValues = EulerComponents * nx;
    const int threads = 256;
    const int faceBlocks = ( faceValues + threads - 1 ) / threads;
    const int cellBlocks = ( cellValues + threads - 1 ) / threads;
    const std::size_t faceBytes = faceValues * sizeof( double );
    const std::size_t cellBytes = cellValues * sizeof( double );

    double * base = nullptr;
    double * current = nullptr;
    double * next = nullptr;
    double * left = nullptr;
    double * right = nullptr;
    double * flux = nullptr;
    double * residual = nullptr;
    hipEvent_t kernelBegin = nullptr;
    hipEvent_t kernelEnd = nullptr;

    try
    {
        const Clock::time_point allocationBegin = Clock::now();
        Check( hipMalloc( &base, cellBytes ), "hipMalloc base" );
        Check( hipMalloc( &current, cellBytes ), "hipMalloc current" );
        Check( hipMalloc( &next, cellBytes ), "hipMalloc next" );
        Check( hipMalloc( &left, faceBytes ), "hipMalloc left" );
        Check( hipMalloc( &right, faceBytes ), "hipMalloc right" );
        Check( hipMalloc( &flux, faceBytes ), "hipMalloc flux" );
        Check( hipMalloc( &residual, cellBytes ), "hipMalloc residual" );
        timing.allocationMilliseconds =
            Milliseconds( allocationBegin, Clock::now() );
        timing.allocationCount = 7;

        const Clock::time_point h2dBegin = Clock::now();
        Check(
            hipMemcpy( base, input, cellBytes, hipMemcpyHostToDevice ),
            "hipMemcpy base" );
        Check(
            hipMemcpy( current, input, cellBytes, hipMemcpyHostToDevice ),
            "hipMemcpy current" );
        timing.h2dMilliseconds = Milliseconds( h2dBegin, Clock::now() );

        Check( hipEventCreate( &kernelBegin ), "hipEventCreate begin" );
        Check( hipEventCreate( &kernelEnd ), "hipEventCreate end" );
        timing.deviceTimingAvailable = true;
        std::copy( input, input + cellValues, trace.state.begin() );

        for ( int stage = 0; stage < EulerRkStages; ++ stage )
        {
            Check( hipEventRecord( kernelBegin ), "hipEventRecord begin" );
            hipLaunchKernelGGL(
                FaceKernel,
                dim3( faceBlocks ),
                dim3( threads ),
                0,
                0,
                current,
                nx,
                gamma,
                static_cast< int >( boundary ),
                left,
                right,
                flux );
            Check( hipGetLastError(), "FaceKernel launch" );
            hipLaunchKernelGGL(
                ResidualKernel,
                dim3( cellBlocks ),
                dim3( threads ),
                0,
                0,
                flux,
                nx,
                1.0 / dx,
                residual );
            Check( hipGetLastError(), "ResidualKernel launch" );
            const double baseWeight = stage == 1 ? 0.75
                : ( stage == 2 ? 1.0 / 3.0 : 0.0 );
            const double updateWeight = stage == 0 ? 1.0
                : ( stage == 1 ? 0.25 : 2.0 / 3.0 );
            hipLaunchKernelGGL(
                UpdateKernel,
                dim3( cellBlocks ),
                dim3( threads ),
                0,
                0,
                base,
                current,
                residual,
                cellValues,
                dt,
                baseWeight,
                updateWeight,
                next );
            Check( hipGetLastError(), "UpdateKernel launch" );
            Check( hipEventRecord( kernelEnd ), "hipEventRecord end" );
            Check( hipEventSynchronize( kernelEnd ), "hipEventSynchronize" );
            float elapsed = 0.0F;
            Check(
                hipEventElapsedTime( &elapsed, kernelBegin, kernelEnd ),
                "hipEventElapsedTime" );
            timing.kernelMilliseconds += elapsed;
            timing.kernelLaunches += 3;

            const Clock::time_point d2hBegin = Clock::now();
            const int faceOffset = stage * faceValues;
            const int cellOffset = stage * cellValues;
            Check(
                hipMemcpy(
                    trace.faceLeft.data() + faceOffset,
                    left,
                    faceBytes,
                    hipMemcpyDeviceToHost ),
                "hipMemcpy face left" );
            Check(
                hipMemcpy(
                    trace.faceRight.data() + faceOffset,
                    right,
                    faceBytes,
                    hipMemcpyDeviceToHost ),
                "hipMemcpy face right" );
            Check(
                hipMemcpy(
                    trace.numericalFlux.data() + faceOffset,
                    flux,
                    faceBytes,
                    hipMemcpyDeviceToHost ),
                "hipMemcpy flux" );
            Check(
                hipMemcpy(
                    trace.residual.data() + cellOffset,
                    residual,
                    cellBytes,
                    hipMemcpyDeviceToHost ),
                "hipMemcpy residual" );
            Check(
                hipMemcpy(
                    trace.state.data() + ( stage + 1 ) * cellValues,
                    next,
                    cellBytes,
                    hipMemcpyDeviceToHost ),
                "hipMemcpy state" );
            timing.d2hMilliseconds += Milliseconds( d2hBegin, Clock::now() );
            std::swap( current, next );
        }
    }
    catch ( ... )
    {
        if ( kernelEnd != nullptr ) static_cast< void >( hipEventDestroy( kernelEnd ) );
        if ( kernelBegin != nullptr )
            static_cast< void >( hipEventDestroy( kernelBegin ) );
        Free( residual );
        Free( flux );
        Free( right );
        Free( left );
        Free( next );
        Free( current );
        Free( base );
        throw;
    }

    Check( hipEventDestroy( kernelEnd ), "hipEventDestroy end" );
    Check( hipEventDestroy( kernelBegin ), "hipEventDestroy begin" );
    Free( residual );
    Free( flux );
    Free( right );
    Free( left );
    Free( next );
    Free( current );
    Free( base );
    timing.hostMilliseconds = Milliseconds( hostBegin, Clock::now() );
}

} // namespace oneflow_1d
