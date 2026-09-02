#include "OneDEulerBackend.h"
#include "OneDEulerProfile.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

constexpr double Gamma = 1.4;
constexpr double DtScale = 0.05;

int ParsePositive( const char * text, const char * name )
{
    char * end = nullptr;
    const long value = std::strtol( text, &end, 10 );
    if ( end == text || *end != '\0' || value <= 0 || value > 100000000 )
    {
        throw std::invalid_argument( std::string( "invalid " ) + name );
    }
    return static_cast< int >( value );
}

std::vector< double > InitialState( int nx )
{
    std::vector< double > state( 3 * nx );
    for ( int cell = 0; cell < nx; ++ cell )
    {
        const double x = ( cell + 0.5 ) / nx;
        const double rho = 1.0 + 0.1 * std::sin( 2.0 * M_PI * x );
        const double velocity = 0.2 + 0.05 * std::cos( 2.0 * M_PI * x );
        const double pressure = 1.0 + 0.08 * std::sin( 4.0 * M_PI * x );
        state[ cell ] = rho;
        state[ nx + cell ] = rho * velocity;
        state[ 2 * nx + cell ] =
            pressure / ( Gamma - 1.0 ) + 0.5 * rho * velocity * velocity;
    }
    return state;
}

double MaxAbsError(
    const std::vector< double > & reference,
    const std::vector< double > & candidate )
{
    double error = 0.0;
    for ( std::size_t i = 0; i < reference.size(); ++ i )
    {
        error = std::max( error, std::abs( reference[ i ] - candidate[ i ] ) );
    }
    return error;
}

double Checksum( const std::vector< double > & values )
{
    double checksum = 0.0;
    for ( std::size_t i = 0; i < values.size(); ++ i )
    {
        checksum += values[ i ] * ( 1.0 + static_cast< double >( i % 17 ) );
    }
    return checksum;
}

void AdvanceCpu(
    oneflow_1d::CpuEulerBackend & backend, std::vector< double > & state,
    int nx, int steps, oneflow_1d::EulerBoundary boundary )
{
    oneflow_1d::EulerTrace trace;
    const double dt = DtScale / nx;
    const double dx = 1.0 / nx;
    for ( int step = 0; step < steps; ++ step )
    {
        backend.Step( state.data(), nx, Gamma, dt, dx, boundary, trace );
        std::copy(
            trace.state.begin() + oneflow_1d::EulerRkStages * oneflow_1d::EulerComponents * nx,
            trace.state.end(),
            state.begin() );
    }
    return;
}

struct HipResult
{
    std::vector< double > state;
    oneflow_1d::EulerProfileTiming timing;
};

HipResult AdvanceHip(
    std::vector< double > state, int nx, int steps,
    oneflow_1d::EulerBoundary boundary )
{
    oneflow_1d::EulerTrace trace;
    oneflow_1d::EulerProfileTiming total;
    const double dt = DtScale / nx;
    const double dx = 1.0 / nx;
    for ( int step = 0; step < steps; ++ step )
    {
        oneflow_1d::EulerProfileTiming timing;
        oneflow_1d::ProfiledHipEulerStep(
            state.data(), nx, Gamma, dt, dx, boundary, trace, timing );
        total.hostMilliseconds += timing.hostMilliseconds;
        total.allocationMilliseconds += timing.allocationMilliseconds;
        total.h2dMilliseconds += timing.h2dMilliseconds;
        total.kernelMilliseconds += timing.kernelMilliseconds;
        total.d2hMilliseconds += timing.d2hMilliseconds;
        total.allocationCount += timing.allocationCount;
        total.kernelLaunches += timing.kernelLaunches;
        total.deviceTimingAvailable =
            total.deviceTimingAvailable || timing.deviceTimingAvailable;
        std::copy(
            trace.state.begin() + oneflow_1d::EulerRkStages * oneflow_1d::EulerComponents * nx,
            trace.state.end(),
            state.begin() );
        if ( step == 0 )
        {
            std::copy(
                timing.deviceName,
                timing.deviceName + sizeof( timing.deviceName ),
                total.deviceName );
            std::copy(
                timing.architecture,
                timing.architecture + sizeof( timing.architecture ),
                total.architecture );
        }
    }
    return { std::move( state ), total };
}

} // namespace

int main( int argc, char ** argv )
{
    try
    {
        const int nx = argc > 1 ? ParsePositive( argv[ 1 ], "nx" ) : 65536;
        const int steps = argc > 2 ? ParsePositive( argv[ 2 ], "steps" ) : 20;
        const int repeats = argc > 3 ? ParsePositive( argv[ 3 ], "repeats" ) : 3;
        const int warmup = argc > 4 ? ParsePositive( argv[ 4 ], "warmup" ) : 2;
        const auto boundary = oneflow_1d::EulerBoundary::Periodic;
        oneflow_1d::CpuEulerBackend cpuBackend;
        const std::vector< double > initial = InitialState( nx );

        for ( int i = 0; i < warmup; ++ i )
        {
            std::vector< double > cpuState = initial;
            AdvanceCpu( cpuBackend, cpuState, nx, steps, boundary );
            AdvanceHip( initial, nx, steps, boundary );
        }

        double cpuMilliseconds = 0.0;
        double hipMilliseconds = 0.0;
        oneflow_1d::EulerProfileTiming hipTiming;
        std::vector< double > cpuFinal;
        std::vector< double > hipFinal;
        for ( int repeat = 0; repeat < repeats; ++ repeat )
        {
            std::vector< double > cpuState = initial;
            const auto cpuBegin = std::chrono::steady_clock::now();
            AdvanceCpu( cpuBackend, cpuState, nx, steps, boundary );
            cpuMilliseconds += std::chrono::duration< double, std::milli >(
                std::chrono::steady_clock::now() - cpuBegin ).count();

            const auto hip = AdvanceHip( initial, nx, steps, boundary );
            hipMilliseconds += hip.timing.hostMilliseconds;
            if ( repeat == 0 ) hipTiming = hip.timing;
            else
            {
                hipTiming.hostMilliseconds += hip.timing.hostMilliseconds;
                hipTiming.allocationMilliseconds += hip.timing.allocationMilliseconds;
                hipTiming.h2dMilliseconds += hip.timing.h2dMilliseconds;
                hipTiming.kernelMilliseconds += hip.timing.kernelMilliseconds;
                hipTiming.d2hMilliseconds += hip.timing.d2hMilliseconds;
                hipTiming.allocationCount += hip.timing.allocationCount;
                hipTiming.kernelLaunches += hip.timing.kernelLaunches;
                hipTiming.deviceTimingAvailable =
                    hipTiming.deviceTimingAvailable || hip.timing.deviceTimingAvailable;
            }
            cpuFinal = std::move( cpuState );
            hipFinal = hip.state;
        }

        const double hipKernelRatio =
            hipTiming.hostMilliseconds > 0.0
            ? hipTiming.kernelMilliseconds / hipTiming.hostMilliseconds
            : 0.0;
        const double hipCopyRatio =
            hipTiming.hostMilliseconds > 0.0
            ? ( hipTiming.h2dMilliseconds + hipTiming.d2hMilliseconds )
                / hipTiming.hostMilliseconds
            : 0.0;
        const double hipAllocationRatio =
            hipTiming.hostMilliseconds > 0.0
            ? hipTiming.allocationMilliseconds / hipTiming.hostMilliseconds
            : 0.0;
        std::cout << std::setprecision( 10 ) << std::fixed;
        std::cout << "OneFLOW 1D Euler backend benchmark\n";
        std::cout << "nx=" << nx << " steps=" << steps
                  << " repeats=" << repeats << " warmup=" << warmup << "\n";
        std::cout << "device=" << hipTiming.deviceName
                  << " arch=" << hipTiming.architecture << "\n";
        std::cout << "cpu_total_ms=" << cpuMilliseconds
                  << " cpu_per_step_ms=" << cpuMilliseconds / repeats / steps
                  << "\n";
        std::cout << "hip_total_ms=" << hipMilliseconds
                  << " hip_per_step_ms=" << hipMilliseconds / repeats / steps
                  << "\n";
        std::cout << "hip_kernel_ms=" << hipTiming.kernelMilliseconds
                  << " hip_h2d_ms=" << hipTiming.h2dMilliseconds
                  << " hip_d2h_ms=" << hipTiming.d2hMilliseconds
                  << " hip_alloc_ms=" << hipTiming.allocationMilliseconds << "\n";
        std::cout << "hip_kernel_ratio=" << hipKernelRatio
                  << " hip_copy_ratio=" << hipCopyRatio
                  << " hip_alloc_ratio=" << hipAllocationRatio << "\n";
        std::cout << "hip_kernel_launches=" << hipTiming.kernelLaunches
                  << " hip_allocations=" << hipTiming.allocationCount << "\n";
        std::cout << "speedup_cpu_over_hip="
                  << cpuMilliseconds / hipMilliseconds << "\n";
        std::cout << "final_max_abs_error="
                  << MaxAbsError( cpuFinal, hipFinal ) << "\n";
        std::cout << "cpu_checksum=" << Checksum( cpuFinal )
                  << " hip_checksum=" << Checksum( hipFinal ) << "\n";
        return 0;
    }
    catch ( const std::exception & error )
    {
        std::cerr << "Euler benchmark failed: " << error.what() << "\n";
        return 1;
    }
}
