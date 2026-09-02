#include "OneDEulerBackend.h"

#include <hip/hip_runtime.h>

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

double Checksum( const std::vector< double > & values )
{
    double checksum = 0.0;
    for ( std::size_t i = 0; i < values.size(); ++ i )
    {
        checksum += values[ i ] * ( 1.0 + static_cast< double >( i % 17 ) );
    }
    return checksum;
}

double MaxAbs( const std::vector< double > & values )
{
    double maximum = 0.0;
    for ( const double value : values ) maximum = std::max( maximum, std::abs( value ) );
    return maximum;
}

} // namespace

int main( int argc, char ** argv )
{
    try
    {
        const int nx = argc > 1 ? ParsePositive( argv[ 1 ], "nx" ) : 1048576;
        const int steps = argc > 2 ? ParsePositive( argv[ 2 ], "steps" ) : 100;
        const int repeats = argc > 3 ? ParsePositive( argv[ 3 ], "repeats" ) : 3;
        const int warmup = argc > 4 ? ParsePositive( argv[ 4 ], "warmup" ) : 1;
        const auto initial = InitialState( nx );
        oneflow_1d::HipEulerBackend backend;

        hipDeviceProp_t properties {};
        const hipError_t deviceResult = hipGetDeviceProperties( &properties, 0 );
        if ( deviceResult != hipSuccess )
        {
            throw std::runtime_error( hipGetErrorString( deviceResult ) );
        }

        double totalAdvanceMs = 0.0;
        double totalKernelMs = 0.0;
        int totalLaunches = 0;
        int totalSyncs = 0;
        std::vector< double > finalState;
        for ( int repeat = 0; repeat < warmup + repeats; ++ repeat )
        {
            const auto begin = std::chrono::steady_clock::now();
            const oneflow_1d::EulerProblem problem{
                nx, Gamma, DtScale / nx, 1.0 / nx,
                oneflow_1d::EulerBoundary::Periodic };
            auto state = backend.CreateState( problem );
            backend.Upload( *state, initial.data() );
            oneflow_1d::EulerRunStats stats;
            backend.Advance(
                *state, steps,
                { oneflow_1d::EulerRunMode::NoTrace, nullptr, &stats } );
            finalState.resize( initial.size() );
            backend.Download( *state, finalState.data() );
            const double elapsed = std::chrono::duration< double, std::milli >(
                std::chrono::steady_clock::now() - begin ).count();
            if ( repeat >= warmup )
            {
                totalAdvanceMs += elapsed;
                totalKernelMs += stats.kernelMilliseconds;
                totalLaunches += stats.kernelLaunches;
                totalSyncs += stats.synchronizationCount;
            }
        }

        std::cout << std::fixed << std::setprecision( 6 );
        std::cout << "OneFLOW 1D Euler GPU-only sustained benchmark\n";
        std::cout << "nx=" << nx << " steps=" << steps
                  << " repeats=" << repeats << " warmup=" << warmup << "\n";
        std::cout << "device=" << properties.name
                  << " arch=" << properties.gcnArchName << "\n";
        std::cout << "gpu_only_total_ms=" << totalAdvanceMs
                  << " gpu_only_per_step_ms="
                  << totalAdvanceMs / repeats / steps << "\n";
        std::cout << "gpu_only_kernel_ms=" << totalKernelMs
                  << " gpu_only_kernel_ratio="
                  << totalKernelMs / totalAdvanceMs << "\n";
        std::cout << "gpu_only_kernel_launches=" << totalLaunches
                  << " gpu_only_syncs=" << totalSyncs << "\n";
        std::cout << "final_state_max_abs=" << MaxAbs( finalState )
                  << " final_state_checksum=" << Checksum( finalState ) << "\n";
        return 0;
    }
    catch ( const std::exception & error )
    {
        std::cerr << "Euler GPU-only benchmark failed: " << error.what() << "\n";
        return 1;
    }
}
