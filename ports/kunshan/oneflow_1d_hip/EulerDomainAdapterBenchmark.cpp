#include "OneDEulerDomainAdapter.h"

#include "AccelRuntime.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{

using Clock = std::chrono::steady_clock;

constexpr double Gamma = 1.4;
constexpr double DtScale = 0.08;

struct Result
{
    std::vector< double > state;
    double createMilliseconds = 0.0;
    double uploadMilliseconds = 0.0;
    double advanceMilliseconds = 0.0;
    double downloadMilliseconds = 0.0;
    oneflow_1d::EulerRunStats stats;
};

int ParsePositive( const char * text, const char * label )
{
    const int value = std::atoi( text );
    if ( value <= 0 )
    {
        throw std::invalid_argument(
            std::string( label ) + " must be positive" );
    }
    return value;
}

std::vector< double > InitialState( int cells )
{
    std::vector< double > state(
        oneflow_1d::EulerComponents * cells );
    const double dx = 1.0 / cells;
    for ( int cell = 0; cell < cells; ++ cell )
    {
        const double x = ( cell + 0.5 ) * dx;
        const double density =
            1.0 + 0.08 * std::sin( 6.283185307179586 * x );
        const double velocity =
            0.2 + 0.04 * std::cos( 6.283185307179586 * x );
        const double pressure =
            1.0 + 0.05 * std::sin( 12.566370614359172 * x );
        state[ cell ] = density;
        state[ cells + cell ] = density * velocity;
        state[ 2 * cells + cell ] =
            pressure / ( Gamma - 1.0 )
            + 0.5 * density * velocity * velocity;
    }
    return state;
}

template< class Duration >
double Milliseconds( Duration duration )
{
    return std::chrono::duration< double, std::milli >( duration ).count();
}

Result RunDirect(
    const oneflow_1d::EulerBackend & backend,
    const std::vector< double > & initial,
    int cells,
    int steps )
{
    const oneflow_1d::EulerProblem problem{
        cells,
        Gamma,
        DtScale / cells,
        1.0 / cells,
        oneflow_1d::EulerBoundary::Periodic,
        oneflow_1d::EulerMethod::Rusanov };

    const auto createBegin = Clock::now();
    auto state = backend.CreateState( problem );
    const auto createEnd = Clock::now();
    backend.Upload( *state, initial.data() );
    const auto uploadEnd = Clock::now();
    oneflow_1d::EulerRunStats stats;
    backend.Advance(
        *state,
        steps,
        { oneflow_1d::EulerRunMode::NoTrace, nullptr, &stats } );
    const auto advanceEnd = Clock::now();
    std::vector< double > finalState( initial.size() );
    backend.Download( *state, finalState.data() );
    const auto downloadEnd = Clock::now();

    return {
        std::move( finalState ),
        Milliseconds( createEnd - createBegin ),
        Milliseconds( uploadEnd - createEnd ),
        Milliseconds( advanceEnd - uploadEnd ),
        Milliseconds( downloadEnd - advanceEnd ),
        stats };
}

Result RunShared(
    const ONEFLOW::EulerDomainBackend & backend,
    const std::vector< double > & initial,
    int cells,
    int steps )
{
    ONEFLOW::EulerDomainProblem problem;
    problem.nCells = cells;
    problem.nGhostCells = 0;
    problem.nEquations = oneflow_1d::EulerComponents;
    problem.gamma = Gamma;
    problem.dt = DtScale / cells;
    problem.dx = 1.0 / cells;
    problem.boundary = ONEFLOW::EulerDomainBoundary::Periodic;
    const ONEFLOW::EulerDomainStateKey key{
        0, 0, 0, ONEFLOW::AccelBackendKind::HIP, 0 };

    const auto createBegin = Clock::now();
    auto state = backend.CreateState( problem, key );
    const auto createEnd = Clock::now();
    const ONEFLOW::EulerDomainConstFieldView input{
        cells, oneflow_1d::EulerComponents, initial.data() };
    backend.Upload( *state, input );
    const auto uploadEnd = Clock::now();
    oneflow_1d::EulerRunStats stats;
    ONEFLOW::EulerDomainRunOptions options;
    options.stats = &stats;
    backend.Advance( *state, steps, options );
    const auto advanceEnd = Clock::now();
    std::vector< double > finalState( initial.size() );
    ONEFLOW::EulerDomainFieldView output{
        cells, oneflow_1d::EulerComponents, finalState.data() };
    backend.Download( *state, output );
    const auto downloadEnd = Clock::now();

    return {
        std::move( finalState ),
        Milliseconds( createEnd - createBegin ),
        Milliseconds( uploadEnd - createEnd ),
        Milliseconds( advanceEnd - uploadEnd ),
        Milliseconds( downloadEnd - advanceEnd ),
        stats };
}

void Accumulate( Result & total, Result value )
{
    total.state = std::move( value.state );
    total.createMilliseconds += value.createMilliseconds;
    total.uploadMilliseconds += value.uploadMilliseconds;
    total.advanceMilliseconds += value.advanceMilliseconds;
    total.downloadMilliseconds += value.downloadMilliseconds;
    total.stats.kernelMilliseconds += value.stats.kernelMilliseconds;
    total.stats.kernelLaunches += value.stats.kernelLaunches;
    total.stats.synchronizationCount += value.stats.synchronizationCount;
}

double MaxAbs(
    const std::vector< double > & left,
    const std::vector< double > & right )
{
    if ( left.size() != right.size() )
    {
        throw std::invalid_argument( "state size mismatch" );
    }
    double result = 0.0;
    for ( std::size_t index = 0; index < left.size(); ++ index )
    {
        result = std::max(
            result, std::abs( left[ index ] - right[ index ] ) );
    }
    return result;
}

long double Checksum( const std::vector< double > & values )
{
    long double result = 0.0;
    for ( double value : values ) result += value;
    return result;
}

}

int main( int argc, char ** argv )
{
    ONEFLOW::InitializeAccelRuntime( 0, 1 );
    try
    {
        const int cells =
            argc > 1 ? ParsePositive( argv[ 1 ], "cells" ) : 1048576;
        const int steps =
            argc > 2 ? ParsePositive( argv[ 2 ], "steps" ) : 100;
        const int repeats =
            argc > 3 ? ParsePositive( argv[ 3 ], "repeats" ) : 3;
        const int warmup =
            argc > 4 ? ParsePositive( argv[ 4 ], "warmup" ) : 1;

        const std::vector< double > initial = InitialState( cells );
        oneflow_1d::HipEulerBackend directBackend;
        oneflow_1d::EulerDomainBackendAdapter sharedBackend(
            directBackend, ONEFLOW::AccelBackendKind::HIP, 0 );

        for ( int repeat = 0; repeat < warmup; ++ repeat )
        {
            static_cast< void >(
                RunDirect( directBackend, initial, cells, steps ) );
            static_cast< void >(
                RunShared( sharedBackend, initial, cells, steps ) );
        }

        Result directTotal;
        Result sharedTotal;
        for ( int repeat = 0; repeat < repeats; ++ repeat )
        {
            if ( repeat % 2 == 0 )
            {
                Accumulate(
                    directTotal,
                    RunDirect( directBackend, initial, cells, steps ) );
                Accumulate(
                    sharedTotal,
                    RunShared( sharedBackend, initial, cells, steps ) );
            }
            else
            {
                Accumulate(
                    sharedTotal,
                    RunShared( sharedBackend, initial, cells, steps ) );
                Accumulate(
                    directTotal,
                    RunDirect( directBackend, initial, cells, steps ) );
            }
        }

        const double maxAbs =
            MaxAbs( directTotal.state, sharedTotal.state );
        const long double directChecksum = Checksum( directTotal.state );
        const long double sharedChecksum = Checksum( sharedTotal.state );
        const bool statsMatch =
            directTotal.stats.kernelLaunches
                == sharedTotal.stats.kernelLaunches
            && directTotal.stats.synchronizationCount
                == sharedTotal.stats.synchronizationCount;

        std::cout << std::fixed << std::setprecision( 6 );
        std::cout
            << "OneFLOW 1D Euler shared adapter benchmark\n"
            << "cells=" << cells
            << " steps=" << steps
            << " repeats=" << repeats
            << " warmup=" << warmup << "\n"
            << "direct_hip_advance_ms="
            << directTotal.advanceMilliseconds
            << " shared_hip_advance_ms="
            << sharedTotal.advanceMilliseconds
            << " shared_over_direct="
            << sharedTotal.advanceMilliseconds
                / directTotal.advanceMilliseconds << "\n"
            << "direct_hip_kernel_ms="
            << directTotal.stats.kernelMilliseconds
            << " shared_hip_kernel_ms="
            << sharedTotal.stats.kernelMilliseconds << "\n"
            << "direct_hip_create_ms="
            << directTotal.createMilliseconds
            << " direct_hip_upload_ms="
            << directTotal.uploadMilliseconds
            << " direct_hip_download_ms="
            << directTotal.downloadMilliseconds << "\n"
            << "shared_hip_create_ms="
            << sharedTotal.createMilliseconds
            << " shared_hip_upload_ms="
            << sharedTotal.uploadMilliseconds
            << " shared_hip_download_ms="
            << sharedTotal.downloadMilliseconds << "\n"
            << "direct_hip_kernel_launches="
            << directTotal.stats.kernelLaunches
            << " shared_hip_kernel_launches="
            << sharedTotal.stats.kernelLaunches
            << " direct_hip_syncs="
            << directTotal.stats.synchronizationCount
            << " shared_hip_syncs="
            << sharedTotal.stats.synchronizationCount << "\n"
            << "final_max_abs_error=" << maxAbs << "\n"
            << std::setprecision( 17 )
            << "direct_checksum=" << directChecksum
            << " shared_checksum=" << sharedChecksum << "\n";

        const bool pass =
            sharedBackend.IsAccelerator()
            && std::string( sharedBackend.Name() ) == directBackend.Name()
            && statsMatch
            && maxAbs <= 1.0e-15
            && directChecksum == sharedChecksum;
        std::cout << "adapter_benchmark=" << ( pass ? "PASS" : "FAIL" )
                  << "\n";
        ONEFLOW::FinalizeAccelRuntime();
        return pass ? 0 : 1;
    }
    catch ( const std::exception & error )
    {
        std::cerr
            << "Euler domain adapter benchmark failed: "
            << error.what() << "\n";
        ONEFLOW::FinalizeAccelRuntime();
        return 1;
    }
}
