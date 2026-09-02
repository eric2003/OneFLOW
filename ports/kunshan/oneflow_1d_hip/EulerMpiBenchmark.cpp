#include "OneDEulerMpi.h"

#include <mpi.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

constexpr double Gamma = 1.4;
constexpr double DtScale = 0.05;
using Clock = std::chrono::steady_clock;

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

void Partition(
    int globalNx, int rank, int ranks, int & localNx, int & globalOffset )
{
    const int quotient = globalNx / ranks;
    const int remainder = globalNx % ranks;
    localNx = quotient + ( rank < remainder ? 1 : 0 );
    globalOffset = rank * quotient + std::min( rank, remainder );
    if ( localNx <= 0 )
    {
        throw std::invalid_argument( "global nx is smaller than MPI rank count" );
    }
}

std::vector< double > InitialState(
    int globalNx, int localNx, int globalOffset )
{
    std::vector< double > state( 3 * localNx );
    for ( int localCell = 0; localCell < localNx; ++ localCell )
    {
        const int globalCell = globalOffset + localCell;
        const double x = ( globalCell + 0.5 ) / globalNx;
        const double rho = 1.0 + 0.1 * std::sin( 2.0 * M_PI * x );
        const double velocity = 0.2 + 0.05 * std::cos( 2.0 * M_PI * x );
        const double pressure = 1.0 + 0.08 * std::sin( 4.0 * M_PI * x );
        state[ localCell ] = rho;
        state[ localNx + localCell ] = rho * velocity;
        state[ 2 * localNx + localCell ] =
            pressure / ( Gamma - 1.0 ) + 0.5 * rho * velocity * velocity;
    }
    return state;
}

double ReduceMax( double local, MPI_Comm communicator )
{
    double maximum = 0.0;
    MPI_Reduce(
        &local, &maximum, 1, MPI_DOUBLE, MPI_MAX, 0, communicator );
    return maximum;
}

long long ReduceSum( int local, MPI_Comm communicator )
{
    const long long value = local;
    long long total = 0;
    MPI_Reduce(
        &value, &total, 1, MPI_LONG_LONG, MPI_SUM, 0, communicator );
    return total;
}

std::uint64_t ReduceHash( std::uint64_t local, MPI_Comm communicator )
{
    unsigned long long value = local;
    unsigned long long result = 0;
    MPI_Reduce(
        &value, &result, 1, MPI_UNSIGNED_LONG_LONG, MPI_BXOR,
        0, communicator );
    return static_cast< std::uint64_t >( result );
}

void PhysicalMinimums(
    const std::vector< double > & state, int nx,
    double & minimumRho, double & minimumPressure )
{
    minimumRho = std::numeric_limits< double >::infinity();
    minimumPressure = std::numeric_limits< double >::infinity();
    for ( int cell = 0; cell < nx; ++ cell )
    {
        const double rho = state[ cell ];
        const double momentum = state[ nx + cell ];
        const double energy = state[ 2 * nx + cell ];
        const double pressure = ( Gamma - 1.0 )
            * ( energy - 0.5 * momentum * momentum / rho );
        minimumRho = std::min( minimumRho, rho );
        minimumPressure = std::min( minimumPressure, pressure );
    }
    return;
}

oneflow_1d::EulerMpiResult Run(
    const oneflow_1d::EulerMpiProblem & problem,
    const std::vector< double > & initial )
{
#ifdef ONEFLOW_1D_MPI_USE_HIP
    return oneflow_1d::RunHipEulerMpi( problem, initial );
#else
    return oneflow_1d::RunCpuEulerMpi( problem, initial );
#endif
}

const char * BackendName()
{
#ifdef ONEFLOW_1D_MPI_USE_HIP
    return "HIP_MPI_SHARED_DEVICE";
#else
    return "CPU_MPI";
#endif
}

} // namespace

int main( int argc, char ** argv )
{
    MPI_Init( &argc, &argv );
    int rank = 0;
    int ranks = 0;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &ranks );

    try
    {
        const int globalNx =
            argc > 1 ? ParsePositive( argv[ 1 ], "global nx" ) : 1048576;
        const int steps =
            argc > 2 ? ParsePositive( argv[ 2 ], "steps" ) : 100;
        const int repeats =
            argc > 3 ? ParsePositive( argv[ 3 ], "repeats" ) : 3;
        const int warmup =
            argc > 4 ? ParsePositive( argv[ 4 ], "warmup" ) : 1;

        int localNx = 0;
        int globalOffset = 0;
        Partition( globalNx, rank, ranks, localNx, globalOffset );
        const auto initial = InitialState( globalNx, localNx, globalOffset );
        const oneflow_1d::EulerMpiProblem problem{
            globalNx, localNx, globalOffset, steps, Gamma,
            DtScale / globalNx, 1.0 / globalNx, MPI_COMM_WORLD };

        MPI_Comm localCommunicator = MPI_COMM_NULL;
        MPI_Comm_split_type(
            MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL,
            &localCommunicator );
        int localRanks = 0;
        MPI_Comm_size( localCommunicator, &localRanks );
        MPI_Comm_free( &localCommunicator );

        for ( int repeat = 0; repeat < warmup; ++ repeat )
        {
            MPI_Barrier( MPI_COMM_WORLD );
            static_cast< void >( Run( problem, initial ) );
            MPI_Barrier( MPI_COMM_WORLD );
        }

        double totalMilliseconds = 0.0;
        double createMilliseconds = 0.0;
        double uploadMilliseconds = 0.0;
        double advanceMilliseconds = 0.0;
        double downloadMilliseconds = 0.0;
        double mpiMilliseconds = 0.0;
        double computeMilliseconds = 0.0;
        double kernelMilliseconds = 0.0;
        long long kernelLaunches = 0;
        long long haloExchanges = 0;
        long long deviceSynchronizations = 0;
        oneflow_1d::EulerMpiResult finalResult;

        for ( int repeat = 0; repeat < repeats; ++ repeat )
        {
            MPI_Barrier( MPI_COMM_WORLD );
            const auto begin = Clock::now();
            auto result = Run( problem, initial );
            MPI_Barrier( MPI_COMM_WORLD );
            const double localTotal =
                std::chrono::duration< double, std::milli >(
                    Clock::now() - begin ).count();

            totalMilliseconds += ReduceMax( localTotal, MPI_COMM_WORLD );
            createMilliseconds += ReduceMax(
                result.timing.createMilliseconds, MPI_COMM_WORLD );
            uploadMilliseconds += ReduceMax(
                result.timing.uploadMilliseconds, MPI_COMM_WORLD );
            advanceMilliseconds += ReduceMax(
                result.timing.advanceMilliseconds, MPI_COMM_WORLD );
            downloadMilliseconds += ReduceMax(
                result.timing.downloadMilliseconds, MPI_COMM_WORLD );
            mpiMilliseconds += ReduceMax(
                result.timing.mpiMilliseconds, MPI_COMM_WORLD );
            computeMilliseconds += ReduceMax(
                result.timing.computeMilliseconds, MPI_COMM_WORLD );
            kernelMilliseconds += ReduceMax(
                result.timing.kernelMilliseconds, MPI_COMM_WORLD );
            kernelLaunches += ReduceSum(
                result.timing.kernelLaunches, MPI_COMM_WORLD );
            haloExchanges += ReduceSum(
                result.timing.haloExchanges, MPI_COMM_WORLD );
            deviceSynchronizations += ReduceSum(
                result.timing.deviceSynchronizations, MPI_COMM_WORLD );
            finalResult = std::move( result );
        }

        const std::uint64_t localHash = oneflow_1d::EulerMpiStateHash(
            finalResult.state, globalNx, localNx, globalOffset );
        const std::uint64_t globalHash =
            ReduceHash( localHash, MPI_COMM_WORLD );
        double localMinimumRho = 0.0;
        double localMinimumPressure = 0.0;
        PhysicalMinimums(
            finalResult.state, localNx,
            localMinimumRho, localMinimumPressure );
        double globalMinimumRho = 0.0;
        double globalMinimumPressure = 0.0;
        MPI_Reduce(
            &localMinimumRho, &globalMinimumRho, 1, MPI_DOUBLE, MPI_MIN,
            0, MPI_COMM_WORLD );
        MPI_Reduce(
            &localMinimumPressure, &globalMinimumPressure, 1, MPI_DOUBLE,
            MPI_MIN, 0, MPI_COMM_WORLD );

        if ( rank == 0 )
        {
            char processorName[ MPI_MAX_PROCESSOR_NAME ];
            int processorNameLength = 0;
            MPI_Get_processor_name( processorName, &processorNameLength );
            std::cout << std::fixed << std::setprecision( 6 );
            std::cout << "OneFLOW 1D Euler MPI backend benchmark\n";
            std::cout << "backend=" << BackendName()
                      << " ranks=" << ranks
                      << " local_ranks=" << localRanks
                      << " host="
                      << std::string( processorName, processorNameLength )
                      << "\n";
            std::cout << "global_nx=" << globalNx
                      << " steps=" << steps
                      << " repeats=" << repeats
                      << " warmup=" << warmup << "\n";
            if ( ! finalResult.deviceName.empty() )
            {
                std::cout << "device=" << finalResult.deviceName
                          << " arch=" << finalResult.architecture
                          << " local_rank=" << finalResult.localRank
                          << " local_ranks=" << finalResult.localRanks
                          << " device_index=" << finalResult.deviceIndex
                          << " visible_devices="
                          << finalResult.visibleDeviceCount << "\n";
            }
            std::cout << "lifecycle_max_ms=" << totalMilliseconds
                      << " lifecycle_per_step_ms="
                      << totalMilliseconds / repeats / steps << "\n";
            std::cout << "create_max_ms=" << createMilliseconds
                      << " upload_max_ms=" << uploadMilliseconds
                      << " advance_max_ms=" << advanceMilliseconds
                      << " download_max_ms=" << downloadMilliseconds << "\n";
            std::cout << "mpi_max_ms=" << mpiMilliseconds
                      << " compute_max_ms=" << computeMilliseconds
                      << " kernel_max_ms=" << kernelMilliseconds << "\n";
            std::cout << "kernel_launches_sum=" << kernelLaunches
                      << " halo_exchanges_sum=" << haloExchanges
                      << " device_syncs_sum=" << deviceSynchronizations
                      << "\n";
            std::cout << "final_hash=0x" << std::hex << globalHash
                      << std::dec
                      << " min_rho=" << globalMinimumRho
                      << " min_pressure=" << globalMinimumPressure << "\n";
        }

        MPI_Finalize();
        return 0;
    }
    catch ( const std::exception & error )
    {
        std::cerr << "rank " << rank
                  << " Euler MPI benchmark failed: " << error.what() << "\n";
        MPI_Abort( MPI_COMM_WORLD, 1 );
        return 1;
    }
}
