#include "OneDEulerMpi.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
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

double MaxScaledError(
    const std::vector< double > & reference,
    const std::vector< double > & candidate,
    long long & violations )
{
    if ( reference.size() != candidate.size() )
    {
        throw std::invalid_argument( "CPU/HIP state sizes differ" );
    }
    double maximum = 0.0;
    violations = 0;
    for ( std::size_t index = 0; index < reference.size(); ++ index )
    {
        const double a = reference[ index ];
        const double b = candidate[ index ];
        const double scale = std::max( 1.0, std::abs( a ) );
        if ( ! std::isfinite( a ) || ! std::isfinite( b ) )
        {
            ++ violations;
            continue;
        }
        const double error = std::abs( b - a );
        maximum = std::max( maximum, error / scale );
        if ( error > 1e-15 + 1e-15 * scale ) ++ violations;
    }
    return maximum;
}

bool Physical( const std::vector< double > & state, int nx )
{
    for ( int cell = 0; cell < nx; ++ cell )
    {
        const double rho = state[ cell ];
        const double momentum = state[ nx + cell ];
        const double energy = state[ 2 * nx + cell ];
        const double pressure = ( Gamma - 1.0 )
            * ( energy - 0.5 * momentum * momentum / rho );
        if ( ! std::isfinite( rho ) || ! std::isfinite( pressure )
             || rho <= 0.0 || pressure <= 0.0 )
        {
            return false;
        }
    }
    return true;
}

std::uint64_t ReduceHash( std::uint64_t local, MPI_Comm communicator )
{
    unsigned long long localValue = local;
    unsigned long long globalValue = 0;
    MPI_Reduce(
        &localValue, &globalValue, 1, MPI_UNSIGNED_LONG_LONG, MPI_BXOR,
        0, communicator );
    return static_cast< std::uint64_t >( globalValue );
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
            argc > 1 ? ParsePositive( argv[ 1 ], "global nx" ) : 65536;
        const int steps =
            argc > 2 ? ParsePositive( argv[ 2 ], "steps" ) : 20;

        int localNx = 0;
        int globalOffset = 0;
        Partition( globalNx, rank, ranks, localNx, globalOffset );
        const auto initial = InitialState( globalNx, localNx, globalOffset );
        const oneflow_1d::EulerMpiProblem problem{
            globalNx, localNx, globalOffset, steps, Gamma,
            DtScale / globalNx, 1.0 / globalNx, MPI_COMM_WORLD };

        MPI_Barrier( MPI_COMM_WORLD );
        const auto cpu = oneflow_1d::RunCpuEulerMpi( problem, initial );
        MPI_Barrier( MPI_COMM_WORLD );
        const auto hip = oneflow_1d::RunHipEulerMpi( problem, initial );
        MPI_Barrier( MPI_COMM_WORLD );

        long long localViolations = 0;
        const double localScaledError =
            MaxScaledError( cpu.state, hip.state, localViolations );
        long long globalViolations = 0;
        MPI_Reduce(
            &localViolations, &globalViolations, 1, MPI_LONG_LONG, MPI_SUM,
            0, MPI_COMM_WORLD );
        double globalScaledError = 0.0;
        MPI_Reduce(
            &localScaledError, &globalScaledError, 1, MPI_DOUBLE, MPI_MAX,
            0, MPI_COMM_WORLD );

        const int localPhysical =
            Physical( cpu.state, localNx ) && Physical( hip.state, localNx );
        int globalPhysical = 0;
        MPI_Reduce(
            &localPhysical, &globalPhysical, 1, MPI_INT, MPI_MIN, 0,
            MPI_COMM_WORLD );

        const std::uint64_t cpuHash = ReduceHash(
            oneflow_1d::EulerMpiStateHash(
                cpu.state, globalNx, localNx, globalOffset ),
            MPI_COMM_WORLD );
        const std::uint64_t hipHash = ReduceHash(
            oneflow_1d::EulerMpiStateHash(
                hip.state, globalNx, localNx, globalOffset ),
            MPI_COMM_WORLD );

        if ( rank == 0 )
        {
            std::cout << std::fixed << std::setprecision( 6 );
            std::cout << "OneFLOW 1D Euler MPI CPU/HIP regression\n";
            std::cout << "ranks=" << ranks
                      << " global_nx=" << globalNx
                      << " steps=" << steps << "\n";
            std::cout << "hip_device=" << hip.deviceName
                      << " arch=" << hip.architecture
                      << " local_rank=" << hip.localRank
                      << " local_ranks=" << hip.localRanks
                      << " device_index=" << hip.deviceIndex
                      << " visible_devices=" << hip.visibleDeviceCount << "\n";
            std::cout << "max_scaled_error=" << globalScaledError
                      << " violations=" << globalViolations
                      << " physical=" << globalPhysical << "\n";
            std::cout << std::hex
                      << "cpu_hash=0x" << cpuHash
                      << " hip_hash=0x" << hipHash
                      << std::dec << "\n";
            const bool pass = globalViolations == 0 && globalPhysical != 0
                && cpuHash == hipHash;
            std::cout << "result=" << ( pass ? "PASS" : "FAIL" ) << "\n";
            MPI_Finalize();
            return pass ? 0 : 1;
        }

        MPI_Finalize();
        return 0;
    }
    catch ( const std::exception & error )
    {
        std::cerr << "rank " << rank
                  << " Euler MPI regression failed: " << error.what() << "\n";
        MPI_Abort( MPI_COMM_WORLD, 1 );
        return 1;
    }
}
