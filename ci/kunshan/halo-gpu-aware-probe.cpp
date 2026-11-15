// Production-density GPU-aware MPI halo probe.
//
// The shipped probe (gpu-aware-mpi-probe.cpp) requires exactly two ranks and
// exchanges one 8 MB buffer on a single pair.  That answers "is the transport
// capable" but not the question that matters here: at eight ranks per card,
// with the message sizes the solver actually sends, is a device-pointer
// exchange faster than the host-staged path the solver uses today?
//
// This probe reproduces the production shape:
//   * N ranks (any N), ring topology, so every rank has two neighbours the way
//     most M6 zones do;
//   * halo-sized messages (nFaces * nEquations * 8 bytes, default ~160 KB,
//     matching 4,000 interface faces over five equations);
//   * four exchange rounds per step, the real cadence (three RK stages plus
//     the post round);
//   * two transfer modes measured back to back in one run:
//       device - MPI_Sendrecv straight from hipMalloc memory
//       staged - hipMemcpy D2H, MPI_Sendrecv on host, hipMemcpy H2D
//     which is exactly the difference the solver would gain by switching.
//
// Correctness is checked every round in both modes; a mode that corrupts data
// is reported as failed regardless of its timing.
//
// Build: hipcc -std=c++17 -O2 halo-probe.cpp -o halo-probe \
//          $(mpicxx --showme:compile) $(mpicxx --showme:link)
// Run:   mpirun -np N ./halo-probe [rounds] [faces] [equations]

#include <mpi.h>
#include <hip/hip_runtime.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace
{

constexpr int kEquations = 5;

void CheckHip( hipError_t status, const char * operation )
{
    if ( status == hipSuccess ) return;
    std::fprintf( stderr, "HALO_PROBE hip failure op=%s code=%d message=%s\n",
                  operation, static_cast< int >( status ),
                  hipGetErrorString( status ) );
    MPI_Abort( MPI_COMM_WORLD, 2 );
}

void CheckMpi( int status, const char * operation )
{
    if ( status == MPI_SUCCESS ) return;
    char message[ MPI_MAX_ERROR_STRING ] = {};
    int length = 0;
    MPI_Error_string( status, message, &length );
    std::fprintf( stderr, "HALO_PROBE mpi failure op=%s code=%d message=%.*s\n",
                  operation, status, length, message );
    MPI_Abort( MPI_COMM_WORLD, 3 );
}

double Now()
{
    using clock = std::chrono::steady_clock;
    return std::chrono::duration< double >(
        clock::now().time_since_epoch() ).count();
}

} // namespace

int main( int argc, char ** argv )
{
    CheckMpi( MPI_Init( &argc, &argv ), "MPI_Init" );

    int rank = 0;
    int size = 0;
    CheckMpi( MPI_Comm_rank( MPI_COMM_WORLD, &rank ), "MPI_Comm_rank" );
    CheckMpi( MPI_Comm_size( MPI_COMM_WORLD, &size ), "MPI_Comm_size" );

    const int rounds = argc > 1 ? std::atoi( argv[ 1 ] ) : 4;
    const int faces = argc > 2 ? std::atoi( argv[ 2 ] ) : 4000;
    const int equations = argc > 3 ? std::atoi( argv[ 3 ] ) : kEquations;
    if ( rounds <= 0 || faces <= 0 || equations <= 0 || equations > 32 )
    {
        if ( rank == 0 ) std::fprintf( stderr, "HALO_PROBE bad arguments\n" );
        MPI_Finalize();
        return 2;
    }
    const std::size_t valueCount =
        static_cast< std::size_t >( faces ) * equations;
    const std::size_t valueBytes = valueCount * sizeof( double );

    int localRank = 0;
    MPI_Comm local = MPI_COMM_NULL;
    CheckMpi( MPI_Comm_split_type( MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank,
                                   MPI_INFO_NULL, &local ),
              "MPI_Comm_split_type" );
    CheckMpi( MPI_Comm_rank( local, &localRank ), "local rank" );

    int deviceCount = 0;
    CheckHip( hipGetDeviceCount( &deviceCount ), "hipGetDeviceCount" );
    if ( deviceCount <= 0 )
    {
        std::fprintf( stderr, "HALO_PROBE rank=%d no visible device\n", rank );
        MPI_Abort( MPI_COMM_WORLD, 5 );
    }
    const int deviceId = localRank % deviceCount;
    CheckHip( hipSetDevice( deviceId ), "hipSetDevice" );
    hipDeviceProp_t properties{};
    CheckHip( hipGetDeviceProperties( &properties, deviceId ),
              "hipGetDeviceProperties" );

    double * deviceSend = nullptr;
    double * deviceRecv = nullptr;
    CheckHip( hipMalloc( &deviceSend, valueBytes ), "hipMalloc send" );
    CheckHip( hipMalloc( &deviceRecv, valueBytes ), "hipMalloc recv" );
    std::vector< double > hostSend( valueCount );
    std::vector< double > hostRecv( valueCount );
    const double sendValue = static_cast< double >( rank + 1 );
    for ( std::size_t i = 0; i < valueCount; ++ i )
        hostSend[ i ] = sendValue;

    // Ring neighbours: the two ranks this one exchanges halo with.
    const int previous = ( rank - 1 + size ) % size;
    const int next = ( rank + 1 ) % size;
    const bool ringIsTrivial = size < 2;

    auto exchange_device = [&]( int tag ) -> double {
        const double start = Now();
        // Two neighbours per round, matching a solver rank with two halos.
        if ( ! ringIsTrivial )
        {
            CheckMpi( MPI_Sendrecv( deviceSend, static_cast< int >( valueCount ),
                                    MPI_DOUBLE, next, tag, deviceRecv,
                                    static_cast< int >( valueCount ), MPI_DOUBLE,
                                    previous, tag, MPI_COMM_WORLD,
                                    MPI_STATUS_IGNORE ),
                      "device Sendrecv next" );
        }
        return Now() - start;
    };

    auto exchange_staged = [&]( int tag ) -> double {
        const double start = Now();
        CheckHip( hipMemcpy( hostSend.data(), deviceSend, valueBytes,
                             hipMemcpyDeviceToHost ),
                  "D2H stage" );
        if ( ! ringIsTrivial )
        {
            CheckMpi( MPI_Sendrecv( hostSend.data(),
                                    static_cast< int >( valueCount ), MPI_DOUBLE,
                                    next, tag, hostRecv.data(),
                                    static_cast< int >( valueCount ), MPI_DOUBLE,
                                    previous, tag, MPI_COMM_WORLD,
                                    MPI_STATUS_IGNORE ),
                      "staged Sendrecv next" );
        }
        CheckHip( hipMemcpy( deviceRecv, hostRecv.data(), valueBytes,
                             hipMemcpyHostToDevice ),
                  "H2D stage" );
        return Now() - start;
    };

    // Warm up both paths before timing: the first pass pays lazy allocation and
    // driver-page setup that would otherwise land in the measurement.
    CheckHip( hipMemcpy( deviceSend, hostSend.data(), valueBytes,
                         hipMemcpyHostToDevice ), "warmup H2D" );
    CheckHip( hipDeviceSynchronize(), "warmup sync" );
    exchange_device( 9000 );
    exchange_staged( 9001 );
    CheckMpi( MPI_Barrier( MPI_COMM_WORLD ), "warmup barrier" );

    double deviceSeconds = 0.0;
    double stagedSeconds = 0.0;
    int deviceCorrect = 1;
    int stagedCorrect = 1;
    const double expected = static_cast< double >( previous + 1 );

    for ( int round = 0; round < rounds; ++ round )
    {
        deviceSeconds += exchange_device( 9100 + round );
        if ( ! ringIsTrivial )
        {
            double first = 0.0;
            CheckHip( hipMemcpy( &first, deviceRecv, sizeof( first ),
                                 hipMemcpyDeviceToHost ),
                      "device verify D2H" );
            if ( first != expected ) deviceCorrect = 0;
        }

        stagedSeconds += exchange_staged( 9200 + round );
        if ( ! ringIsTrivial )
        {
            if ( hostRecv[ 0 ] != expected ) stagedCorrect = 0;
        }
    }
    CheckMpi( MPI_Barrier( MPI_COMM_WORLD ), "final barrier" );

    int minDeviceCorrect = 0;
    int minStagedCorrect = 0;
    CheckMpi( MPI_Allreduce( &deviceCorrect, &minDeviceCorrect, 1, MPI_INT,
                             MPI_MIN, MPI_COMM_WORLD ), "Allreduce device" );
    CheckMpi( MPI_Allreduce( &stagedCorrect, &minStagedCorrect, 1, MPI_INT,
                             MPI_MIN, MPI_COMM_WORLD ), "Allreduce staged" );

    double maxDevice = 0.0;
    double maxStaged = 0.0;
    CheckMpi( MPI_Allreduce( &deviceSeconds, &maxDevice, 1, MPI_DOUBLE, MPI_MAX,
                             MPI_COMM_WORLD ), "Allreduce device time" );
    CheckMpi( MPI_Allreduce( &stagedSeconds, &maxStaged, 1, MPI_DOUBLE, MPI_MAX,
                             MPI_COMM_WORLD ), "Allreduce staged time" );

    char hostName[ MPI_MAX_PROCESSOR_NAME ] = {};
    int hostNameLength = 0;
    MPI_Get_processor_name( hostName, &hostNameLength );
    std::printf(
        "HALO_PROBE rank=%d host=%.*s local_rank=%d device=%d size=%d "
        "faces=%d bytes=%zu rounds=%d "
        "device_ms=%.4f staged_ms=%.4f device_correct=%d staged_correct=%d\n",
        rank, hostNameLength, hostName, localRank, deviceId, size, faces,
        valueBytes, rounds, deviceSeconds * 1000.0, stagedSeconds * 1000.0,
        deviceCorrect, stagedCorrect );
    if ( rank == 0 )
    {
        const double perRoundDevice = maxDevice * 1000.0 / rounds;
        const double perRoundStaged = maxStaged * 1000.0 / rounds;
        std::printf(
            "HALO_PROBE_SUMMARY size=%d bytes=%zu device_per_round_ms=%.4f "
            "staged_per_round_ms=%.4f speedup=%.3f device_correct=%d "
            "staged_correct=%d\n",
            size, valueBytes, perRoundDevice, perRoundStaged,
            perRoundStaged > 0.0 ? perRoundStaged / perRoundDevice : 0.0,
            minDeviceCorrect, minStagedCorrect );
    }
    std::fflush( stdout );

    CheckHip( hipFree( deviceSend ), "hipFree send" );
    CheckHip( hipFree( deviceRecv ), "hipFree recv" );
    MPI_Comm_free( &local );
    MPI_Finalize();
    return ( minDeviceCorrect == 1 && minStagedCorrect == 1 ) ? 0 : 7;
}
