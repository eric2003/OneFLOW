#include <mpi.h>
#include <hip/hip_runtime.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

namespace
{

void CheckHip( hipError_t status, const char * operation )
{
    if ( status == hipSuccess ) return;
    std::fprintf(
        stderr,
        "GPU_AWARE_MPI_PROBE hip failure op=%s code=%d message=%s\n",
        operation,
        static_cast< int >( status ),
        hipGetErrorString( status ) );
    MPI_Abort( MPI_COMM_WORLD, 2 );
}

void CheckMpi( int status, const char * operation )
{
    if ( status == MPI_SUCCESS ) return;
    char message[ MPI_MAX_ERROR_STRING ] = {};
    int length = 0;
    MPI_Error_string( status, message, &length );
    std::fprintf(
        stderr,
        "GPU_AWARE_MPI_PROBE mpi failure op=%s code=%d message=%.*s\n",
        operation,
        status,
        length,
        message );
    MPI_Abort( MPI_COMM_WORLD, 3 );
}

} // namespace

int main( int argc, char ** argv )
{
    CheckMpi( MPI_Init( &argc, &argv ), "MPI_Init" );

    int rank = 0;
    int size = 0;
    CheckMpi( MPI_Comm_rank( MPI_COMM_WORLD, &rank ), "MPI_Comm_rank" );
    CheckMpi( MPI_Comm_size( MPI_COMM_WORLD, &size ), "MPI_Comm_size" );
    if ( size != 2 )
    {
        if ( rank == 0 )
            std::fprintf(
                stderr,
                "GPU_AWARE_MPI_PROBE requires exactly two MPI ranks\n" );
        MPI_Finalize();
        return 4;
    }

    int localRank = 0;
    MPI_Comm local = MPI_COMM_NULL;
    CheckMpi(
        MPI_Comm_split_type(
            MPI_COMM_WORLD,
            MPI_COMM_TYPE_SHARED,
            rank,
            MPI_INFO_NULL,
            &local ),
        "MPI_Comm_split_type" );
    CheckMpi( MPI_Comm_rank( local, &localRank ), "local rank" );

    int deviceCount = 0;
    CheckHip( hipGetDeviceCount( &deviceCount ), "hipGetDeviceCount" );
    if ( deviceCount <= 0 )
    {
        std::fprintf(
            stderr,
            "GPU_AWARE_MPI_PROBE rank=%d has no visible HIP device\n",
            rank );
        MPI_Abort( MPI_COMM_WORLD, 5 );
    }
    CheckHip( hipSetDevice( localRank % deviceCount ), "hipSetDevice" );

    hipDeviceProp_t properties {};
    CheckHip(
        hipGetDeviceProperties(
            &properties,
            localRank % deviceCount ),
        "hipGetDeviceProperties" );

    constexpr std::size_t valueCount = 1024 * 1024;
    constexpr std::size_t valueBytes = valueCount * sizeof( double );
    double * deviceSend = nullptr;
    double * deviceReceive = nullptr;
    CheckHip( hipMalloc( &deviceSend, valueBytes ), "hipMalloc send" );
    CheckHip(
        hipMalloc( &deviceReceive, valueBytes ),
        "hipMalloc receive" );

    const double sendValue = static_cast< double >( rank + 1 );
    CheckHip( hipMemset( deviceSend, 0, valueBytes ), "hipMemset" );
    CheckHip(
        hipMemcpy(
            deviceSend,
            &sendValue,
            sizeof( sendValue ),
            hipMemcpyHostToDevice ),
        "hipMemcpyHostToDevice" );
    CheckHip( hipDeviceSynchronize(), "hipDeviceSynchronize before MPI" );

    CheckMpi( MPI_Barrier( MPI_COMM_WORLD ), "MPI_Barrier" );
    const int peer = 1 - rank;
    const int tag = 9100;
    const int exchangeStatus = MPI_Sendrecv(
        deviceSend,
        static_cast< int >( valueCount ),
        MPI_DOUBLE,
        peer,
        tag,
        deviceReceive,
        static_cast< int >( valueCount ),
        MPI_DOUBLE,
        peer,
        tag,
        MPI_COMM_WORLD,
        MPI_STATUS_IGNORE );
    CheckMpi( exchangeStatus, "MPI_Sendrecv(device pointer)" );

    double received = 0.0;
    CheckHip(
        hipMemcpy(
            &received,
            deviceReceive,
            sizeof( received ),
            hipMemcpyDeviceToHost ),
        "hipMemcpyDeviceToHost" );
    const double expected = static_cast< double >( peer + 1 );
    const int localPass = received == expected ? 1 : 0;
    int globalPass = 0;
    CheckMpi(
        MPI_Allreduce(
            &localPass,
            &globalPass,
            1,
            MPI_INT,
            MPI_MIN,
            MPI_COMM_WORLD ),
        "MPI_Allreduce" );

    char hostName[MPI_MAX_PROCESSOR_NAME] = {};
    int hostNameLength = 0;
    MPI_Get_processor_name( hostName, &hostNameLength );
    std::printf(
        "GPU_AWARE_MPI_PROBE rank=%d host=%.*s local_rank=%d device=%d "
        "arch=%s received=%.17g expected=%.17g pass=%d\n",
        rank,
        hostNameLength,
        hostName,
        localRank,
        localRank % deviceCount,
        properties.gcnArchName,
        received,
        expected,
        globalPass );
    std::fflush( stdout );

    CheckHip( hipFree( deviceSend ), "hipFree send" );
    CheckHip( hipFree( deviceReceive ), "hipFree receive" );
    MPI_Comm_free( &local );
    MPI_Finalize();
    return globalPass == 1 ? 0 : 6;
}
