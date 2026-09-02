#include "OneDEulerMpi.h"

#include <hip/hip_runtime.h>

#include <chrono>
#include <cmath>
#include <cstdlib>
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

void CheckProblem(
    const EulerMpiProblem & problem,
    const std::vector< double > & initialState )
{
    if ( problem.communicator == MPI_COMM_NULL || problem.globalNx <= 0
         || problem.localNx <= 0 || problem.steps <= 0
         || problem.gamma <= 1.0 || problem.dt <= 0.0 || problem.dx <= 0.0 )
    {
        throw std::invalid_argument( "invalid HIP MPI Euler problem" );
    }
    if ( initialState.size()
         != static_cast< std::size_t >( 3 * problem.localNx ) )
    {
        throw std::invalid_argument( "invalid HIP MPI Euler initial state size" );
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
    const double state[ 3 ], double gamma, double flux[ 3 ],
    double & waveSpeed )
{
    const double rho = state[ 0 ];
    const double momentum = state[ 1 ];
    const double energy = state[ 2 ];
    const double velocity = momentum / rho;
    const double pressure = ( gamma - 1.0 )
        * ( energy - 0.5 * momentum * momentum / rho );
    flux[ 0 ] = momentum;
    flux[ 1 ] = momentum * velocity + pressure;
    flux[ 2 ] = ( energy + pressure ) * velocity;
    waveSpeed = fabs( velocity ) + sqrt( gamma * pressure / rho );
    return;
}

__global__ void PackBoundaryKernel(
    const double * current, int nx, double * boundary )
{
    const int component = threadIdx.x;
    if ( component < 3 )
    {
        boundary[ component ] = current[ CellIndex( component, 0, nx ) ];
        boundary[ 3 + component ] =
            current[ CellIndex( component, nx - 1, nx ) ];
    }
    return;
}

__global__ void FaceKernel(
    const double * current, const double * halos, int nx, double gamma,
    double * flux )
{
    const int face = blockIdx.x * blockDim.x + threadIdx.x;
    if ( face > nx ) return;

    double leftState[ 3 ];
    double rightState[ 3 ];
    for ( int component = 0; component < 3; ++ component )
    {
        leftState[ component ] = face == 0
            ? halos[ component ]
            : current[ CellIndex( component, face - 1, nx ) ];
        rightState[ component ] = face == nx
            ? halos[ 3 + component ]
            : current[ CellIndex( component, face, nx ) ];
    }

    double leftFlux[ 3 ];
    double rightFlux[ 3 ];
    double leftSpeed = 0.0;
    double rightSpeed = 0.0;
    PhysicalFlux( leftState, gamma, leftFlux, leftSpeed );
    PhysicalFlux( rightState, gamma, rightFlux, rightSpeed );
    const double speed = leftSpeed > rightSpeed ? leftSpeed : rightSpeed;
    for ( int component = 0; component < 3; ++ component )
    {
        flux[ FaceIndex( component, face, nx ) ] =
            0.5 * ( leftFlux[ component ] + rightFlux[ component ] )
            - 0.5 * speed
                * ( rightState[ component ] - leftState[ component ] );
    }
    return;
}

__global__ void ResidualKernel(
    const double * flux, int nx, double inverseDx, double * residual )
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    const int count = 3 * nx;
    if ( index >= count ) return;
    const int component = index / nx;
    const int cell = index - component * nx;
    residual[ index ] = -(
        flux[ FaceIndex( component, cell + 1, nx ) ]
        - flux[ FaceIndex( component, cell, nx ) ] ) * inverseDx;
    return;
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
    return;
}

struct DeviceState
{
    double * base = nullptr;
    double * work = nullptr;
    double * scratch = nullptr;
    double * flux = nullptr;
    double * residual = nullptr;
    double * boundary = nullptr;
    double * halos = nullptr;
    hipStream_t stream = nullptr;
    hipEvent_t kernelBegin = nullptr;
    hipEvent_t kernelEnd = nullptr;

    void Release() noexcept
    {
        if ( stream != nullptr )
        {
            static_cast< void >( hipStreamSynchronize( stream ) );
        }
        for ( double * pointer :
              { halos, boundary, residual, flux, scratch, work, base } )
        {
            if ( pointer != nullptr ) static_cast< void >( hipFree( pointer ) );
        }
        if ( kernelEnd != nullptr )
        {
            static_cast< void >( hipEventDestroy( kernelEnd ) );
        }
        if ( kernelBegin != nullptr )
        {
            static_cast< void >( hipEventDestroy( kernelBegin ) );
        }
        if ( stream != nullptr )
        {
            static_cast< void >( hipStreamDestroy( stream ) );
        }
        base = work = scratch = flux = residual = boundary = halos = nullptr;
        stream = nullptr;
        kernelBegin = kernelEnd = nullptr;
        return;
    }

    ~DeviceState()
    {
        Release();
    }
};

} // namespace

EulerMpiResult RunHipEulerMpi(
    const EulerMpiProblem & problem,
    const std::vector< double > & initialState )
{
    CheckProblem( problem, initialState );
    EulerMpiResult result;
    const int nx = problem.localNx;
    const int cellValues = 3 * nx;
    const int faceValues = 3 * ( nx + 1 );
    const std::size_t cellBytes =
        static_cast< std::size_t >( cellValues ) * sizeof( double );
    const std::size_t faceBytes =
        static_cast< std::size_t >( faceValues ) * sizeof( double );
    const std::size_t haloBytes = 6 * sizeof( double );

    int rank = 0;
    int ranks = 0;
    MPI_Comm_rank( problem.communicator, &rank );
    MPI_Comm_size( problem.communicator, &ranks );
    const int leftRank = ( rank + ranks - 1 ) % ranks;
    const int rightRank = ( rank + 1 ) % ranks;

    MPI_Comm localCommunicator = MPI_COMM_NULL;
    MPI_Comm_split_type(
        problem.communicator, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL,
        &localCommunicator );
    int localRank = 0;
    int localRanks = 0;
    MPI_Comm_rank( localCommunicator, &localRank );
    MPI_Comm_size( localCommunicator, &localRanks );
    MPI_Comm_free( &localCommunicator );

    int visibleDeviceCount = 0;
    Check( hipGetDeviceCount( &visibleDeviceCount ), "hipGetDeviceCount" );
    if ( visibleDeviceCount <= 0 )
    {
        throw std::runtime_error( "no visible HIP device for MPI rank" );
    }
    const char * requireOneRankPerDevice =
        std::getenv( "ONEFLOW_MPI_REQUIRE_ONE_RANK_PER_DEVICE" );
    if ( requireOneRankPerDevice != nullptr
         && std::string( requireOneRankPerDevice ) == "1"
         && localRanks > visibleDeviceCount )
    {
        throw std::runtime_error(
            "more local MPI ranks than visible HIP devices" );
    }
    const int deviceIndex = localRank % visibleDeviceCount;

    const auto createBegin = Clock::now();
    Check( hipSetDevice( deviceIndex ), "hipSetDevice" );
    hipDeviceProp_t properties {};
    Check(
        hipGetDeviceProperties( &properties, deviceIndex ),
        "hipGetDeviceProperties" );
    result.deviceName = properties.name;
    result.architecture = properties.gcnArchName;
    result.localRank = localRank;
    result.localRanks = localRanks;
    result.deviceIndex = deviceIndex;
    result.visibleDeviceCount = visibleDeviceCount;

    DeviceState device;
    try
    {
        Check(
            hipStreamCreateWithFlags( &device.stream, hipStreamNonBlocking ),
            "hipStreamCreateWithFlags" );
        Check( hipEventCreate( &device.kernelBegin ), "hipEventCreate begin" );
        Check( hipEventCreate( &device.kernelEnd ), "hipEventCreate end" );
        Check( hipMalloc( &device.base, cellBytes ), "hipMalloc base" );
        Check( hipMalloc( &device.work, cellBytes ), "hipMalloc work" );
        Check( hipMalloc( &device.scratch, cellBytes ), "hipMalloc scratch" );
        Check( hipMalloc( &device.flux, faceBytes ), "hipMalloc flux" );
        Check( hipMalloc( &device.residual, cellBytes ), "hipMalloc residual" );
        Check( hipMalloc( &device.boundary, haloBytes ), "hipMalloc boundary" );
        Check( hipMalloc( &device.halos, haloBytes ), "hipMalloc halos" );
    }
    catch ( ... )
    {
        device.Release();
        throw;
    }
    result.timing.createMilliseconds =
        Milliseconds( createBegin, Clock::now() );

    const auto uploadBegin = Clock::now();
    Check(
        hipMemcpyAsync(
            device.base, initialState.data(), cellBytes,
            hipMemcpyHostToDevice, device.stream ),
        "initial H2D" );
    Check( hipStreamSynchronize( device.stream ), "initial H2D sync" );
    ++ result.timing.deviceSynchronizations;
    result.timing.uploadMilliseconds =
        Milliseconds( uploadBegin, Clock::now() );

    const int threads = 256;
    const int faceBlocks = ( nx + 1 + threads - 1 ) / threads;
    const int cellBlocks = ( cellValues + threads - 1 ) / threads;
    const double inverseDx = 1.0 / problem.dx;
    double boundaryHost[ 6 ];
    double haloHost[ 6 ];

    const auto advanceBegin = Clock::now();
    for ( int step = 0; step < problem.steps; ++ step )
    {
        for ( int stage = 0; stage < 3; ++ stage )
        {
            double * current = stage == 0 ? device.base
                : ( stage == 1 ? device.work : device.scratch );
            double * next = stage == 0 ? device.work
                : ( stage == 1 ? device.scratch : device.work );

            hipLaunchKernelGGL(
                PackBoundaryKernel, dim3( 1 ), dim3( 3 ), 0, device.stream,
                current, nx, device.boundary );
            Check( hipGetLastError(), "PackBoundaryKernel launch" );
            ++ result.timing.kernelLaunches;
            const auto boundaryBegin = Clock::now();
            Check(
                hipMemcpyAsync(
                    boundaryHost, device.boundary, haloBytes,
                    hipMemcpyDeviceToHost, device.stream ),
                "boundary D2H" );
            Check(
                hipStreamSynchronize( device.stream ), "boundary D2H sync" );
            ++ result.timing.deviceSynchronizations;
            result.timing.computeMilliseconds +=
                Milliseconds( boundaryBegin, Clock::now() );

            const auto mpiBegin = Clock::now();
            MPI_Sendrecv(
                boundaryHost, 3, MPI_DOUBLE, leftRank, 100,
                haloHost + 3, 3, MPI_DOUBLE, rightRank, 100,
                problem.communicator, MPI_STATUS_IGNORE );
            MPI_Sendrecv(
                boundaryHost + 3, 3, MPI_DOUBLE, rightRank, 101,
                haloHost, 3, MPI_DOUBLE, leftRank, 101,
                problem.communicator, MPI_STATUS_IGNORE );
            result.timing.mpiMilliseconds +=
                Milliseconds( mpiBegin, Clock::now() );
            ++ result.timing.haloExchanges;

            Check(
                hipMemcpyAsync(
                    device.halos, haloHost, haloBytes,
                    hipMemcpyHostToDevice, device.stream ),
                "halo H2D" );
            Check(
                hipEventRecord( device.kernelBegin, device.stream ),
                "kernel begin event" );
            hipLaunchKernelGGL(
                FaceKernel, dim3( faceBlocks ), dim3( threads ), 0,
                device.stream, current, device.halos, nx, problem.gamma,
                device.flux );
            Check( hipGetLastError(), "FaceKernel launch" );
            hipLaunchKernelGGL(
                ResidualKernel, dim3( cellBlocks ), dim3( threads ), 0,
                device.stream, device.flux, nx, inverseDx, device.residual );
            Check( hipGetLastError(), "ResidualKernel launch" );
            const double baseWeight = stage == 1 ? 0.75
                : ( stage == 2 ? 1.0 / 3.0 : 0.0 );
            const double updateWeight = stage == 0 ? 1.0
                : ( stage == 1 ? 0.25 : 2.0 / 3.0 );
            hipLaunchKernelGGL(
                UpdateKernel, dim3( cellBlocks ), dim3( threads ), 0,
                device.stream, device.base, current, device.residual,
                cellValues, problem.dt, baseWeight, updateWeight, next );
            Check( hipGetLastError(), "UpdateKernel launch" );
            result.timing.kernelLaunches += 3;
            Check(
                hipEventRecord( device.kernelEnd, device.stream ),
                "kernel end event" );
            Check(
                hipEventSynchronize( device.kernelEnd ),
                "kernel event synchronize" );
            ++ result.timing.deviceSynchronizations;
            float kernelMilliseconds = 0.0F;
            Check(
                hipEventElapsedTime(
                    &kernelMilliseconds,
                    device.kernelBegin, device.kernelEnd ),
                "kernel elapsed time" );
            result.timing.kernelMilliseconds += kernelMilliseconds;
        }
        std::swap( device.base, device.work );
    }
    result.timing.advanceMilliseconds =
        Milliseconds( advanceBegin, Clock::now() );

    const auto downloadBegin = Clock::now();
    result.state.resize( cellValues );
    Check(
        hipMemcpy(
            result.state.data(), device.base, cellBytes,
            hipMemcpyDeviceToHost ),
        "final D2H" );
    ++ result.timing.deviceSynchronizations;
    result.timing.downloadMilliseconds =
        Milliseconds( downloadBegin, Clock::now() );
    result.timing.computeMilliseconds =
        result.timing.advanceMilliseconds - result.timing.mpiMilliseconds;
    return result;
}

} // namespace oneflow_1d
