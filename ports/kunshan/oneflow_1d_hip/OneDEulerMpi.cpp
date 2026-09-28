#include "OneDEulerMpi.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <stdexcept>
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

int CellIndex( int component, int cell, int nx )
{
    return component * nx + cell;
}

int FaceIndex( int component, int face, int nx )
{
    return component * ( nx + 1 ) + face;
}

void CheckProblem(
    const EulerMpiProblem & problem,
    const std::vector< double > & initialState )
{
    if ( problem.communicator == MPI_COMM_NULL || problem.globalNx <= 0
         || problem.localNx <= 0 || problem.globalOffset < 0
         || problem.globalOffset + problem.localNx > problem.globalNx
         || problem.steps <= 0 || problem.gamma <= 1.0
         || problem.dt <= 0.0 || problem.dx <= 0.0 )
    {
        throw std::invalid_argument( "invalid MPI Euler problem" );
    }
    if ( initialState.size()
         != static_cast< std::size_t >( 3 * problem.localNx ) )
    {
        throw std::invalid_argument( "invalid MPI Euler initial state size" );
    }
}

void PhysicalFlux(
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
    waveSpeed = std::abs( velocity ) + std::sqrt( gamma * pressure / rho );
    return;
}

void ExchangeHalos(
    const double * current, int nx, int rank, int ranks, MPI_Comm communicator,
    double leftHalo[ 3 ], double rightHalo[ 3 ] )
{
    double sendLeft[ 3 ];
    double sendRight[ 3 ];
    for ( int component = 0; component < 3; ++ component )
    {
        sendLeft[ component ] = current[ CellIndex( component, 0, nx ) ];
        sendRight[ component ] =
            current[ CellIndex( component, nx - 1, nx ) ];
    }

    const int leftRank = ( rank + ranks - 1 ) % ranks;
    const int rightRank = ( rank + 1 ) % ranks;
    MPI_Sendrecv(
        sendLeft, 3, MPI_DOUBLE, leftRank, 100,
        rightHalo, 3, MPI_DOUBLE, rightRank, 100,
        communicator, MPI_STATUS_IGNORE );
    MPI_Sendrecv(
        sendRight, 3, MPI_DOUBLE, rightRank, 101,
        leftHalo, 3, MPI_DOUBLE, leftRank, 101,
        communicator, MPI_STATUS_IGNORE );
    return;
}

void CalculateFlux(
    const double * current, int nx, double gamma, const double leftHalo[ 3 ],
    const double rightHalo[ 3 ], double * flux )
{
    for ( int face = 0; face <= nx; ++ face )
    {
        double leftState[ 3 ];
        double rightState[ 3 ];
        for ( int component = 0; component < 3; ++ component )
        {
            leftState[ component ] = face == 0
                ? leftHalo[ component ]
                : current[ CellIndex( component, face - 1, nx ) ];
            rightState[ component ] = face == nx
                ? rightHalo[ component ]
                : current[ CellIndex( component, face, nx ) ];
        }

        double leftFlux[ 3 ];
        double rightFlux[ 3 ];
        double leftSpeed = 0.0;
        double rightSpeed = 0.0;
        PhysicalFlux( leftState, gamma, leftFlux, leftSpeed );
        PhysicalFlux( rightState, gamma, rightFlux, rightSpeed );
        const double speed = std::max( leftSpeed, rightSpeed );
        for ( int component = 0; component < 3; ++ component )
        {
            flux[ FaceIndex( component, face, nx ) ] =
                0.5 * ( leftFlux[ component ] + rightFlux[ component ] )
                - 0.5 * speed
                    * ( rightState[ component ] - leftState[ component ] );
        }
    }
    return;
}

void CalculateResidual(
    const double * flux, int nx, double inverseDx, double * residual )
{
    for ( int component = 0; component < 3; ++ component )
    {
        for ( int cell = 0; cell < nx; ++ cell )
        {
            residual[ CellIndex( component, cell, nx ) ] = -(
                flux[ FaceIndex( component, cell + 1, nx ) ]
                - flux[ FaceIndex( component, cell, nx ) ] ) * inverseDx;
        }
    }
    return;
}

void Update(
    const double * base, const double * current, const double * residual,
    int count, double dt, double baseWeight, double updateWeight, double * next )
{
    for ( int index = 0; index < count; ++ index )
    {
        next[ index ] = baseWeight * base[ index ]
            + updateWeight * ( current[ index ] + dt * residual[ index ] );
    }
    return;
}

std::uint64_t Mix( std::uint64_t value )
{
    value += UINT64_C( 0x9e3779b97f4a7c15 );
    value = ( value ^ ( value >> 30 ) ) * UINT64_C( 0xbf58476d1ce4e5b9 );
    value = ( value ^ ( value >> 27 ) ) * UINT64_C( 0x94d049bb133111eb );
    return value ^ ( value >> 31 );
}

} // namespace

EulerMpiResult RunCpuEulerMpi(
    const EulerMpiProblem & problem,
    const std::vector< double > & initialState )
{
    CheckProblem( problem, initialState );
    EulerMpiResult result;
    const int nx = problem.localNx;
    const int cellValues = 3 * nx;
    const int faceValues = 3 * ( nx + 1 );

    const auto createBegin = Clock::now();
    std::vector< double > base( cellValues );
    std::vector< double > work( cellValues );
    std::vector< double > scratch( cellValues );
    std::vector< double > flux( faceValues );
    std::vector< double > residual( cellValues );
    result.timing.createMilliseconds =
        Milliseconds( createBegin, Clock::now() );

    const auto uploadBegin = Clock::now();
    std::copy( initialState.begin(), initialState.end(), base.begin() );
    result.timing.uploadMilliseconds =
        Milliseconds( uploadBegin, Clock::now() );

    int rank = 0;
    int ranks = 0;
    MPI_Comm_rank( problem.communicator, &rank );
    MPI_Comm_size( problem.communicator, &ranks );
    const double inverseDx = 1.0 / problem.dx;

    const auto advanceBegin = Clock::now();
    for ( int step = 0; step < problem.steps; ++ step )
    {
        for ( int stage = 0; stage < 3; ++ stage )
        {
            const double * current = stage == 0 ? base.data()
                : ( stage == 1 ? work.data() : scratch.data() );
            double * next = stage == 0 ? work.data()
                : ( stage == 1 ? scratch.data() : work.data() );

            double leftHalo[ 3 ];
            double rightHalo[ 3 ];
            const auto mpiBegin = Clock::now();
            ExchangeHalos(
                current, nx, rank, ranks, problem.communicator,
                leftHalo, rightHalo );
            result.timing.mpiMilliseconds +=
                Milliseconds( mpiBegin, Clock::now() );
            ++ result.timing.haloExchanges;

            const auto computeBegin = Clock::now();
            CalculateFlux(
                current, nx, problem.gamma, leftHalo, rightHalo,
                flux.data() );
            CalculateResidual(
                flux.data(), nx, inverseDx, residual.data() );
            const double baseWeight = stage == 1 ? 0.75
                : ( stage == 2 ? 1.0 / 3.0 : 0.0 );
            const double updateWeight = stage == 0 ? 1.0
                : ( stage == 1 ? 0.25 : 2.0 / 3.0 );
            Update(
                base.data(), current, residual.data(), cellValues,
                problem.dt, baseWeight, updateWeight, next );
            result.timing.computeMilliseconds +=
                Milliseconds( computeBegin, Clock::now() );
        }
        base.swap( work );
    }
    result.timing.advanceMilliseconds =
        Milliseconds( advanceBegin, Clock::now() );

    const auto downloadBegin = Clock::now();
    result.state = base;
    result.timing.downloadMilliseconds =
        Milliseconds( downloadBegin, Clock::now() );
    return result;
}

std::uint64_t EulerMpiStateHash(
    const std::vector< double > & localState, int globalNx, int localNx,
    int globalOffset )
{
    if ( localState.size() != static_cast< std::size_t >( 3 * localNx ) )
    {
        throw std::invalid_argument( "invalid local state for Euler MPI hash" );
    }
    std::uint64_t hash = 0;
    for ( int component = 0; component < 3; ++ component )
    {
        for ( int cell = 0; cell < localNx; ++ cell )
        {
            std::uint64_t bits = 0;
            const double value = localState[ CellIndex( component, cell, localNx ) ];
            std::memcpy( &bits, &value, sizeof( bits ) );
            const std::uint64_t globalIndex = static_cast< std::uint64_t >(
                component * globalNx + globalOffset + cell );
            hash ^= Mix( bits ^ Mix( globalIndex ) );
        }
    }
    return hash;
}

} // namespace oneflow_1d
