#include "OneDEuler.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace oneflow_1d
{
namespace
{

int CellIndex( int component, int cell, int nx )
{
    return component * nx + cell;
}

int FaceIndex( int component, int face, int nx )
{
    return component * ( nx + 1 ) + face;
}

void CheckInput(
    const double * state, int nx, double gamma, double dt, double dx )
{
    if ( state == nullptr )
    {
        throw std::invalid_argument( "OneD Euler received a null state" );
    }
    if ( nx <= 0 || gamma <= 1.0 || dt <= 0.0 || dx <= 0.0 )
    {
        throw std::invalid_argument( "OneD Euler received invalid parameters" );
    }
}

void PhysicalFlux(
    const double * q, int stride, int index, double gamma, double flux[ 3 ],
    double & waveSpeed )
{
    const double rho = q[ 0 * stride + index ];
    const double momentum = q[ 1 * stride + index ];
    const double energy = q[ 2 * stride + index ];
    const double velocity = momentum / rho;
    const double pressure = ( gamma - 1.0 )
        * ( energy - 0.5 * momentum * momentum / rho );
    const double soundSpeed = std::sqrt( gamma * pressure / rho );

    flux[ 0 ] = momentum;
    flux[ 1 ] = momentum * velocity + pressure;
    flux[ 2 ] = ( energy + pressure ) * velocity;
    waveSpeed = std::abs( velocity ) + soundSpeed;
}

void CalculateRhs(
    const double * state, int nx, double gamma, double dx,
    EulerBoundary boundary, double * left, double * right, double * flux,
    double * residual )
{
    for ( int face = 0; face <= nx; ++ face )
    {
        int leftCell = face - 1;
        int rightCell = face;
        if ( boundary == EulerBoundary::Periodic )
        {
            leftCell = ( face + nx - 1 ) % nx;
            rightCell = face % nx;
        }
        else
        {
            leftCell = std::max( 0, leftCell );
            rightCell = std::min( nx - 1, rightCell );
        }

        double fluxLeft[ 3 ];
        double fluxRight[ 3 ];
        double speedLeft = 0.0;
        double speedRight = 0.0;
        PhysicalFlux( state, nx, leftCell, gamma, fluxLeft, speedLeft );
        PhysicalFlux( state, nx, rightCell, gamma, fluxRight, speedRight );
        const double speed = std::max( speedLeft, speedRight );

        for ( int component = 0; component < EulerComponents; ++ component )
        {
            const double qLeft = state[ CellIndex( component, leftCell, nx ) ];
            const double qRight = state[ CellIndex( component, rightCell, nx ) ];
            const int faceIndex = FaceIndex( component, face, nx );
            left[ faceIndex ] = qLeft;
            right[ faceIndex ] = qRight;
            flux[ faceIndex ] = 0.5 * ( fluxLeft[ component ] + fluxRight[ component ] )
                - 0.5 * speed * ( qRight - qLeft );
        }
    }

    const double inverseDx = 1.0 / dx;
    for ( int component = 0; component < EulerComponents; ++ component )
    {
        for ( int cell = 0; cell < nx; ++ cell )
        {
            residual[ CellIndex( component, cell, nx ) ] = -(
                flux[ FaceIndex( component, cell + 1, nx ) ]
                - flux[ FaceIndex( component, cell, nx ) ] ) * inverseDx;
        }
    }
}

} // namespace

void ResizeEulerTrace( EulerTrace & trace, int nx )
{
    trace.nx = nx;
    const int faceValues = EulerComponents * ( nx + 1 );
    const int cellValues = EulerComponents * nx;
    trace.faceLeft.resize( EulerRkStages * faceValues );
    trace.faceRight.resize( EulerRkStages * faceValues );
    trace.numericalFlux.resize( EulerRkStages * faceValues );
    trace.residual.resize( EulerRkStages * cellValues );
    trace.state.resize( ( EulerRkStages + 1 ) * cellValues );
}

void OneDCpuEulerStep(
    const double * input, int nx, double gamma, double dt, double dx,
    EulerBoundary boundary, EulerTrace & trace )
{
    CheckInput( input, nx, gamma, dt, dx );
    ResizeEulerTrace( trace, nx );

    const int faceValues = EulerComponents * ( nx + 1 );
    const int cellValues = EulerComponents * nx;
    std::copy( input, input + cellValues, trace.state.begin() );
    std::vector<double> current( input, input + cellValues );
    std::vector<double> next( cellValues );

    for ( int stage = 0; stage < EulerRkStages; ++ stage )
    {
        double * left = trace.faceLeft.data() + stage * faceValues;
        double * right = trace.faceRight.data() + stage * faceValues;
        double * flux = trace.numericalFlux.data() + stage * faceValues;
        double * residual = trace.residual.data() + stage * cellValues;
        CalculateRhs( current.data(), nx, gamma, dx, boundary,
                      left, right, flux, residual );

        double baseWeight = 0.0;
        double updateWeight = 1.0;
        if ( stage == 1 )
        {
            baseWeight = 0.75;
            updateWeight = 0.25;
        }
        else if ( stage == 2 )
        {
            baseWeight = 1.0 / 3.0;
            updateWeight = 2.0 / 3.0;
        }

        for ( int i = 0; i < cellValues; ++ i )
        {
            next[ i ] = baseWeight * input[ i ]
                + updateWeight * ( current[ i ] + dt * residual[ i ] );
        }
        std::copy( next.begin(), next.end(),
                   trace.state.begin() + ( stage + 1 ) * cellValues );
        current.swap( next );
    }
}

} // namespace oneflow_1d
