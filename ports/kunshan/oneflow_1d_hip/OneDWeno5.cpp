#include "OneDWeno5.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace oneflow_1d
{
namespace
{

int PaddedIndex( int component, int cell, int nx )
{
    return component * ( nx + 2 * WenoGhostCells )
        + cell + WenoGhostCells;
}

int CellIndex( int component, int cell, int nx )
{
    return component * nx + cell;
}

int FaceIndex( int component, int face, int nx )
{
    return component * ( nx + 1 ) + face;
}

double Square( double value )
{
    return value * value;
}

double WenoLeft(
    double v1, double v2, double v3, double v4, double v5 )
{
    constexpr double epsilon = 1.0e-6;
    const double s1 = 13.0 / 12.0 * Square( v1 - 2.0 * v2 + v3 )
        + 0.25 * Square( v1 - 4.0 * v2 + 3.0 * v3 );
    const double s2 = 13.0 / 12.0 * Square( v2 - 2.0 * v3 + v4 )
        + 0.25 * Square( v2 - v4 );
    const double s3 = 13.0 / 12.0 * Square( v3 - 2.0 * v4 + v5 )
        + 0.25 * Square( 3.0 * v3 - 4.0 * v4 + v5 );
    const double c1 = 0.1 / Square( epsilon + s1 );
    const double c2 = 0.6 / Square( epsilon + s2 );
    const double c3 = 0.3 / Square( epsilon + s3 );
    const double sum = c1 + c2 + c3;
    const double q1 = v1 / 3.0 - 7.0 / 6.0 * v2 + 11.0 / 6.0 * v3;
    const double q2 = -v2 / 6.0 + 5.0 / 6.0 * v3 + v4 / 3.0;
    const double q3 = v3 / 3.0 + 5.0 / 6.0 * v4 - v5 / 6.0;
    return c1 / sum * q1 + c2 / sum * q2 + c3 / sum * q3;
}

double WenoRight(
    double v1, double v2, double v3, double v4, double v5 )
{
    constexpr double epsilon = 1.0e-6;
    const double s1 = 13.0 / 12.0 * Square( v1 - 2.0 * v2 + v3 )
        + 0.25 * Square( v1 - 4.0 * v2 + 3.0 * v3 );
    const double s2 = 13.0 / 12.0 * Square( v2 - 2.0 * v3 + v4 )
        + 0.25 * Square( v2 - v4 );
    const double s3 = 13.0 / 12.0 * Square( v3 - 2.0 * v4 + v5 )
        + 0.25 * Square( 3.0 * v3 - 4.0 * v4 + v5 );
    const double c1 = 0.3 / Square( epsilon + s1 );
    const double c2 = 0.6 / Square( epsilon + s2 );
    const double c3 = 0.1 / Square( epsilon + s3 );
    const double sum = c1 + c2 + c3;
    const double q1 = -v1 / 6.0 + 5.0 / 6.0 * v2 + v3 / 3.0;
    const double q2 = v2 / 3.0 + 5.0 / 6.0 * v3 - v4 / 6.0;
    const double q3 = 11.0 / 6.0 * v3 - 7.0 / 6.0 * v4 + v5 / 3.0;
    return c1 / sum * q1 + c2 / sum * q2 + c3 / sum * q3;
}

void FillPadded(
    const double * state, int nx, EulerBoundary boundary,
    std::vector<double> & padded )
{
    const int stride = nx + 2 * WenoGhostCells;
    padded.resize( EulerComponents * stride );
    for ( int component = 0; component < EulerComponents; ++ component )
    {
        for ( int cell = 0; cell < nx; ++ cell )
        {
            padded[ PaddedIndex( component, cell, nx ) ] =
                state[ CellIndex( component, cell, nx ) ];
        }
        for ( int ghost = 1; ghost <= WenoGhostCells; ++ ghost )
        {
            const int leftSource = boundary == EulerBoundary::Periodic
                ? nx - ghost : 0;
            const int rightSource = boundary == EulerBoundary::Periodic
                ? ghost - 1 : nx - 1;
            padded[ PaddedIndex( component, -ghost, nx ) ] =
                state[ CellIndex( component, leftSource, nx ) ];
            padded[ PaddedIndex( component, nx - 1 + ghost, nx ) ] =
                state[ CellIndex( component, rightSource, nx ) ];
        }
    }
}

void CalculateRhs(
    const double * state, int nx, double gamma, double dx,
    EulerBoundary boundary, double * splitPositive, double * splitNegative,
    double * reconstructedPositive, double * reconstructedNegative,
    double * numericalFlux, double * residual )
{
    std::vector<double> padded;
    FillPadded( state, nx, boundary, padded );
    for ( int cell = -WenoGhostCells;
          cell < nx + WenoGhostCells; ++ cell )
    {
        const double rho = padded[ PaddedIndex( 0, cell, nx ) ];
        const double momentum = padded[ PaddedIndex( 1, cell, nx ) ];
        const double energy = padded[ PaddedIndex( 2, cell, nx ) ];
        const double velocity = momentum / rho;
        const double specificEnergy = energy / rho;
        const double pressure = ( gamma - 1.0 )
            * ( rho * specificEnergy - 0.5 * rho * velocity * velocity );
        const double enthalpy = specificEnergy + pressure / rho;
        const double weightedVelocity = std::sqrt( std::abs( rho ) )
            * velocity;
        const double weightedEnthalpy = std::sqrt( std::abs( rho ) )
            * enthalpy;
        const double soundSpeed = std::sqrt( std::abs(
            ( gamma - 1.0 )
            * ( weightedEnthalpy
                - 0.5 * weightedVelocity * weightedVelocity ) ) );
        const double speed = std::abs( soundSpeed + weightedVelocity );
        const double physical[ 3 ]{
            momentum,
            momentum * velocity + pressure,
            momentum * energy / rho + pressure * momentum / rho
        };
        for ( int component = 0; component < EulerComponents; ++ component )
        {
            const int index = PaddedIndex( component, cell, nx );
            splitPositive[ index ] = 0.5
                * ( physical[ component ] + speed * padded[ index ] );
            splitNegative[ index ] = 0.5
                * ( physical[ component ] - speed * padded[ index ] );
        }
    }

    for ( int component = 0; component < EulerComponents; ++ component )
    {
        for ( int face = 0; face <= nx; ++ face )
        {
            const int center = face - 1;
            const int faceIndex = FaceIndex( component, face, nx );
            reconstructedPositive[ faceIndex ] = WenoLeft(
                splitPositive[ PaddedIndex( component, center - 2, nx ) ],
                splitPositive[ PaddedIndex( component, center - 1, nx ) ],
                splitPositive[ PaddedIndex( component, center, nx ) ],
                splitPositive[ PaddedIndex( component, center + 1, nx ) ],
                splitPositive[ PaddedIndex( component, center + 2, nx ) ] );
            reconstructedNegative[ faceIndex ] = WenoRight(
                splitNegative[ PaddedIndex( component, center - 1, nx ) ],
                splitNegative[ PaddedIndex( component, center, nx ) ],
                splitNegative[ PaddedIndex( component, center + 1, nx ) ],
                splitNegative[ PaddedIndex( component, center + 2, nx ) ],
                splitNegative[ PaddedIndex( component, center + 3, nx ) ] );
            numericalFlux[ faceIndex ] =
                reconstructedPositive[ faceIndex ]
                + reconstructedNegative[ faceIndex ];
        }
        for ( int cell = 0; cell < nx; ++ cell )
        {
            residual[ CellIndex( component, cell, nx ) ] = -(
                ( reconstructedPositive[ FaceIndex( component, cell + 1, nx ) ]
                    - reconstructedPositive[ FaceIndex( component, cell, nx ) ] ) / dx
                + ( reconstructedNegative[ FaceIndex( component, cell + 1, nx ) ]
                    - reconstructedNegative[ FaceIndex( component, cell, nx ) ] ) / dx );
        }
    }
}

} // namespace

void ResizeWeno5Trace( Weno5Trace & trace, int nx )
{
    trace.nx = nx;
    const int paddedValues = EulerComponents
        * ( nx + 2 * WenoGhostCells );
    const int faceValues = EulerComponents * ( nx + 1 );
    const int cellValues = EulerComponents * nx;
    trace.splitPositive.resize( EulerRkStages * paddedValues );
    trace.splitNegative.resize( EulerRkStages * paddedValues );
    trace.reconstructedPositive.resize( EulerRkStages * faceValues );
    trace.reconstructedNegative.resize( EulerRkStages * faceValues );
    trace.numericalFlux.resize( EulerRkStages * faceValues );
    trace.residual.resize( EulerRkStages * cellValues );
    trace.state.resize( ( EulerRkStages + 1 ) * cellValues );
}

void OneDCpuLaxWeno5Step(
    const double * input, int nx, double gamma, double dt, double dx,
    EulerBoundary boundary, Weno5Trace & trace )
{
    if ( input == nullptr || nx < WenoGhostCells || gamma <= 1.0
        || dt <= 0.0 || dx <= 0.0 )
    {
        throw std::invalid_argument( "invalid CPU Lax-WENO5 input" );
    }
    ResizeWeno5Trace( trace, nx );
    const int paddedValues = EulerComponents
        * ( nx + 2 * WenoGhostCells );
    const int faceValues = EulerComponents * ( nx + 1 );
    const int cellValues = EulerComponents * nx;
    std::copy( input, input + cellValues, trace.state.begin() );
    std::vector<double> current( input, input + cellValues );
    std::vector<double> next( cellValues );

    for ( int stage = 0; stage < EulerRkStages; ++ stage )
    {
        CalculateRhs(
            current.data(), nx, gamma, dx, boundary,
            trace.splitPositive.data() + stage * paddedValues,
            trace.splitNegative.data() + stage * paddedValues,
            trace.reconstructedPositive.data() + stage * faceValues,
            trace.reconstructedNegative.data() + stage * faceValues,
            trace.numericalFlux.data() + stage * faceValues,
            trace.residual.data() + stage * cellValues );
        const double baseWeight = stage == 1 ? 0.75
            : ( stage == 2 ? 1.0 / 3.0 : 0.0 );
        const double updateWeight = stage == 0 ? 1.0
            : ( stage == 1 ? 0.25 : 2.0 / 3.0 );
        const double * residual =
            trace.residual.data() + stage * cellValues;
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
