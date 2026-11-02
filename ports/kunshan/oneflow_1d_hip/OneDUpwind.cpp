#include "OneDUpwind.h"
#include <cmath>
#include <stdexcept>

void OneDCpuUpwindResidual(
    const double * u, int uSize, int uStart, int nx,
    double c, double dx, double * left, double * right, double * residual )
{
    if ( u == nullptr || left == nullptr || right == nullptr || residual == nullptr )
    {
        throw std::invalid_argument( "OneD CPU upwind received a null pointer" );
    }
    if ( nx <= 0 || dx <= 0.0 || uStart > -1 || uSize < nx - uStart + 1 )
    {
        throw std::invalid_argument( "OneD CPU upwind received an invalid extent" );
    }

    for ( int i = 0; i <= nx; ++ i )
    {
        const int ii = i - 1;
        left[ i ] = u[ ii - uStart ];
        right[ i ] = u[ ii + 1 - uStart ];
    }

    const double cp = 0.5 * ( c + std::abs( c ) );
    const double cn = 0.5 * ( c - std::abs( c ) );
    for ( int i = 0; i < nx; ++ i )
    {
        const double fip = cp * left[ i + 1 ] + cn * right[ i + 1 ];
        const double fim = cp * left[ i ] + cn * right[ i ];
        residual[ i ] = -( fip - fim ) / dx;
    }
}
