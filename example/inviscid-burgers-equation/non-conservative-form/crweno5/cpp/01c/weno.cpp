#include "weno.h"
#include "hxmath.h"
#include "post.h"
#include <cmath>
#include <vector>
#include <print>
#include <numbers>
#include <fstream>

double wcL( double v1, double v2, double v3, double v4, double v5 )
{
    double eps = 1.0e-6;

    // smoothness indicators
    double s1 = ( 13.0 / 12.0 ) * SQR( v1 - 2.0 * v2 + v3 ) + 0.25 * SQR( v1 - 4.0 * v2 + 3.0 * v3 );
    double s2 = ( 13.0 / 12.0 ) * SQR( v2 - 2.0 * v3 + v4 ) + 0.25 * SQR( v2 - v4 );
    double s3 = ( 13.0 / 12.0 ) * SQR( v3 - 2.0 * v4 + v5 ) + 0.25 * SQR( 3.0 * v3 - 4.0 * v4 + v5 );

    // computing nonlinear weights w1, w2, w3
    double c1 = 1.0e-1 / ( SQR( eps + s1 ) );
    double c2 = 6.0e-1 / ( SQR( eps + s2 ) );
    double c3 = 3.0e-1 / ( SQR( eps + s3 ) );

    double w1 = c1 / ( c1 + c2 + c3 );
    double w2 = c2 / ( c1 + c2 + c3 );
    double w3 = c3 / ( c1 + c2 + c3 );

    // candiate stencils
    double q1 = v1 / 3.0 - 7.0 / 6.0 * v2 + 11.0 / 6.0 * v3;
    double q2 = -v2 / 6.0 + 5.0 / 6.0 * v3 + v4 / 3.0;
    double q3 = v3 / 3.0 + 5.0 / 6.0 * v4 - v5 / 6.0;

    // reconstructed value at interface
    double f = ( w1 * q1 + w2 * q2 + w3 * q3 );

    return f;
}

double wcR( double v1, double v2, double v3, double v4, double v5 )
{
    double eps = 1.0e-6;

    // smoothness indicators
    double s1 = ( 13.0 / 12.0 ) * SQR( v1 - 2.0 * v2 + v3 ) + 0.25 * SQR( v1 - 4.0 * v2 + 3.0 * v3 );
    double s2 = ( 13.0 / 12.0 ) * SQR( v2 - 2.0 * v3 + v4 ) + 0.25 * SQR( v2 - v4 );
    double s3 = ( 13.0 / 12.0 ) * SQR( v3 - 2.0 * v4 + v5 ) + 0.25 * SQR( 3.0 * v3 - 4.0 * v4 + v5 );

    // computing nonlinear weights w1, w2, w3
    double c1 = 3.0e-1 / SQR( eps + s1 );
    double c2 = 6.0e-1 / SQR( eps + s2 );
    double c3 = 1.0e-1 / SQR( eps + s3 );

    double w1 = c1 / ( c1 + c2 + c3 );
    double w2 = c2 / ( c1 + c2 + c3 );
    double w3 = c3 / ( c1 + c2 + c3 );

    // candiate stencils;
    double q1 = -v1 / 6.0 + 5.0 / 6.0 * v2 + v3 / 3.0;
    double q2 = v2 / 3.0 + 5.0 / 6.0 * v3 - v4 / 6.0;
    double q3 = 11.0 / 6.0 * v3 - 7.0 / 6.0 * v4 + v5 / 3.0;

    // reconstructed value at interface
    double f = ( w1 * q1 + w2 * q2 + w3 * q3 );

    return f;
}

//-----------------------------------------------------------------------------
// WENO reconstruction for upwind direction (positive; left to right)
// u(i): solution values at finite difference grid nodes i = 1,...,N+1
// f(j): reconstructed values at nodes j = i+1/2; j = 1,...,N
//-----------------------------------------------------------------------------
void wenoL( int nx, Vec1d & u, std::vector<double> & f )
{
    int i;
    double v1, v2, v3, v4, v5;

    i = 0;
    v1 = 3.0 * u[ i ] - 2.0 * u[ i + 1 ];
    v2 = 2.0 * u[ i ] - u[ i + 1 ];
    v3 = u[ i ];
    v4 = u[ i + 1 ];
    v5 = u[ i + 2 ];
    f[ i ] = wcL( v1, v2, v3, v4, v5 );

    i = 1;
    v1 = 2.0 * u[ i - 1 ] - u[ i ];
    v2 = u[ i - 1 ];
    v3 = u[ i ];
    v4 = u[ i + 1 ];
    v5 = u[ i + 2 ];
    f[ i ] = wcL( v1, v2, v3, v4, v5 );

    for ( int i = 2; i < nx - 1; ++ i )
    {
        v1 = u[ i - 2 ];
        v2 = u[ i - 1 ];
        v3 = u[ i ];
        v4 = u[ i + 1 ];
        v5 = u[ i + 2 ];
        f[ i ] = wcL( v1, v2, v3, v4, v5 );
    }

    i = nx - 1;
    v1 = u[ i - 2 ];
    v2 = u[ i - 1 ];
    v3 = u[ i ];
    v4 = u[ i + 1 ];
    v5 = 2.0 * u[ i + 1 ] - u[ i ];
    f[ i ] = wcL( v1, v2, v3, v4, v5 );
}

//-----------------------------------------------------------------------------
// CRWENO reconstruction for downwind direction (negative; right to left)
// u(i): solution values at finite difference grid nodes i = 1,...,N+1
// f(j): reconstructed values at nodes j = i-1/2; j = 2,...,N+1
//-----------------------------------------------------------------------------
void wenoR( int nx, Vec1d & u, std::vector<double> & f )
{
    int i;
    double v1, v2, v3, v4, v5;

    i = 1;
    v1 = 2.0 * u[ i - 1 ] - u[ i ];
    v2 = u[ i - 1 ];
    v3 = u[ i ];
    v4 = u[ i + 1 ];
    v5 = u[ i + 2 ];
    f[ i - 1 ] = wcR( v1, v2, v3, v4, v5 );
    for ( int i = 2; i < nx - 1; ++ i )
    {
        v1 = u[ i - 2 ];
        v2 = u[ i - 1 ];
        v3 = u[ i ];
        v4 = u[ i + 1 ];
        v5 = u[ i + 2 ];
        f[ i - 1 ] = wcR( v1, v2, v3, v4, v5 );
    }

    i = nx - 1;
    v1 = u[ i - 2 ];
    v2 = u[ i - 1 ];
    v3 = u[ i ];
    v4 = u[ i + 1 ];
    v5 = 2.0 * u[ i + 1 ] - u[ i ];
    f[ i - 1 ] = wcR( v1, v2, v3, v4, v5 );

    i = nx;
    v1 = u[ i - 2 ];
    v2 = u[ i - 1 ];
    v3 = u[ i ];
    v4 = 2.0 * u[ i ] - u[ i - 1 ];
    v5 = 3.0 * u[ i ] - 2.0 * u[ i - 1 ];
    f[ i - 1 ] = wcR( v1, v2, v3, v4, v5 );
}

//-----------------------------------------------------------------------------
// Calculate right hand term of the inviscid Burgers equation
//-----------------------------------------------------------------------------
void rhs( int nx, double dx, Vec1d & u, Vec1d & r )
{
    std::vector<double> uL( nx, 0 );
    std::vector<double> uR( nx, 0 );

    wenoL( nx, u, uL );
    wenoR( nx, u, uR );

    for ( int i = 1; i < nx; ++ i )
    {
        if ( u[ i ] >= 0.0 )
        {
            r[ i ] = - u[ i ] * ( uL[ i ] - uL[ i - 1 ] ) / dx;
        }
        else
        {
            r[ i ] = - u[ i ] * ( uR[ i ] - uR[ i - 1 ] ) / dx;
        }
    }
}

