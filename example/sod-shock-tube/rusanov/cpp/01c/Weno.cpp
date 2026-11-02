#include "Weno.h"
#include "hxmath.h"
#include <format>
#include <numbers>
#include <print>
#include <fstream>
#include <iostream>

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
void wenoL( int N, Vec1d & u, Vec1d & f )
{
    for ( int i = 0; i <= N; ++ i )
    {
        int ii = i - 1;
        double v1 = u[ ii - 2 ];
        double v2 = u[ ii - 1 ];
        double v3 = u[ ii ];
        double v4 = u[ ii + 1 ];
        double v5 = u[ ii + 2 ];
        f[ i ] = wcL( v1, v2, v3, v4, v5 );
    }
}

//-----------------------------------------------------------------------------
// CRWENO reconstruction for downwind direction (negative; right to left)
// u(i): solution values at finite difference grid nodes i = 1,...,N+1
// f(j): reconstructed values at nodes j = i-1/2; j = 2,...,N+1
//-----------------------------------------------------------------------------
void wenoR( int N, Vec1d & u, Vec1d & f )
{
    for ( int i = 0; i <= N; ++ i )
    {
        int ii = i - 1;
        double v1 = u[ ii - 1 ];
        double v2 = u[ ii ];
        double v3 = u[ ii + 1 ];
        double v4 = u[ ii + 2 ];
        double v5 = u[ ii + 3 ];
        f[ i ] = wcR( v1, v2, v3, v4, v5 );
    }
}

void crabL( double w1, double w2, double w3,
    double & a1, double & a2, double & a3,
    double & b1, double & b2, double & b3 )
{
    a1 = ( 2.0 * w1 + w2 ) / 3.0;
    a2 = ( w1 + 2.0 * w2 + 2.0 * w3 ) / 3.0;
    a3 = w3 / 3.0;

    b1 = w1 / 6.0;
    b2 = ( 5.0 * w1 + 5.0 * w2 + w3 ) / 6.0;
    b3 = ( w2 + 5.0 * w3 ) / 6.0;
}

void crabR( double w1, double w2, double w3,
    double & a1, double & a2, double & a3,
    double & b1, double & b2, double & b3 )
{
    a1 = w1 / 3.0;
    a2 = ( w3 + 2.0 * w2 + 2.0 * w1 ) / 3.0;
    a3 = ( 2.0 * w3 + w2 ) / 3.0;

    b1 = ( w2 + 5.0 * w1 ) / 6.0;
    b2 = ( 5.0 * w3 + 5.0 * w2 + w1 ) / 6.0;
    b3 = w3 / 6.0;
}

void crwcL( double v1, double v2, double v3, double v4, double v5,
    double & a1, double & a2, double & a3,
    double & b1, double & b2, double & b3 )
{
    double eps = 1.0e-6;

    // smoothness indicators
    double s1 = ( 13.0 / 12.0 ) * SQR( v1 - 2.0 * v2 + v3 ) + 0.25 * SQR( v1 - 4.0 * v2 + 3.0 * v3 );
    double s2 = ( 13.0 / 12.0 ) * SQR( v2 - 2.0 * v3 + v4 ) + 0.25 * SQR( v2 - v4 );
    double s3 = ( 13.0 / 12.0 ) * SQR( v3 - 2.0 * v4 + v5 ) + 0.25 * SQR( 3.0 * v3 - 4.0 * v4 + v5 );

    // computing nonlinear weights w1, w2, w3
    double c1 = 2.0e-1 / ( SQR( eps + s1 ) );
    double c2 = 5.0e-1 / ( SQR( eps + s2 ) );
    double c3 = 3.0e-1 / ( SQR( eps + s3 ) );

    double w1 = c1 / ( c1 + c2 + c3 );
    double w2 = c2 / ( c1 + c2 + c3 );
    double w3 = c3 / ( c1 + c2 + c3 );

    crabL( w1, w2, w3, a1, a2, a3, b1, b2, b3 );
}

void crwcR( double v1, double v2, double v3, double v4, double v5,
    double & a1, double & a2, double & a3,
    double & b1, double & b2, double & b3 )
{
    double eps = 1.0e-6;

    // smoothness indicators
    double s1 = ( 13.0 / 12.0 ) * SQR( v1 - 2.0 * v2 + v3 ) + 0.25 * SQR( v1 - 4.0 * v2 + 3.0 * v3 );
    double s2 = ( 13.0 / 12.0 ) * SQR( v2 - 2.0 * v3 + v4 ) + 0.25 * SQR( v2 - v4 );
    double s3 = ( 13.0 / 12.0 ) * SQR( v3 - 2.0 * v4 + v5 ) + 0.25 * SQR( 3.0 * v3 - 4.0 * v4 + v5 );

    // computing nonlinear weights w1, w2, w3
    double c1 = 3.0e-1 / SQR( eps + s1 );
    double c2 = 5.0e-1 / SQR( eps + s2 );
    double c3 = 2.0e-1 / SQR( eps + s3 );

    double w1 = c1 / ( c1 + c2 + c3 );
    double w2 = c2 / ( c1 + c2 + c3 );
    double w3 = c3 / ( c1 + c2 + c3 );

    crabR( w1, w2, w3, a1, a2, a3, b1, b2, b3 );
}

void crwenoL( int ni, Vec1d & u, Vec1d & f )
{
    std::vector<double> a( ni + 1 );
    std::vector<double> b( ni + 1 );
    std::vector<double> c( ni + 1 );
    std::vector<double> r( ni + 1 );
    std::vector<double> y( ni + 1 );

    int i, ii;
    double v1, v2, v3, v4, v5;
    double a1, a2, a3, b1, b2, b3;

    int ist = -1;

    i = -1;
    ii = i - ist;
    crabL( 0, 0, 1, a1, a2, a3, b1, b2, b3 );
    a[ ii ] = a1;
    b[ ii ] = a2;
    c[ ii ] = a3;
    r[ ii ] = b1 * u[ i - 1 ] + b2 * u[ i ] + b3 * u[ i + 1 ];

    for ( int i = 0; i < ni - 1; ++ i )
    {
        ii = i - ist;
        v1 = u[ i - 2 ];
        v2 = u[ i - 1 ];
        v3 = u[ i ];
        v4 = u[ i + 1 ];
        v5 = u[ i + 2 ];
        crwcL( v1, v2, v3, v4, v5, a1, a2, a3, b1, b2, b3 );
        a[ ii ] = a1;
        b[ ii ] = a2;
        c[ ii ] = a3;
        r[ ii ] = b1 * u[ i - 1 ] + b2 * u[ i ] + b3 * u[ i + 1 ];
    }

    i = ni - 1;
    ii = i - ist;
    crabL( 1, 0, 0, a1, a2, a3, b1, b2, b3 );
    a[ ii ] = a1;
    b[ ii ] = a2;
    c[ ii ] = a3;
    r[ ii ] = b1 * u[ i - 1 ] + b2 * u[ i ] + b3 * u[ i + 1 ];

    thomas_algorithm( a, b, c, r, y );

    for ( int i = 0; i < ni + 1; ++ i )
    {
        f[ i ] = y[ i ];
    }

}

void crwenoR( int ni, Vec1d & u, Vec1d & f )
{
    std::vector<double> a( ni + 1 );
    std::vector<double> b( ni + 1 );
    std::vector<double> c( ni + 1 );
    std::vector<double> r( ni + 1 );
    std::vector<double> y( ni + 1 );

    int i, ii;
    double v1, v2, v3, v4, v5;
    double a1, a2, a3, b1, b2, b3;

    int ist = -1;

    i = -1;
    ii = i - ist;
    crabR( 0, 0, 1, a1, a2, a3, b1, b2, b3 );
    a[ ii ] = a1;
    b[ ii ] = a2;
    c[ ii ] = a3;
    r[ ii ] = b1 * u[ i ] + b2 * u[ i + 1 ] + b3 * u[ i + 2 ];

    for ( int i = 0; i < ni - 1; ++ i )
    {
        ii = i - ist;
        v1 = u[ i - 1 ];
        v2 = u[ i ];
        v3 = u[ i + 1 ];
        v4 = u[ i + 2 ];
        v5 = u[ i + 3 ];

        crwcR( v1, v2, v3, v4, v5, a1, a2, a3, b1, b2, b3 );

        a[ ii ] = a1;
        b[ ii ] = a2;
        c[ ii ] = a3;
        r[ ii ] = b1 * u[ i ] + b2 * u[ i + 1 ] + b3 * u[ i + 2 ];
    }

    i = ni - 1;
    ii = i - ist;
    crabR( 1, 0, 0, a1, a2, a3, b1, b2, b3 );
    a[ ii ] = a1;
    b[ ii ] = a2;
    c[ ii ] = a3;
    r[ ii ] = b1 * u[ i ] + b2 * u[ i + 1 ] + b3 * u[ i + 2 ];

    thomas_algorithm( a, b, c, r, y );

    for ( int i = 0; i < ni + 1; ++ i )
    {
        f[ i ] = y[ i ];
    }
}


void crwenoL( int ni, VecWrap & u, VecWrap & f )
{
    int nequ = u.get_nequ();
    for ( int m = 0; m < nequ; ++ m )
    {
        crwenoL( ni, u.vec( m ), f.vec( m ) );
    }
}

void crwenoR( int ni, VecWrap & u, VecWrap & f )
{
    int nequ = u.get_nequ();
    for ( int m = 0; m < nequ; ++ m )
    {
        crwenoR( ni, u.vec( m ), f.vec( m ) );
    }
}

void wenoL( int ni, VecWrap & u, VecWrap & f )
{
    int nequ = u.get_nequ();
    for ( int m = 0; m < nequ; ++ m )
    {
        wenoL( ni, u.vec( m ), f.vec( m ) );
    }
}

void wenoR( int ni, VecWrap & u, VecWrap & f )
{
    int nequ = u.get_nequ();
    for ( int m = 0; m < nequ; ++ m )
    {
        wenoR( ni, u.vec( m ), f.vec( m ) );
    }
}