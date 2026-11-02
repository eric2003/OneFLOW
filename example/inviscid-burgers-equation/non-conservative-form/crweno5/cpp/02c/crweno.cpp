#include "crweno.h"
#include "hxmath.h"
#include "vec1d.h"
#include "post.h"
#include <cmath>
#include <vector>
#include <print>
#include <numbers>
#include <fstream>

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

void crwenoL( int ni, Vec1d & u, std::vector<double> & f )
{
    std::vector<double> a( ni + 1 );
    std::vector<double> b( ni + 1 );
    std::vector<double> c( ni + 1 );
    std::vector<double> r( ni + 1 );

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

    thomas_algorithm( a, b, c, r, f );
}

void crwenoR( int ni, Vec1d & u, std::vector<double> & f )
{
    std::vector<double> a( ni + 1 );
    std::vector<double> b( ni + 1 );
    std::vector<double> c( ni + 1 );
    std::vector<double> r( ni + 1 );

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

    thomas_algorithm( a, b, c, r, f );
}

//-----------------------------------------------------------------------------
// Calculate right hand term of the inviscid Burgers equation
//-----------------------------------------------------------------------------
void rhs_crweno( int nx, double dx, Vec1d & u, Vec1d & r )
{
    int ni = nx + 1;
    std::vector<double> uL( nx + 2, 0 );
    std::vector<double> uR( nx + 2, 0 );

    crwenoL( ni, u, uL );
    crwenoR( ni, u, uR );

    for ( int i = 0; i < ni; ++ i )
    {
        if ( u[ i ] >= 0.0 )
        {
            r[ i ] = - u[ i ] * ( uL[ i + 1 ] - uL[ i ] ) / dx;
        }
        else
        {
            r[ i ] = - u[ i ] * ( uR[ i + 1 ] - uR[ i ] ) / dx;
        }
    }
}

void boundary( int ni, Vec1d & u )
{
//         x=0                    ui+1/2,L                 x=L
// o  o  o  o  o  o ...  *   *   *  |   *    *       ...    o    o   o   o
//-3 -2 -1 i=0 1  2     i-2 i-1  i  |  i+1  i+2 i+3       i=N-1  N  N+1 N+2
//
//         x=0                    ui+1/2,R                 x=L
// o  o  o  o  o  o ...      *   *  |   *    *    *  ...    o    o   o   o
//-3 -2 -1 i=0 1  2     i-2 i-1  i  |  i+1  i+2  i+3      i=N-1  N  N+1 N+2

    //left bc
    int i = 0;
    u[ i ] = 0.0;
    u[ i - 1 ] = 2.0 * u[ i ] - u[ i + 1 ];
    u[ i - 2 ] = 3.0 * u[ i ] - 2.0 * u[ i + 1 ];
    u[ i - 3 ] = 4.0 * u[ i ] - 3.0 * u[ i + 1 ];

    //right bc
    i = ni - 1;
    u[ i ] = 0.0;
    u[ i + 1 ] = 2.0 * u[ i ] - u[ i - 1 ];
    u[ i + 2 ] = 3.0 * u[ i ] - 2.0 * u[ i - 1 ];
    u[ i + 3 ] = 4.0 * u[ i ] - 3.0 * u[ i - 1 ];
}

void numerical( RhsPtr rhs )
{
    int nx = 200;
    int ns = 10;
    double dt = 0.0001;
    double tm = 0.25;

    double dx = 1.0 / nx;
    int nt = std::round( tm / dt );
    double ds = tm / ns;

    std::print( "nx={}\n", nx );
    std::print( "ns={}\n", ns );
    std::print( "dt={}\n", dt );
    std::print( "tm={}\n", tm );
    std::print( "dx={}\n", dx );
    std::print( "nt={}\n", nt );

    int ni = nx + 1;

    int ighost = 3;
    int ist = 0 - ighost;
    int ied = ni - 1 + ighost;

    Vec1d un;
    un.Allocate( ist, ied, 0 ); // numerical solsution at every time step

    Vec1d ut;
    ut.Allocate( ist, ied, 0 ); // temporary array during RK3 integration

    Vec1d res;
    res.Allocate( 0, ni, 0 ); //N+1

    Vec1d x;
    x.Allocate( 0, ni - 1, 0 ); // npoints = N = nx + 1

    for ( int i = 0; i < ni; ++ i )
    {
        x[ i ] = dx * ( i );
        ut[ i ] = std::sin( 2.0 * std::numbers::pi * x[ i ] );
    }

    int k = 0; // record index
    int freq = std::round( nt / ns );
    std::print( "freq = {}\n", freq );

    boundary( ni, ut );
    un = ut;

    int iter = 0;
    DumpField( iter, x, un );

    for ( int it = 0; it < nt; ++ it )
    {
        rhs( nx, dx, un, res );
        for ( int i = 0; i < ni; ++ i )
        {
            ut[ i ] = un[ i ] + dt * res[ i ];
        }

        boundary( ni, ut );
        rhs( nx, dx, ut, res );

        for ( int i = 0; i < ni; ++ i )
        {
            ut[ i ] = 0.75 * un[ i ] + 0.25 * ut[ i ] + 0.25 * dt * res[ i ];
        }
        boundary( ni, ut );
        rhs( nx, dx, ut, res );

        for ( int i = 0; i < ni; ++ i )
        {
            un[ i ] = ( 1.0 / 3.0 ) * un[ i ] + ( 2.0 / 3.0 ) * ut[ i ] + ( 2.0 / 3.0 ) * dt * res[ i ];
        }

        if ( ( it + 1 ) % freq == 0 )
        {
            k = k + 1;
            std::print( "k={}, ns={}\n", k, ns );
            DumpField( it + 1, x, un );
        }
    }
}


