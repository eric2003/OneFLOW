#include "crweno.h"
#include "hxmath.h"
#include "post.h"
#include <cmath>
#include <vector>
#include <print>
#include <numbers>
#include <fstream>

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

    a1 = ( 2.0 * w1 + w2 ) / 3.0;
    a2 = ( w1 + 2.0 * w2 + 2.0 * w3 ) / 3.0;
    a3 = w3 / 3.0;

    b1 = w1 / 6.0;
    b2 = ( 5.0 * w1 + 5.0 * w2 + w3 ) / 6.0;
    b3 = ( w2 + 5.0 * w3 ) / 6.0;
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

    a1 = w1 / 3.0;
    a2 = ( w3 + 2.0 * w2 + 2.0 * w1 ) / 3.0;
    a3 = ( 2.0 * w3 + w2 ) / 3.0;

    b1 = ( w2 + 5.0 * w1 ) / 6.0;
    b2 = ( 5.0 * w3 + 5.0 * w2 + w1 ) / 6.0;
    b3 = w3 / 6.0;
}

void crwenoL( int nx, std::vector<double> & u, std::vector<double> & f )
{
    std::vector<double> a( nx );
    std::vector<double> b( nx );
    std::vector<double> c( nx );
    std::vector<double> r( nx );

    int i;
    double v1, v2, v3, v4, v5;
    double a1, a2, a3, b1, b2, b3;

    i = 0;
    b[ i ] = 2.0 / 3.0;
    c[ i ] = 1.0 / 3.0;
    r[ i ] = ( u[ i ] + 5.0 * u[ i + 1 ] ) / 6.0;

    i = 1;
    v1 = 2.0 * u[ i - 1 ] - u[ i ];
    v2 = u[ i - 1 ];
    v3 = u[ i ];
    v4 = u[ i + 1 ];
    v5 = u[ i + 2 ];

    crwcL( v1, v2, v3, v4, v5, a1, a2, a3, b1, b2, b3 );
    a[ i ] = a1;
    b[ i ] = a2;
    c[ i ] = a3;
    r[ i ] = b1 * u[ i - 1 ] + b2 * u[ i ] + b3 * u[ i + 1 ];

    for ( int i = 2; i < nx - 1; ++ i )
    {
        v1 = u[ i - 2 ];
        v2 = u[ i - 1 ];
        v3 = u[ i ];
        v4 = u[ i + 1 ];
        v5 = u[ i + 2 ];
        crwcL( v1, v2, v3, v4, v5, a1, a2, a3, b1, b2, b3 );
        a[ i ] = a1;
        b[ i ] = a2;
        c[ i ] = a3;
        r[ i ] = b1 * u[ i - 1 ] + b2 * u[ i ] + b3 * u[ i + 1 ];
    }

    i = nx - 1;
    a[ i ] = 1.0 / 3.0;
    b[ i ] = 2.0 / 3.0;
    r[ i ] = ( 5.0 * u[ i ] + u[ i + 1 ] ) / 6.0;

    thomas_algorithm( a, b, c, r, f );
}

void crwenoR( int nx, std::vector<double> & u, std::vector<double> & f )
{
    std::vector<double> a( nx );
    std::vector<double> b( nx );
    std::vector<double> c( nx );
    std::vector<double> r( nx );

    int i;
    double v1, v2, v3, v4, v5;
    double a1, a2, a3, b1, b2, b3;

    i = 1;
    b[ i - 1 ] = 2.0 / 3.0;
    c[ i - 1 ] = 1.0 / 3.0;
    r[ i - 1 ] = ( u[ i - 1 ] + 5.0 * u[ i ] ) / 6.0;

    for ( int i = 2; i < nx - 1; ++ i )
    {
        v1 = u[ i - 2 ];
        v2 = u[ i - 1 ];
        v3 = u[ i ];
        v4 = u[ i + 1 ];
        v5 = u[ i + 2 ];

        crwcR( v1, v2, v3, v4, v5, a1, a2, a3, b1, b2, b3 );

        a[ i - 1 ] = a1;
        b[ i - 1 ] = a2;
        c[ i - 1 ] = a3;
        r[ i - 1 ] = b1 * u[ i - 1 ] + b2 * u[ i ] + b3 * u[ i + 1 ];
    }

    i = nx - 1;
    v1 = u[ i - 2 ];
    v2 = u[ i - 1 ];
    v3 = u[ i ];
    v4 = u[ i + 1 ];
    v5 = 2.0 * u[ i + 1 ] - u[ i ];

    crwcR( v1, v2, v3, v4, v5, a1, a2, a3, b1, b2, b3 );
    a[ i - 1 ] = a1;
    b[ i - 1 ] = a2;
    c[ i - 1 ] = a3;
    r[ i - 1 ] = b1 * u[ i - 1 ] + b2 * u[ i ] + b3 * u[ i + 1 ];

    i = nx;
    a[ i - 1 ] = 1.0 / 3.0;
    b[ i - 1 ] = 2.0 / 3.0;
    r[ i - 1 ] = ( 5.0 * u[ i - 1 ] + u[ i ] ) / 6.0;

    thomas_algorithm( a, b, c, r, f );
}

//-----------------------------------------------------------------------------
// Calculate right hand term of the inviscid Burgers equation
//-----------------------------------------------------------------------------
void rhs_crweno( int nx, double dx, std::vector<double> & u, std::vector<double> & r )
{
    std::vector<double> uL( nx, 0 );
    std::vector<double> uR( nx, 0 );

    crwenoL( nx, u, uL );
    crwenoR( nx, u, uR );

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

void mynumerical( RhsPtr rhs )
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

    std::vector<std::vector<double>> u( nx + 1 );
    for ( int i = 0; i < u.size(); ++ i )
    {
        u[ i ].resize( ns + 1 );
    }

    std::vector<double> x( nx + 1, 0 );
    std::vector<double> un( nx + 1, 0 ); // numerical solsution at every time step
    std::vector<double> ut( nx + 1, 0 ); // temporary array during RK3 integration
    std::vector<double> r( nx, 0 );

    int k = 0; // record index
    int freq = std::round( nt / ns );
    std::print( "freq = {}\n", freq );

    for ( int i = 0; i < nx + 1; ++ i )
    {
        x[ i ] = dx * ( i );
        un[ i ] = std::sin( 2.0 * std::numbers::pi * x[ i ] );
        u[ i ][ k ] = un[ i ]; //store solution at t = 0;
    }

    //dirichlet boundary condition
    u[ 0 ][ k ] = 0.0;
    u[ nx ][ k ] = 0.0;

    un[ 0 ] = 0.0;
    un[ nx ] = 0.0;

    //dirichlet boundary condition for temporary array
    ut[ 0 ] = 0.0;
    ut[ nx ] = 0.0;

    for ( int j = 1; j < nt + 1; ++ j )
    {
        rhs( nx, dx, un, r );
        for ( int i = 1; i < nx; ++ i )
        {
            ut[ i ] = un[ i ] + dt * r[ i ];
        }
        rhs( nx, dx, ut, r );

        for ( int i = 1; i < nx; ++ i )
        {
            ut[ i ] = 0.75 * un[ i ] + 0.25 * ut[ i ] + 0.25 * dt * r[ i ];
        }
        rhs( nx, dx, ut, r );

        for ( int i = 1; i < nx; ++ i )
        {
            un[ i ] = ( 1.0 / 3.0 ) * un[ i ] + ( 2.0 / 3.0 ) * ut[ i ] + ( 2.0 / 3.0 ) * dt * r[ i ];
        }

        if ( j % freq == 0 )
        {
            k = k + 1;
            for ( int i = 0; i < nx + 1; ++ i )
            {
                u[ i ][ k ] = un[ i ]; 
            }
            std::print( "k={}, ns={}\n", k, ns );
        }
    }

    DumpCsvFile( "field_final.csv", x, u );
}


