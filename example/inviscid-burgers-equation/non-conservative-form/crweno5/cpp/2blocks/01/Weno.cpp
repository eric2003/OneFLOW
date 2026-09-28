#include "Weno.h"
#include "hxmath.h"
#include "Global.h"
#include "Grid.h"
#include "cgnslib.h"
#include <format>
#include <numbers>
#include <print>
#include <fstream>
#include <iostream>

void WenoField::Init( Grid * grid )
{
    this->ni = grid->x.size();
    Vec1d & x = grid->x;
    this->dx = std::abs( x[ 1 ] - x[ 0 ] );
    this->dt = 0.0001;
    this->total_time = 0.25;
    this->nt = std::round( total_time / dt );
    this->ns = 10;

    std::print( "ni={}\n",ni );
    std::print( "ns={}\n", ns );
    std::print( "dt={}\n", dt );
    std::print( "total_time={}\n", total_time );
    std::print( "dx={}\n", dx );
    std::print( "nt={}\n", nt );
    std::cout << "this->dt   = " << this->dt << "\n";
    std::cout << "this->total_time  = " << this->total_time << "\n";
    std::cout << "this->nt   = " << this->nt << "\n";
    std::cout << "this->ni   = " << this->ni << "\n";
    std::cout << "nt * dt = " << nt * dt << "\n";

    int freq = std::round( nt / ns );
    std::print( "freq = {}\n", freq );
    Global::nt = nt;

    int ighost = 3;
    int ist = 0 - ighost;
    int ied = this->ni - 1 + ighost;

    std::vector<Vec1d> ulist( ns + 1 );

    for ( int i = 0; i < ulist.size(); ++ i )
    {
        ulist[ i ].Allocate( ist, ied, 0 );
    }

    un.Allocate( ist, ied, 0 ); // numerical solsution at every time step
    u.Allocate( ist, ied, 0 ); // temporary array during RK3 integration
    res.Allocate( 0, this->ni, 0 ); //N+1

    for ( int i = 0; i < ni; ++ i )
    {
        u[ i ] = std::sin( 2.0 * std::numbers::pi * x[ i ] );
    }

    int kkk = 1;
}

void WenoField::PhysicalBoundary( Zone * zone )
{
    int nbccos = zone->bccos.size();
    for ( int ibcco = 0; ibcco < nbccos; ++ ibcco )
    {
        ZoneBc * zonebc = zone->bccos[ ibcco ];
        Region region;
        region.SetRegion( zonebc->pnts );
        Boundary( region, zonebc->bcType );
    }
}

void WenoField::Boundary( Region &region, int bcType )
{
    if ( bcType == BCInflow )
    {
        this->InflowBc( region );
    }
    else if ( bcType == BCExtrapolate || bcType == BCOutflow )
    {
        this->OutflowBc( region );
    }
}

void WenoField::InflowBc( Region &region )
{
    int index_dim = region.start.size();
    if ( index_dim != 1 ) return;
    int st = region.start[ 0 ];
    int ed = region.end[ 0 ];
    for ( int i = st; i <= ed; ++ i )
    {
        int idir = 1;
        if ( i == 1 )
        {
            idir = -1;
        }
        int ib = i - 1; //index from 0
        int in = ib - idir;

        int ig1 = ib + idir;
        int ig2 = ig1 + idir;
        int ig3 = ig2 + idir;

        this->u[ ib ] = 0.0;
        this->u[ ig1 ] = 2.0 * u[ ib ] - 1.0 * u[ in ]; 
        this->u[ ig2 ] = 3.0 * u[ ib ] - 2.0 * u[ in ];   
        this->u[ ig3 ] = 4.0 * u[ ib ] - 3.0 * u[ in ];
    }
}

void WenoField::OutflowBc( Region &region )
{
    int index_dim = region.start.size();
    if ( index_dim != 1 ) return;
    int st = region.start[ 0 ];
    int ed = region.end[ 0 ];
    for ( int i = st; i <= ed; ++ i )
    {
        int idir = 1;
        if ( i == 1 )
        {
            idir = -1;
        }
        int ib = i - 1; //index from 0
        int in = ib - idir;

        int ig1 = ib + idir;
        int ig2 = ig1 + idir;
        int ig3 = ig2 + idir;

        this->u[ ib ] = 0.0;
        this->u[ ig1 ] = 2.0 * u[ ib ] - 1.0 * u[ in ]; 
        this->u[ ig2 ] = 3.0 * u[ ib ] - 2.0 * u[ in ];   
        this->u[ ig3 ] = 4.0 * u[ ib ] - 3.0 * u[ in ];
    }
}

void WenoField::RungeKutta( Zone * zone, int istage )
{
    if ( istage == 0 )
    {
        this->RungeKutta3Stage0( zone );
        return;
    }

    if ( istage == 1 )
    {
        this->RungeKutta3Stage1( zone );
        return;
    }

    if ( istage == 2 )
    {
        this->RungeKutta3Stage2( zone );
        return;
    }
}

void WenoField::RungeKutta3Stage0( Zone * zone )
{
    this->Rhs( this->u, this->res );

    for ( int i = 0; i < ni; ++ i )
    {
        u[ i ] = u[ i ] + dt * res[ i ];
    }
}

void WenoField::RungeKutta3Stage1( Zone * zone )
{
    this->Rhs( this->u, this->res );

    for ( int i = 0; i < ni; ++ i )
    {
        u[ i ] = 0.75 * un[ i ] + 0.25 * u[ i ] + 0.25 * dt * res[ i ];
    }
}

void WenoField::RungeKutta3Stage2( Zone * zone )
{
    this->Rhs( this->u, this->res );

    double c1 = 1.0 / 3.0;
    double c2 = 2.0 / 3.0;
    double c3 = 2.0 / 3.0;

    for ( int i = 0; i < ni; ++ i )
    {
        u[ i ] = c1 * un[ i ] + c2 * u[ i ] + c3 * dt * res[ i ];
    }
}

void WenoField::Rhs( Vec1d & u, Vec1d & res )
{
    Vec1d uL;
    uL.Allocate( 0, ni, 0 );

    Vec1d uR;
    uR.Allocate( 0, ni, 0 );
    if ( Global::scheme == Scheme::CRWENO )
    {
        crwenoL( ni, u, uL );
        crwenoR( ni, u, uR );
    }
    else
    {
        wenoL( ni, u, uL );
        wenoR( ni, u, uR );
    }

    for ( int i = 0; i < ni; ++ i )
    {
        if ( u[ i ] >= 0.0 )
        {
            res[ i ] = - u[ i ] * ( uL[ i + 1 ] - uL[ i ] ) / dx;
        }
        else
        {
            res[ i ] = - u[ i ] * ( uR[ i + 1 ] - uR[ i ] ) / dx;
        }
    }
}

void WenoField::UpdateOldField()
{
    this->un = this->u;
}

void WenoField::PostProcess( Grid * grid )
{
    this->DumpField( grid->x, u );
}

void WenoField::DumpField( Vec1d & x, Vec1d & u )
{
    for ( int i = 0; i < x.size(); ++ i )
    {
        Global::file_string += std::format( "{:.16f} {:.16f}\n", x[ i ], u[ i ] );
    }
}

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
