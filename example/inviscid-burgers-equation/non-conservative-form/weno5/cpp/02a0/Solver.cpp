#include "Solver.h"

Solver::Solver()
{
    ;
}

Solver::~Solver()
{
    ;
}

void Solver::Solve()
{
    numerical();
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
void wenoL( int N, vec1d & u, vec1d & f )
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
void wenoR( int N, vec1d & u, vec1d & f )
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

void boundary( int N, vec1d & u )
{
    //left bc
    int i = 0;
    u[ i ] = 0.0;
    u[ i - 1 ] = 2.0 * u[ i ] - u[ i + 1 ]; 
    u[ i - 2 ] = 3.0 * u[ i ] - 2.0 * u[ i + 1 ];   
    u[ i - 3 ] = 4.0 * u[ i ] - 3.0 * u[ i + 1 ];

    //right bc
    i = N - 1;
    u[ i ] = 0.0;
    u[ i + 1 ] = 2.0 * u[ i ] - u[ i - 1 ];
    u[ i + 2 ] = 3.0 * u[ i ] - 2.0 * u[ i - 1 ];
    u[ i + 3 ] = 4.0 * u[ i ] - 3.0 * u[ i - 1 ];
}

//-----------------------------------------------------------------------------
// Calculate right hand term of the inviscid Burgers equation
//-----------------------------------------------------------------------------
void rhs( int N, double dx, vec1d & u, vec1d & res )
{
    vec1d uL;
    uL.Allocate( 0, N, 0 );

    vec1d uR;
    uR.Allocate( 0, N, 0 );

    wenoL( N, u, uL );
    wenoR( N, u, uR );

    for ( int i = 0; i < N; ++ i )
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

void runge_kutta( int N, vec1d & un, vec1d & ut, vec1d & res, double dx, double dt )
{
    runge_kutta_stage1( N, un, ut, res, dx, dt );
    runge_kutta_stage2( N, un, ut, res, dx, dt );
    runge_kutta_stage3( N, un, ut, res, dx, dt );
    un = ut;
}

void runge_kutta_stage1( int N, vec1d & un, vec1d & ut, vec1d & res, double dx, double dt )
{
    rhs( N, dx, ut, res );
    for ( int i = 0; i < N; ++ i )
    {
        ut[ i ] = un[ i ] + dt * res[ i ];
    }
    boundary( N, ut );
}

void runge_kutta_stage2( int N, vec1d & un, vec1d & ut, vec1d & res, double dx, double dt )
{
    rhs( N, dx, ut, res );

    for ( int i = 0; i < N; ++ i )
    {
        ut[ i ] = 0.75 * un[ i ] + 0.25 * ut[ i ] + 0.25 * dt * res[ i ];
    }
    boundary( N, ut );
}

void runge_kutta_stage3( int N, vec1d & un, vec1d & ut, vec1d & res, double dx, double dt )
{
    rhs( N, dx, ut, res );

    for ( int i = 0; i < N; ++ i )
    {
        ut[ i ] = ( 1.0 / 3.0 ) * un[ i ] + ( 2.0 / 3.0 ) * ut[ i ] + ( 2.0 / 3.0 ) * dt * res[ i ];
    }
    boundary( N, ut );
}

//-----------------------------------------------------------------------------
// Compute numerical solution
//   - Time integration using Runge-Kutta third order
//   - 5th-order Compact WENO scheme for spatial terms
//-----------------------------------------------------------------------------
void numerical()
{
    int N = 201;
    int ns = 10;
    double dt = 0.0001;
    double tm = 0.25;

    double dx = 1.0 / ( N - 1 );
    int nt = std::round( tm / dt );
    double ds = tm / ns;

    std::print( "Npoints={}\n", N );
    std::print( "ns={}\n", ns );
    std::print( "dt={}\n", dt );
    std::print( "tm={}\n", tm );
    std::print( "dx={}\n", dx );
    std::print( "nt={}\n", nt );

    vec1d x;
    x.Allocate( 0, N - 1, 0 ); // npoints = N = nx + 1

    int ighost = 3;
    int ist = 0 - ighost;
    int ied = N - 1 + ighost;

    std::vector<vec1d> ulist( ns + 1 );

    for ( int i = 0; i < ulist.size(); ++ i )
    {
        ulist[ i ].Allocate( ist, ied, 0 );
    }

    vec1d un;
    un.Allocate( ist, ied, 0 ); // numerical solsution at every time step

    vec1d ut;
    ut.Allocate( ist, ied, 0 ); // temporary array during RK3 integration

    vec1d res;
    res.Allocate( 0, N, 0 ); //N+1

    int k = 0; // record index
    int freq = std::round( nt / ns );
    std::print( "freq = {}\n", freq );

    for ( int i = 0; i < N; ++ i )
    {
        x[ i ] = dx * ( i );
        un[ i ] = std::sin( 2.0 * std::numbers::pi * x[ i ] );
    }

    boundary( N, un );
    ut = un;
    ulist[ k ] = un; //store solution at t = 0;

    for ( int it = 1; it <= nt; ++ it )
    {
        runge_kutta( N, un, ut, res, dx, dt );

        if ( it % freq == 0 )
        {
            k += 1;
            ulist[ k ] = un; 
            std::print( "k={}, ns={}\n", k, ns );
        }
    }

    DumpCsvFile( "field_final.csv", x, ulist );
}

void DumpCsvFile( const std::string &filename, vec1d &x, std::vector<vec1d> &u )
{
    std::fstream file;
    file.open( filename.c_str(), std::fstream::out );
    for ( int i = 0; i < x.size(); ++ i )
    {
        std::string str = {};
        str += std::format( "{:.16f} ", x[ i ] );
        int nj = u.size();
        for ( int j = 0; j < nj; ++ j )
        {
            double value = u[ j ][ i ];
            str += std::format( "{:.16f}", value );
            if ( j != nj - 1 )
            {
                str += " ";
            }
        }
        std::format_to(std::ostream_iterator<char>(file), "{}\n", str );
    }
    file.close();
}
