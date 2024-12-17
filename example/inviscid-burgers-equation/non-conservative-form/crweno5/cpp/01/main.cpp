#include <cmath>
#include <vector>
#include <print>
#include <numbers>
#include <fstream>

template<typename... Args>  
auto SQR(Args... args) {  
    return (... + (args * args));
}

void DumpCsvFile( const std::string & filename, std::vector<double> & x, std::vector<std::vector<double>> & u );
double wcL( double v1, double v2, double v3, double v4, double v5 );
double wcR( double v1, double v2, double v3, double v4, double v5 );
void rhs( int nx, double dx, std::vector<double> & u, std::vector<double> & r );
void numerical( int nx, int ns, int nt, double dx, double dt, std::vector<std::vector<double>> & u );

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
void wenoL( int nx, std::vector<double> & u, std::vector<double> & f )
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
void wenoR( int nx, std::vector<double> & u, std::vector<double> & f )
{
    int i;
    double v1, v2, v3, v4, v5;

    i = 1;
    v1 = 2.0 * u[ i - 1 ] - u[ i ];
    v2 = u[ i - 1 ];
    v3 = u[ i ];
    v4 = u[ i + 1 ];
    v5 = u[ i + 2 ];
    f[ i ] = wcR( v1, v2, v3, v4, v5 );
    for ( int i = 2; i < nx - 1; ++ i )
    {
        v1 = u[ i - 2 ];
        v2 = u[ i - 1 ];
        v3 = u[ i ];
        v4 = u[ i + 1 ];
        v5 = u[ i + 2 ];
        f[ i ] = wcR( v1, v2, v3, v4, v5 );
    }

    i = nx - 1;
    v1 = u[ i - 2 ];
    v2 = u[ i - 1 ];
    v3 = u[ i ];
    v4 = u[ i + 1 ];
    v5 = 2.0 * u[ i + 1 ] - u[ i ];
    f[ i ] = wcR( v1, v2, v3, v4, v5 );

    i = nx;
    v1 = u[ i - 2 ];
    v2 = u[ i - 1 ];
    v3 = u[ i ];
    v4 = 2.0 * u[ i ] - u[ i - 1 ];
    v5 = 3.0 * u[ i ] - 2.0 * u[ i - 1 ];
    f[ i ] = wcR( v1, v2, v3, v4, v5 );
}

//-----------------------------------------------------------------------------
// Calculate right hand term of the inviscid Burgers equation
//-----------------------------------------------------------------------------
void rhs( int nx, double dx, std::vector<double> & u, std::vector<double> & r )
{
    std::vector<double> uL( nx, 0 );
    std::vector<double> uR( nx + 1, 0 );

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
            r[ i ] = - u[ i ] * ( uR[ i + 1 ] - uR[ i ] ) / dx;
        }
    }
}

//-----------------------------------------------------------------------------
// Compute numerical solution
//   - Time integration using Runge-Kutta third order
//   - 5th-order Compact WENO scheme for spatial terms
//-----------------------------------------------------------------------------
void numerical( int nx, int ns, int nt, double dx, double dt, std::vector<std::vector<double>> &u )
{
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

void DumpCsvFile( const std::string &filename, std::vector<double> &x, std::vector<std::vector<double>> &u )
{
    std::fstream file;
    file.open(filename.c_str(), std::fstream::out);
    for ( int i = 0; i < u.size(); ++ i )
    {
        std::string str = {};
        str += std::format( "{:.16f} ", x[ i ] );
        int nj = u[ i ].size();
        for ( int j = 0; j < nj; ++ j )
        {
            double value = u[ i ][ j ];
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


int main( int argc, char ** argv )
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

    numerical( nx, ns, nt, dx, dt, u );

    return 0;
}
