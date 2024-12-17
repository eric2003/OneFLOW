#include "Field.h"
#include "Grid.h"
#include "Global.h"
#include "cgnslib.h"
#include <iostream>
#include <numbers>

void thomas_algorithm( const std::vector<double> & a,
    const std::vector<double> & b,
    const std::vector<double> & c,
    const std::vector<double> & d,
    std::vector<double> & x )
{
    size_t N = d.size();

    std::vector<double> c_star( N, 0.0 );
    std::vector<double> d_star( N, 0.0 );

    c_star[ 0 ] = c[ 0 ] / b[ 0 ];
    d_star[ 0 ] = d[ 0 ] / b[ 0 ];

    for ( int i = 1; i < N; ++ i )
    {
        double coef = 1.0 / ( b[ i ] - a[ i ] * c_star[ i - 1 ] );
        c_star[ i ] = c[ i ] * coef;
        d_star[ i ] = ( d[ i ] - a[ i ] * d_star[ i - 1 ] ) * coef;
    }

    x[ N - 1 ] = d_star[ N - 1 ];

    for ( int i = N - 2; i >= 0; -- i )
    {
        x[ i ] = d_star[ i ] - c_star[ i ] * x[ i + 1 ];
    }
}

void Field::Init( Grid * grid )
{
    this->ni = grid->x.size();
    std::cout << "ni = " << ni << "\n";

    std::vector<double> & x = grid->x;
    this->dx = std::abs( x[ 1 ] - x[ 0 ] );
    this->dt = dx / 10.0;
    this->t = 1.0;
    this->nt = std::round( t / dt );

    double total_time = nt * dt;

    std::cout << "this->dt   = " << this->dt << "\n";
    std::cout << "this->t    = " << this->t << "\n";
    std::cout << "this->nt   = " << this->nt << "\n";
    std::cout << "this->ni   = " << this->ni << "\n";
    std::cout << "total_time = " << total_time << "\n";

    Global::nt = nt;

    this->alpha = 1 / ( std::numbers::pi * std::numbers::pi );
    this->beta = this->alpha * dt / ( dx * dx );

    int nghost = 2;
    int ni_total = ni + nghost;

    int N = this->ni;

    u_e.resize( ni_total );
    u.Allocate( 0, ni_total - 1 );
    un.Allocate( 0, ni_total - 1 );
    r.Allocate( 0, ni_total - 1 );

    // 0 1 2 ... ni-1 ni ni+1
    int ist = 1;
    int ied = ni;

    for ( int i = ist; i <= ied; ++ i )
    {
        double xm = x[ i - ist ];
        u_e[ i ] = - std::exp( -total_time ) * std::sin( std::numbers::pi * xm ); //theory solution
        u[ i ] = - std::sin( std::numbers::pi * xm ); //initial condition @ t=0
    }
    int kkk = 1;
}

void Field::FTCS( Zone * zone )
{
    int ist = 1;
    int ied = ni;

    this->Rhs( this->u, this->r );

    for ( int i = ist; i <= ied; ++ i )
    {
        u[ i ] = u[ i ] + dt * r[ i ];
    }
}

void Field::CN( Zone * zone )
{
    int ist = 1;
    int ied = ni;

    double rr = 0.5 * this->alpha * dt / ( dx * dx );
    std::vector<double> a( ni );//0:ni-1
    std::vector<double> b( ni );//0:ni-1
    std::vector<double> c( ni );//0:ni-1
    std::vector<double> d( ni );//0:ni-1

    for ( int i = 0; i < ni; ++ i )
    {
        a[ i ] = - rr;
        b[ i ] = 1.0 + 2.0 * rr;
        c[ i ] = - rr;
    }

    for ( int i = ist; i <= ied; ++ i )
    {
        int ii = i - 1;
        d[ ii ] = rr * u[ i - 1 ] + ( 1.0 - 2.0 * rr ) * u[ i ] + rr * u[ i + 1 ];
    }

    a[ 0 ] = 0;
    c[ ni - 1 ] = 0;

    d[ 0 ] -= ( - rr ) * u[ 0 ];
    d[ ni - 1 ] -= ( - rr ) * u[ ni + 1 ];

    std::vector<double> values( d.size() );

    thomas_algorithm( a, b, c, d, values );

    for ( int i = ist; i <= ied; ++ i )
    {
        int ii = i - 1;
        u[ i ] = values[ ii ];
    }
}

void Field::ICP( Zone * zone )
{
    int ist = 1;
    int ied = ni;

    double rr = 0.5 * this->alpha * dt / ( dx * dx );
    std::vector<double> a( ni );//0:ni-1
    std::vector<double> b( ni );//0:ni-1
    std::vector<double> c( ni );//0:ni-1
    std::vector<double> d( ni );//0:ni-1

    for ( int i = 0; i < ni; ++ i )
    {
        a[ i ] = 1.0 / 12.0 - rr;
        b[ i ] = 10.0 / 12.0 + 2.0 * rr;
        c[ i ] = 1.0 / 12.0 - rr;
    }

    a[ 0 ] = 0;
    b[ 0 ] = 1;
    c[ 0 ] = 0;

    a[ ni - 1 ] = 0;
    b[ ni - 1 ] = 1;
    c[ ni - 1 ] = 0;

    for ( int i = ist; i <= ied; ++ i )
    {
        double aa = 1.0 / 12.0 + rr;
        double bb = 10.0 / 12.0 - 2.0 * rr;
        double cc = 1.0 / 12.0 + rr;
        int ii = i - 1;
        d[ ii ] = aa * u[ i - 1 ] + bb * u[ i ] + cc * u[ i + 1 ];
    }

    d[ 0 ] = 0;
    d[ ni - 1 ] = 0;

    std::vector<double> values( d.size() );

    thomas_algorithm( a, b, c, d, values );

    for ( int i = ist; i <= ied; ++ i )
    {
        int ii = i - 1;
        u[ i ] = values[ ii ];
    }
}

void Field::RungeKutta( Zone * zone, int istage )
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

void Field::UpdateOldField()
{
    this->un = this->u;
}

void Field::RungeKutta3Stage0( Zone * zone )
{
    int ist = 1;
    int ied = ni;

    this->Rhs( this->u, this->r );

    for ( int i = ist; i <= ied; ++ i )
    {
        u[ i ] = u[ i ] + dt * r[ i ];
    }
}

void Field::RungeKutta3Stage1( Zone * zone )
{
    int ist = 1;
    int ied = ni;

    this->Rhs( this->u, this->r );

    for ( int i = ist; i <= ied; ++ i )
    {
        u[ i ] = 0.75 * un[ i ] + 0.25 * u[ i ] + 0.25 * dt * r[ i ];
    }
}

void Field::RungeKutta3Stage2( Zone * zone )
{
    int ist = 1;
    int ied = ni;

    this->Rhs( this->u, this->r );

    double c1 = 1.0 / 3.0;
    double c2 = 2.0 / 3.0;
    double c3 = 2.0 / 3.0;

    for ( int i = ist; i <= ied; ++ i )
    {
        u[ i ] = c1 * un[ i ] + c2 * u[ i ] + c3 * dt * r[ i ];
    }
}


void Field::Rhs( Vec1d & u, Vec1d & r )
{
    int ist = 1;
    int ied = ni;

    double coef = this->alpha / ( dx * dx );

    for ( int i = ist; i <= ied; ++ i )
    {
        r[ i ] = coef * ( u[ i + 1 ] - 2.0 * u[ i ] + u[ i - 1 ] );
    }
}

void Field::PhysicalBoundary( Zone * zone )
{
    int nbccos = zone->bccos.size();
    for ( int ibcco = 0; ibcco < nbccos; ++ ibcco )
    {
        ZoneBc * zonebc = zone->bccos[ ibcco ];
        //std::cout << "zonebc->bcType = " << zonebc->bcType << "\n";
        Region region;
        region.SetRegion( zonebc->pnts );
        //region.Print();
        Boundary( region, zonebc->bcType );
    }
}

void Field::Boundary( Region &region, int bcType )
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

void Field::InflowBc( Region &region )
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
        int ighost = i + idir;
        int iinner = i - idir;
        //this->u[ ighost ] = 2 * this->u[ i ] - this->u[ iinner ];
        this->u[ ighost ] = - this->u[ iinner ];
        int kkk = 1;
    }
}

void Field::OutflowBc( Region &region )
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
        int ighost = i + idir;
        int iinner = i - idir;
        //this->u[ ighost ] = 2 * this->u[ i ] - this->u[ iinner ];
        this->u[ ighost ] = - this->u[ iinner ];
    }
}
