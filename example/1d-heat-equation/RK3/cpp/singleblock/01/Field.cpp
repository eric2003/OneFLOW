#include "Field.h"
#include "CgnsGrid.h"

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

    u_e.resize( ni_total );
    u.resize( ni_total );
    un.resize( ni_total );
    u1.resize( ni_total );
    r.resize( ni_total );

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

void Field::Solve( Zone * zone )
{
    //this->FTCS( zone );
    this->RungeKutta3( zone );
}

void Field::FTCS( Zone * zone )
{
    int ist = 1;
    int ied = ni;

    this->Rhs( this->un, this->r );

    for ( int i = ist; i <= ied; ++ i )
    {
        u[ i ] = un[ i ] + dt * r[ i ];
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

void Field::RungeKutta3Stage0( Zone * zone )
{
    int ist = 1;
    int ied = ni;

    this->Rhs( this->un, this->r );

    for ( int i = ist; i <= ied; ++ i )
    {
        u[ i ] = un[ i ] + dt * r[ i ];
    }
}

void Field::UpdateOldField()
{
    for ( int i = 0; i < this->u.size(); ++ i )
    {
        this->un[ i ] = this->u[ i ];
    }
}

void Field::UpdateRungeKuttaOldField( int istage )
{
    if ( istage == 0 )
    {
        for ( int i = 0; i < this->u.size(); ++ i )
        {
            this->u1[ i ] = this->u[ i ];
        }
        return;
    }

    if ( istage == 1 )
    {
        return;
    }

    if ( istage == 2 )
    {
        for ( int i = 0; i < this->u.size(); ++ i )
        {
            this->un[ i ] = this->u[ i ];
        }
        return;
    }
}

void Field::RungeKutta3Stage1( Zone * zone )
{
    int ist = 1;
    int ied = ni;

    this->Rhs( this->u1, this->r );

    for ( int i = ist; i <= ied; ++ i )
    {
        u[ i ] = 0.75 * un[ i ] + 0.25 * u1[ i ] + 0.25 * dt * r[ i ];
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


void Field::RungeKutta3( Zone * zone )
{
    int ist = 1;
    int ied = ni;

    this->Rhs( this->un, this->r );

    for ( int i = ist; i <= ied; ++ i )
    {
        u1[ i ] = un[ i ] + dt * r[ i ];
    }

    this->Rhs( this->u1, this->r );

    for ( int i = ist; i <= ied; ++ i )
    {
        u[ i ] = 0.75 * un[ i ] + 0.25 * u1[ i ] + 0.25 * dt * r[ i ];
    }

    this->Rhs( this->u, this->r );

    double c1 = 1.0 / 3.0;
    double c2 = 2.0 / 3.0;
    double c3 = 2.0 / 3.0;

    for ( int i = ist; i <= ied; ++ i )
    {
        u[ i ] = c1 * un[ i ] + c2 * u[ i ] + c3 * dt * r[ i ];
    }

}

void Field::Rhs( std::vector<double> & u, std::vector<double> & r )
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
        this->u[ ighost ] = 2 * this->u[ i ] - this->u[ iinner ];
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
        this->u[ ighost ] = 2 * this->u[ i ] - this->u[ iinner ];
    }
}


void Field::PostProcess( Grid * grid )
{
    std::vector<double> & x = grid->x;
    //compute L2 norm of the error
    std::vector<double> u_error( ni );
    for ( int i = 0; i < ni; ++ i )
    {
        u_error[ i ] = un[ i ] - u_e[ i ];
    }
}
