#include "Field.h"
#include "CgnsGrid.h"

void Field::Init( Grid * grid )
{
    this->ni = grid->x.size();
    std::cout << "ni = " << ni << "\n";

    std::vector<double> & x = grid->x;
    this->dx = x[ 1 ] - x[ 0 ];
    this->dt = dx / 10.0;
    this->t = 1.0;
    this->nt = static_cast<int>( t / dt ) + 1;
    std::cout << "nt = " << nt << "\n";

    Global::nt = nt;

    this->alpha = 1 / ( std::numbers::pi * std::numbers::pi );
    this->beta = this->alpha * dt / ( dx * dx );
    std::cout << "alpha = " << std::setprecision( 15 ) << this->alpha << "\n";
    std::cout << "beta = " << std::setprecision( 15 ) << this->beta << "\n";

    int nghost = 2;
    int ni_total = ni + nghost;

    u_e.resize( ni_total );
    u.resize( ni_total );
    un.resize( ni_total );

    // 0 1 2 ... ni-1 ni ni+1

    int ist = 1;
    int ied = ni;

    // 0(bc) 1(ist) 2 3 ... ni

    for ( int i = ist; i <= ied; ++ i )
    {
        double xm = x[ i - ist ];
        u_e[ i ] = - std::exp( -t ) * std::sin( std::numbers::pi * xm ); //theory solution
        u[ i ] = - std::sin( std::numbers::pi * xm ); //initial condition @ t=0
    }
}

void Field::Solve( Zone * zone )
{
    int nghost = 2;
    int ni_total = ni + nghost;

    int ist = 1;
    int ied = ni;

    for ( int i = ist; i <= ied; ++ i )
    {
        u[ i ] = un[ i ] + beta * ( un[ i + 1 ] - 2.0 * un[ i ] + un[ i - 1 ] );
    }

    this->PhysicalBoundary( zone );

    this->Update( un, u );
}

void Field::PhysicalBoundary( Zone * zone )
{
    int nbccos = zone->bccos.size();
    for ( int ibcco = 0; ibcco < nbccos; ++ ibcco )
    {
        ZoneBc * zonebc = zone->bccos[ ibcco ];
        std::cout << "zonebc->bcType = " << zonebc->bcType << "\n";
        Region region;
        region.SetRegion( zonebc->pnts );
        region.Print();
        Boundary( region, zonebc->bcType );
    }
}

void Field::Boundary( Region &region, int bcType )
{
    if ( bcType == BCInflow )
    {
        this->InflowBc( region );
    }
    else if ( bcType == BCExtrapolate )
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

void Field::AddData( PointFactory &ptfactory, Grid * grid, std::vector<double> &global_x, std::vector<double> &global_ue, std::vector<double> &global_un )
{
    int ni = grid->x.size();
    for ( int i = 0; i < ni; ++ i )
    {
        double xm = grid->x[ i ];
        bool flag = ptfactory.FindPoint( Point( xm ) );
        if ( i == 0 || ( i == ni - 1 ) )
        {
            if ( ! flag )
            {
                ptfactory.AddPoint( Point( xm ) );
                global_x.push_back( grid->x[ i ] );
                global_ue.push_back( this->u_e[ i ] );
                global_un.push_back( this->un[ i ] );
            }
        }
        else
        {
            global_x.push_back( grid->x[ i ] );
            global_ue.push_back( this->u_e[ i ] );
            global_un.push_back( this->un[ i ] );
        }
    }
}

void Field::Update( std::vector<double> &un, std::vector<double> &u )
{
    for ( int i = 0; i < u.size(); ++ i )
    {
        un[ i ] = u[ i ];
    }
}