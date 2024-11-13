#include "Field.h"
#include "CgnsGrid.h"

int Para::nt = -1;
double Para::dx = 0;
double Para::dt = 0;
double Para::t = 0;
double Para::alpha = 0;
double Para::beta = 0;

void Para::Init()
{
    Para::dx = 0.025;
    Para::dt = dx / 10.0;
    Para::t = 1.0;
    Para::nt = std::round( t / dt );
    //Para::nt = 2;
    //Para::nt = 1;

}

void GenGrid( std::vector<double> & gx );
void GenGrid( std::vector<double> & gx )
{
    const int ni = 81;
    gx.resize( ni );

    double x_l = -1.0;
    double x_r = 1.0;

    double dx = ( x_r - x_l ) / ( ni - 1 );

    for ( int i = 0; i < ni; ++ i )
    {
        gx[ i ] = x_l + i * dx;
    }
}

void Field::ModifyGrid( Grid * grid )
{
    std::vector<double> globalx;
    GenGrid( globalx );

    this->ni = grid->x.size();

    int ishift = 0;

    if ( grid->zoneIndex == 1 )
    {
        ishift = ni - 1;
    }

    for ( int i = 0; i < ni; ++ i )
    {
        grid->x[ i ] = globalx[ i + ishift ];
    }
}

void Field::Init( Grid * grid )
{
    this->ni = grid->x.size();
    std::cout << "ni = " << ni << "\n";
    this->ModifyGrid( grid );

    std::vector<double> & x = grid->x;
    //this->dx = std::abs( x[ 1 ] - x[ 0 ] );
    //this->dt = dx / 10.0;
    //this->t = 1.0;
    //this->nt = std::round( t / dt );

    this->dx = Para::dx;
    this->dt = Para::dt;
    this->t = Para::t;
    this->nt = Para::nt;

    double total_time = nt * dt;

    std::cout << "this->dt = " << this->dt << "\n";
    std::cout << "this->t = " << this->t << "\n";
    std::cout << "this->nt = " << this->nt << "\n";
    std::cout << "this->ni = " << this->ni << "\n";
    std::cout << "total_time = " << std::setprecision( 25 ) << total_time << "\n";

    Global::nt = nt;

    this->alpha = 1 / ( std::numbers::pi * std::numbers::pi );
    this->beta = this->alpha * dt / ( dx * dx );
    std::cout << "alpha = " << std::setprecision( 25 ) << this->alpha << "\n";
    std::cout << "beta = " << std::setprecision( 25 ) << this->beta << "\n";

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
        u_e[ i ] = - std::exp( -total_time ) * std::sin( std::numbers::pi * xm ); //theory solution
        if ( std::abs( xm - 0.2 ) < 1.0e-5 )
        {
            double term1 = xm;
            double term2 = std::sin( std::numbers::pi * xm );
            double term3 = std::exp( -total_time );
            double term4 = term2 * term3;
            std::cout << "this->dx = " << std::setprecision( 25 ) << this->dx << "\n";
            std::cout << "this->dt = " << std::setprecision( 25 ) << this->dt << "\n";
            std::cout << "total_time = " << std::setprecision( 25 ) << total_time << "\n";
            std::cout << "term1 = " << std::setprecision( 25 ) << term1 << "\n";
            std::cout << "term2 = " << std::setprecision( 25 ) << term2 << "\n";
            std::cout << "term3 = " << std::setprecision( 25 ) << term3 << "\n";
            std::cout << "term4 = " << std::setprecision( 25 ) << term4 << "\n";
        }
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
        double term1 = un[ i - 1 ];
        double term2 = un[ i ];
        double term3 = un[ i + 1 ];
        double term4 = term1 -2.0 * term2 + term3;
        double term5 = beta * term4;
        double term6 = term2 + term5;
        if ( ( zone->zoneIndex == 0 && i== ied ) || 
             ( zone->zoneIndex == 1 && i== ist ) )
        {
            std::cout << "zone->zoneIndex = " << zone->zoneIndex << "\n";
            std::cout << "term1 = " << std::setprecision( 25 ) << term1 << "\n";
            std::cout << "term2 = " << std::setprecision( 25 ) << term2 << "\n";
            std::cout << "term3 = " << std::setprecision( 25 ) << term3 << "\n";
            std::cout << "term4 = " << std::setprecision( 25 ) << term4 << "\n";
            std::cout << "term5 = " << std::setprecision( 25 ) << term5 << "\n";
            std::cout << "term6 = " << std::setprecision( 25 ) << term6 << "\n";
        }

        //u[ i ] = un[ i ] + beta * ( un[ i + 1 ] - 2.0 * un[ i ] + un[ i - 1 ] );
        u[ i ] = term6;
    }

    //this->PhysicalBoundary( zone );

    //this->Update( un, u );
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
        //this->u[ i ] = 0.0;
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
        //this->u[ i ] = 0.0;
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

void Field::UpdateOldField()
{
    for ( int i = 0; i < this->u.size(); ++ i )
    {
        this->un[ i ] = this->u[ i ];
    }
}

void Field::Update( std::vector<double> &un, std::vector<double> &u )
{
    for ( int i = 0; i < u.size(); ++ i )
    {
        un[ i ] = u[ i ];
    }
}