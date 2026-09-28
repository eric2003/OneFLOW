#include "Field.h"
#include "hxmath.h"
#include "Grid.h"
#include "Global.h"
#include "cgnslib.h"
#include <iostream>
#include <numbers>

void Field::CrankNicolsonSeries( Zone * zone )
{
    BasicScheme time_scheme = to_BasicScheme( Global::scheme.time_scheme );
    switch ( time_scheme ) {
    case BasicScheme::CN:
        this->CN( zone );
        break;
    case BasicScheme::ICP:
        this->ICP( zone );
        break;
    default:
        this->CN( zone );
    }
}

void Field::RungeKutta( Zone * zone, int nStage, int istage )
{
    if ( nStage == 1 )
    {
        this->RungeKutta1( zone, istage );
    }
    else if ( nStage == 3 )
    {
        this->RungeKutta3( zone, istage );
    }
}

void Field::RungeKutta1( Zone * zone, int istage )
{
    this->RungeKutta3Stage0( zone );
}

void Field::RungeKutta3( Zone * zone, int istage )
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
    this->Rhs( this->u, this->res );
    for ( int m = 0; m < nequ; ++ m )
    {
        Vec1d & u = this->u.vec( m );
        Vec1d & res = this->res.vec( m );
        for ( int i = 0; i < ni; ++ i )
        {
            u[ i ] = u[ i ] + dt * res[ i ];
        }
    }
}

void Field::RungeKutta3Stage1( Zone * zone )
{
    this->Rhs( this->u, this->res );
    for ( int m = 0; m < nequ; ++ m )
    {
        Vec1d & u = this->u.vec( m );
        Vec1d & un = this->un.vec( m );
        Vec1d & res = this->res.vec( m );
        for ( int i = 0; i < ni; ++ i )
        {
            u[ i ] = 0.75 * un[ i ] + 0.25 * u[ i ] + 0.25 * dt * res[ i ];
        }
    }
}

void Field::RungeKutta3Stage2( Zone * zone )
{
    this->Rhs( this->u, this->res );

    double c1 = 1.0 / 3.0;
    double c2 = 2.0 / 3.0;
    double c3 = 2.0 / 3.0;
    for ( int m = 0; m < nequ; ++ m )
    {
        Vec1d & u = this->u.vec( m );
        Vec1d & un = this->un.vec( m );
        Vec1d & res = this->res.vec( m );

        for ( int i = 0; i < ni; ++ i )
        {
            u[ i ] = c1 * un[ i ] + c2 * u[ i ] + c3 * dt * res[ i ];
        }
    }
}

void Field::PhysicalBoundary( Zone * zone )
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

void Field::InterfaceBoundary( Zone * zone )
{
    int nbc1to1s = zone->bc1to1s.size();

    for ( int ibc1to1 = 0; ibc1to1 < nbc1to1s; ++ ibc1to1 )
    {
        ZoneBc1To1 * bc1to1 = zone->bc1to1s[ ibc1to1 ];
        Region region;
        region.SetRegion( bc1to1->pnts );
        this->InterfaceBc( region );
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

    Vec1d & u = this->u.vec();

    for ( int i = st; i <= ed; ++ i )
    {
        int idir = 1;
        int ib = i - 1; //index from 0
        if ( i == 1 )
        {
            idir = -1;
        }
        else
        {
            ib -= Global::ifinite_volume;
        }
        int in = ib - idir;

        int ig1 = ib + idir;
        double ub = 0.0;
        double uin = u[ in ];

        if ( Global::ifinite_volume == 0 )
        {
            u[ ib ] = ub;
        }

        u[ ig1 ] = 2.0 * ub - 1.0 * uin; 

        if ( Global::nghost >= 2 )
        {
            int ig2 = ig1 + idir;
            u[ ig2 ] = 3.0 * ub - 2.0 * uin;
            if ( Global::nghost >= 3 )
            {
                int ig3 = ig2 + idir;
                u[ ig3 ] = 4.0 * ub - 3.0 * uin;
            }
        }
    }
}

void Field::OutflowBc( Region &region )
{
    int index_dim = region.start.size();
    if ( index_dim != 1 ) return;
    int st = region.start[ 0 ];
    int ed = region.end[ 0 ];
    Vec1d & u = this->u.vec();
    for ( int i = st; i <= ed; ++ i )
    {
        int idir = 1;
        int ib = i - 1; //index from 0
        if ( i == 1 )
        {
            idir = -1;
        }
        else
        {
            ib -= Global::ifinite_volume;
        }
        int in = ib - idir;

        int ig1 = ib + idir;
        double ub = 0.0;
        double uin = u[ in ];

        if ( Global::ifinite_volume == 0 )
        {
            u[ ib ] = ub;
        }

        u[ ig1 ] = 2.0 * ub - 1.0 * uin; 

        if ( Global::nghost >= 2 )
        {
            int ig2 = ig1 + idir;
            u[ ig2 ] = 3.0 * ub - 2.0 * uin;
            if ( Global::nghost >= 3 )
            {
                int ig3 = ig2 + idir;
                u[ ig3 ] = 4.0 * ub - 3.0 * uin;
            }
        }
    }
}

void Field::InterfaceBc( Region & region )
{
    //int index_dim = region.start.size();
    //if ( index_dim != 1 ) return;
    //int st = region.start[ 0 ];
    //int ed = region.end[ 0 ];
    //for ( int i = st; i <= ed; ++ i )
    //{
    //    int ib = i - 1; //index from 0

    //    double value = 0.25 * ( 2 * this->u[ ib ] + this->u[ ib + 1 ] + this->u[ ib - 1 ] );
    //    this->u[ ib ] = value;
    //}
}