#include "BurgersField.h"
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

void BurgersField::Init( Grid * grid )
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

    int ist = 0 - Global::nghost;
    int ied = this->ni - 1 + Global::nghost;

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

//void BurgersField::PhysicalBoundary( Zone * zone )
//{
//    int nbccos = zone->bccos.size();
//    for ( int ibcco = 0; ibcco < nbccos; ++ ibcco )
//    {
//        ZoneBc * zonebc = zone->bccos[ ibcco ];
//        Region region;
//        region.SetRegion( zonebc->pnts );
//        Boundary( region, zonebc->bcType );
//    }
//}
//
//void BurgersField::InterfaceBoundary( Zone * zone )
//{
//    int nbc1to1s = zone->bc1to1s.size();
//
//    for ( int ibc1to1 = 0; ibc1to1 < nbc1to1s; ++ ibc1to1 )
//    {
//        ZoneBc1To1 * bc1to1 = zone->bc1to1s[ ibc1to1 ];
//        Region region;
//        region.SetRegion( bc1to1->pnts );
//        this->InterfaceBc( region );
//    }
//}
//
//void BurgersField::Boundary( Region &region, int bcType )
//{
//    if ( bcType == BCInflow )
//    {
//        this->InflowBc( region );
//    }
//    else if ( bcType == BCExtrapolate || bcType == BCOutflow )
//    {
//        this->OutflowBc( region );
//    }
//}
//
//void BurgersField::InflowBc( Region &region )
//{
//    int index_dim = region.start.size();
//    if ( index_dim != 1 ) return;
//    int st = region.start[ 0 ];
//    int ed = region.end[ 0 ];
//    for ( int i = st; i <= ed; ++ i )
//    {
//        int idir = 1;
//        if ( i == 1 )
//        {
//            idir = -1;
//        }
//        int ib = i - 1; //index from 0
//        int in = ib - idir;
//
//        int ig1 = ib + idir;
//        this->u[ ib ] = 0.0;
//        this->u[ ig1 ] = 2.0 * u[ ib ] - 1.0 * u[ in ]; 
//
//        if ( Global::nghost >= 2 )
//        {
//            int ig2 = ig1 + idir;
//            this->u[ ig2 ] = 3.0 * u[ ib ] - 2.0 * u[ in ];
//            if ( Global::nghost >= 3 )
//            {
//                int ig3 = ig2 + idir;
//                this->u[ ig3 ] = 4.0 * u[ ib ] - 3.0 * u[ in ];
//            }
//        }
//    }
//}
//
//void BurgersField::OutflowBc( Region &region )
//{
//    int index_dim = region.start.size();
//    if ( index_dim != 1 ) return;
//    int st = region.start[ 0 ];
//    int ed = region.end[ 0 ];
//    for ( int i = st; i <= ed; ++ i )
//    {
//        int idir = 1;
//        if ( i == 1 )
//        {
//            idir = -1;
//        }
//        int ib = i - 1; //index from 0
//        int in = ib - idir;
//
//        int ig1 = ib + idir;
//        this->u[ ib ] = 0.0;
//        this->u[ ig1 ] = 2.0 * u[ ib ] - 1.0 * u[ in ]; 
//
//        if ( Global::nghost >= 2 )
//        {
//            int ig2 = ig1 + idir;
//            this->u[ ig2 ] = 3.0 * u[ ib ] - 2.0 * u[ in ];
//            if ( Global::nghost >= 3 )
//            {
//                int ig3 = ig2 + idir;
//                this->u[ ig3 ] = 4.0 * u[ ib ] - 3.0 * u[ in ];
//            }
//        }
//    }
//}
//
//void BurgersField::InterfaceBc( Region & region )
//{
//    //int index_dim = region.start.size();
//    //if ( index_dim != 1 ) return;
//    //int st = region.start[ 0 ];
//    //int ed = region.end[ 0 ];
//    //for ( int i = st; i <= ed; ++ i )
//    //{
//    //    int ib = i - 1; //index from 0
//
//    //    double value = 0.25 * ( 2 * this->u[ ib ] + this->u[ ib + 1 ] + this->u[ ib - 1 ] );
//    //    this->u[ ib ] = value;
//    //}
//}

void BurgersField::InviscidResidual( Vec1d & u, Vec1d & res )
{
    Vec1d uL;
    uL.Allocate( 0, ni, 0 );

    Vec1d uR;
    uR.Allocate( 0, ni, 0 );

    if ( Global::scheme.inviscid == to_int( BasicScheme::CRWENO ) )
    {
        crwenoL( ni, u, uL );
        crwenoR( ni, u, uR );
    }
    else if ( Global::scheme.inviscid == to_int( BasicScheme::WENO ) )
    {
        wenoL( ni, u, uL );
        wenoR( ni, u, uR );
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
            res[ i ] += ( - u[i] * ( uL[i + 1] - uL[i] ) / dx );
        }
        else
        {
            res[ i ] += ( - u[i] * ( uR[i + 1] - uR[i] ) / dx );
        }
    }
}
void BurgersField::ViscousResidual( Vec1d & u, Vec1d & res )
{
    ;
}

void BurgersField::Rhs( Vec1d & u, Vec1d & res )
{
    res = 0;
    InviscidResidual( u, res );
    if ( Global::iviscous > 0 )
    {
        ViscousResidual( u, res );
    }
 }

void BurgersField::UpdateOldField()
{
    this->un = this->u;
}

void BurgersField::PostProcess( Grid * grid )
{
    this->DumpField( grid->x, u );
}

void BurgersField::DumpField( Vec1d & x, Vec1d & u )
{
    for ( int i = 0; i < x.size(); ++ i )
    {
        Global::file_string += std::format( "{:.16f} {:.16f}\n", x[ i ], u[ i ] );
    }
}
