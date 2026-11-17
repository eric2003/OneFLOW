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
    this->nt = std::round(  Global::total_time / dt );
    this->ns = 10;

    std::print( "ni={}\n", ni );
    std::print( "ns={}\n", ns );
    std::print( "dt={}\n", dt );
    std::print( "dx={}\n", dx );
    std::print( "nt={}\n", nt );
    std::cout << "this->dt   = " << this->dt << "\n";
    std::cout << "this->nt   = " << this->nt << "\n";
    std::cout << "this->ni   = " << this->ni << "\n";
    std::cout << "nt * dt = " << nt * dt << "\n";

    int freq = std::round( nt / ns );
    std::print( "freq = {}\n", freq );
    Global::nt = nt;

    int ist = 0 - Global::nghost;
    int ied = this->ni - 1 + Global::nghost;

    un.Allocate( ist, ied, 0 ); // numerical solsution at every time step
    u.Allocate( ist, ied, 0 ); // temporary array during RK3 integration
    res.Allocate( 0, this->ni, 0 ); //N+1

    for ( int i = 0; i < ni; ++ i )
    {
        u[ i ] = std::sin( 2.0 * std::numbers::pi * x[ i ] );
    }

    int kkk = 1;
}

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

    if ( Global::iconservation == 0 )
    {
        for ( int i = 0; i < ni; ++ i )
        {
            if ( u[ i ] >= 0.0 )
            {
                res[ i ] += ( - u[ i ] * ( uL[ i + 1 ] - uL[ i ] ) / dx );
            }
            else
            {
                res[ i ] += ( - u[ i ] * ( uR[ i + 1 ] - uR[ i ] ) / dx );
            }
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

void BurgersField::DumpField( Grid * grid )
{
    this->DumpField( grid->x, u );
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
