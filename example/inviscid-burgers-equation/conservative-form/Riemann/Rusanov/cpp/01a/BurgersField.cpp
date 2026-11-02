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
    this->ni = grid->ni;
    this->nic = grid->nic;
    grid->CalcMetrics();
    if ( Global::ifinite_volume == 1 )
    {
        this->nx = this->nic;
    }
    else
    {
        // finite difference
        this->nx = this->ni;
    }

    Vec1d & x = grid->x;
    this->dx = std::abs( x[ 1 ] - x[ 0 ] );
    this->dt = 0.0001;
    this->nt = std::round(  Global::total_time / dt );

    std::print( "ni={}\n", ni );
    std::print( "dt={}\n", dt );
    std::print( "dx={}\n", dx );
    std::print( "nt={}\n", nt );
    std::cout << "this->dt   = " << this->dt << "\n";
    std::cout << "this->nt   = " << this->nt << "\n";
    std::cout << "this->ni   = " << this->ni << "\n";
    std::cout << "nt * dt = " << nt * dt << "\n";

    Global::nt = nt;

    int ist = 0 - Global::nghost;
    int ied = this->nx - 1 + Global::nghost;

    un.Allocate( ist, ied, 0 ); // numerical solsution at every time step
    u.Allocate( ist, ied, 0 ); // temporary array during RK3 integration
    res.Allocate( 0, this->nx, 0 ); //N+1

    if ( Global::ifinite_volume == 0 )
    {
        //node
        for ( int i = 0; i < ni; ++ i )
        {
            u[ i ] = std::sin( 2.0 * std::numbers::pi * x[ i ] );
        }
    }
    else
    {
        //cell center
        Vec1d & xcc = grid->xcc;
        for ( int i = 0; i < nic; ++ i )
        {
            u[ i ] = std::sin( 2.0 * std::numbers::pi * xcc[ i ] );
        }

    }


    int kkk = 1;
}

void BurgersField::InviscidResidual( Vec1d & u, Vec1d & res )
{
    if ( Global::iconservation == 0 )
    {
        this->InviscidNonConservative( u, res );
    }
    else
    {
        this->InviscidConservative( u, res );
    }
}

void BurgersField::InviscidNonConservative( Vec1d & u, Vec1d & res )
{
    Vec1d uL;
    uL.Allocate( 0, nx, 0 );

    Vec1d uR;
    uR.Allocate( 0, nx, 0 );

    if ( Global::scheme.reconstruction == to_int( BasicScheme::CRWENO ) )
    {
        crwenoL( nx, u, uL );
        crwenoR( nx, u, uR );
    }
    else if ( Global::scheme.reconstruction == to_int( BasicScheme::WENO ) )
    {
        wenoL( nx, u, uL );
        wenoR( nx, u, uR );
    }

    if ( Global::scheme.inviscid == to_int( BasicScheme::CENTER ) )
    {
        for ( int i = 0; i < nx; ++ i )
        {
            res[ i ] += ( - u[ i ] * ( u[ i + 1 ] - u[ i - 1 ] ) / dx );
        }
    }
    else
    {
        for ( int i = 0; i < nx; ++ i )
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

void BurgersField::WaveSpeed( Vec1d & u, Vec1d & ps )
{
    for ( int i = 0; i <= nx; ++ i )
    {
        ps[ i ] = std::max( { std::abs( u[ i - 2 ] ), std::abs( u[ i - 1 ] ), std::abs( u[ i ] ), std::abs( u[ i + 1 ] ), std::abs( u[ i + 2 ] ) } );
    }

    for ( int i = ps.ist; i < 0; ++ i )
    {
        ps[ i ] = ps[ 0 ];
    }

    for ( int i = nx + 1; i <= ps.ied; ++ i )
    {
        ps[ i ] = ps[ nx ];
    }
}

void BurgersField::LaxFriedrichs( Vec1d & u, Vec1d & res )
{
    Vec1d fL;
    fL.Allocate( 0, nx, 0 );

    Vec1d fR;
    fR.Allocate( 0, nx, 0 );

    Vec1d f, fP, fN;

    int ist = 0 - Global::nghost;
    int ied = this->nx - 1 + Global::nghost;

    f.Allocate( ist, ied, 0 );
    fP.Allocate( ist, ied, 0 );
    fN.Allocate( ist, ied, 0 );

    burgers_fluxes( ist, ied, u, f );

    Vec1d ps;
    ps.Allocate( ist, ied, 0 );

    WaveSpeed( u, ps );

    // left and right side fluxes at the interface
    for ( int i = ist; i <= ied; ++ i )
    {
        fP[ i ] = 0.5 * ( f[ i ] + ps[ i ] * u[ i ] );
        fN[ i ] = 0.5 * ( f[ i ] - ps[ i ] * u[ i ] );
    }

    if ( Global::scheme.reconstruction == to_int( BasicScheme::CRWENO ) )
    {
        crwenoL( nx, fP, fL );
        crwenoR( nx, fN, fR );
    }
    else if ( Global::scheme.reconstruction == to_int( BasicScheme::WENO ) )
    {
        wenoL( nx, fP, fL );
        wenoR( nx, fN, fR );
    }

    for ( int i = 0; i < nx; ++ i )
    {
        res[ i ] -= ( fL[ i + 1 ] - fL[ i ] ) / dx + ( fR[ i + 1 ] - fR[ i ] ) / dx;
    }
}

void BurgersField::burgers_fluxes( int ist, int ied, Vec1d & u, Vec1d & f )
{
    for ( int i = ist; i <= ied; ++ i )
    {
        f[ i ] = 0.5 * u[ i ] * u[ i ];
    }
}

void BurgersField::Rusanov( Vec1d & u, Vec1d & res )
{
    Vec1d uL, uR;
    uL.Allocate( 0, nx, 0 );
    uR.Allocate( 0, nx, 0 );

    //WENO Reconstruction
    if ( Global::scheme.reconstruction == to_int( BasicScheme::CRWENO ) )
    {
        crwenoL( nx, u, uL );
        crwenoR( nx, u, uR );
    }
    else if ( Global::scheme.reconstruction == to_int( BasicScheme::WENO ) )
    {
        wenoL( nx, u, uL );
        wenoR( nx, u, uR );
    }

    //left and right side fluxes at the interface
    Vec1d fL, fR;
    fL.Allocate( 0, nx, 0 );
    fR.Allocate( 0, nx, 0 );

    int ist = 0 - Global::nghost;
    int ied = this->nx - 1 + Global::nghost;

    //Computing fluxes
    burgers_fluxes( 0, nx, uL, fL );
    burgers_fluxes( 0, nx, uR, fR );

    Vec1d ps;
    ps.Allocate( ist, ied, 0 );

    WaveSpeed( u, ps );

    //fluxes at the interface
    Vec1d f;
    f.Allocate( 0, nx, 0 );

    //Interface fluxes (Rusanov)
    for ( int i = 0; i <= nx; ++ i )
    {
        f[ i ] = 0.5 * ( fR[ i ] + fL[ i ] ) - 0.5 * ps[ i ] * ( uR[ i ] - uL[ i ] );
    }

    for ( int i = 0; i < nx; ++ i )
    {
        res[ i ] -= ( f[ i + 1 ] - f[ i ] ) / dx;
    }
}

void BurgersField::InviscidConservative( Vec1d & u, Vec1d & res )
{
    if ( Global::scheme.inviscid == to_int( BasicScheme::LAX ) )
    {
        this->LaxFriedrichs( u, res );
    }
    else if ( Global::scheme.inviscid == to_int( BasicScheme::Rusanov ) )
    {
        this->Rusanov( u, res );
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
    if ( Global::ifinite_volume == 1 )
    {
        this->DumpField( grid->xcc, u );
    }
    else
    {
        this->DumpField( grid->x, u );
    }
    
}

void BurgersField::PostProcess( Grid * grid )
{
    if ( Global::ifinite_volume == 1 )
    {
        this->DumpField( grid->xcc, u );
    }
    else
    {
        this->DumpField( grid->x, u );
    }
}

void BurgersField::DumpField( Vec1d & x, Vec1d & u )
{
    for ( int i = 0; i < x.size(); ++ i )
    {
        Global::file_string += std::format( "{:.16f} {:.16f}\n", x[ i ], u[ i ] );
    }
}
