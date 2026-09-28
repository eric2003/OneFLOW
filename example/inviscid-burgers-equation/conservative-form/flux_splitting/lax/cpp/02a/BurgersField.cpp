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

    this->u.Allocate( this->nequ, ist, ied );
    this->un.Allocate( this->nequ, ist, ied );
    this->res.Allocate( this->nequ, 0, this->nx ); //N+1

    Vec1d &u = this->u.vec();

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

void BurgersField::InviscidResidual( VecWrap & u, VecWrap & res )
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

void BurgersField::InviscidNonConservative( VecWrap & u, VecWrap & res )
{
    VecWrap uL, uR;
    uL.Allocate( this->nequ, 0, nx, 0 );
    uR.Allocate( this->nequ, 0, nx, 0 );

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
        for ( int m = 0; m < nequ; ++ m )
        {
            Vec1d & U = u.vec( m );
            Vec1d & Res = res.vec( m );
            for ( int i = 0; i < nx; ++ i )
            {
                Res[ i ] += ( - U[ i ] * ( U[ i + 1 ] - U[ i - 1 ] ) / dx );
            }
        }
    }
    else
    {
        for ( int m = 0; m < nequ; ++ m )
        {
            Vec1d & U = u.vec( m );
            Vec1d & Res = res.vec( m );
            Vec1d & UL = uL.vec( m );
            Vec1d & UR = uR.vec( m );
            for ( int i = 0; i < nx; ++ i )
            {
                if ( U[ i ] >= 0.0 )
                {
                    Res[ i ] += ( - U[ i ] * ( UL[ i + 1 ] - UL[ i ] ) / dx );
                }
                else
                {
                    Res[ i ] += ( - U[ i ] * ( UR[ i + 1 ] - UR[ i ] ) / dx );
                }
            }
        }
    }
}

void BurgersField::WaveSpeed( VecWrap & um, VecWrap & psm )
{
    for ( int m = 0; m < nequ; ++ m )
    {
        Vec1d & u = um.vec( m );
        Vec1d & ps = psm.vec( m );

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
}

void BurgersField::LaxFriedrichs( VecWrap & u, VecWrap & res )
{
    VecWrap fL, fR;
    fL.Allocate( this->nequ, 0, nx, 0 );
    fR.Allocate( this->nequ, 0, nx, 0 );

    int ist = 0 - Global::nghost;
    int ied = this->nx - 1 + Global::nghost;

    VecWrap f, fP, fN;

    f.Allocate( this->nequ, ist, ied, 0 );
    fP.Allocate( this->nequ, ist, ied, 0 );
    fN.Allocate( this->nequ, ist, ied, 0 );

    burgers_fluxes( ist, ied, u, f );

    VecWrap psm;
    psm.Allocate( this->nequ, ist, ied, 0 );

    WaveSpeed( u, psm );

    // left and right side fluxes at the interface
    for ( int m = 0; m < nequ; ++ m )
    {
        Vec1d & u = this->u.vec( m );
        Vec1d & FP = fP.vec( m );
        Vec1d & FN = fN.vec( m );
        Vec1d & F = f.vec( m );
        Vec1d & ps = psm.vec( m );
        for ( int i = ist; i <= ied; ++ i )
        {
            FP[ i ] = 0.5 * ( F[ i ] + ps[ i ] * u[ i ] );
            FN[ i ] = 0.5 * ( F[ i ] - ps[ i ] * u[ i ] );
        }
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

    for ( int m = 0; m < nequ; ++ m )
    {
        Vec1d & FL = fL.vec( m );
        Vec1d & FR = fR.vec( m );
        Vec1d & res = this->res.vec( m );
        for ( int i = 0; i < nx; ++ i )
        {
            res[ i ] -= ( FL[ i + 1 ] - FL[ i ] ) / dx + ( FR[ i + 1 ] - FR[ i ] ) / dx;
        }
    }
}

void BurgersField::burgers_fluxes( int ist, int ied, VecWrap & u, VecWrap & f )
{
    Vec1d & F = f.vec();
    Vec1d & U = u.vec();
    for ( int i = ist; i <= ied; ++ i )
    {
        F[ i ] = 0.5 * U[ i ] * U[ i ];
    }
}

void BurgersField::Rusanov( VecWrap & u, VecWrap & res )
{
    VecWrap uL, uR;
    uL.Allocate( this->nequ, 0, nx, 0 );
    uR.Allocate( this->nequ, 0, nx, 0 );

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
    VecWrap fL, fR;
    fL.Allocate( this->nequ, 0, nx, 0 );
    fR.Allocate( this->nequ, 0, nx, 0 );

    int ist = 0 - Global::nghost;
    int ied = this->nx - 1 + Global::nghost;

    //Computing fluxes
    burgers_fluxes( 0, nx, uL, fL );
    burgers_fluxes( 0, nx, uR, fR );

    VecWrap psm;
    psm.Allocate( this->nequ, ist, ied, 0 );

    WaveSpeed( u, psm );

    //fluxes at the interface
    Vec1d f;
    f.Allocate( 0, nx, 0 );

    //Interface fluxes (Rusanov)
    for ( int m = 0; m < nequ; ++ m )
    {
        Vec1d & u = this->u.vec( m );
        Vec1d & res = this->res.vec( m );
        Vec1d & qL = uL.vec( m );
        Vec1d & qR = uR.vec( m );
        Vec1d & FL = fL.vec( m );
        Vec1d & FR = fR.vec( m );
        Vec1d & ps = psm.vec( m );

        for ( int i = 0; i <= nx; ++ i )
        {
            f[ i ] = 0.5 * ( FR[ i ] + FL[ i ] ) - 0.5 * ps[ i ] * ( qR[ i ] - qL[ i ] );
        }
        for ( int i = 0; i < nx; ++ i )
        {
            res[ i ] -= ( f[ i + 1 ] - f[ i ] ) / dx;
        }
    }

}


void BurgersField::InviscidConservative( VecWrap & u, VecWrap & res )
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

void BurgersField::ViscousResidual( VecWrap & u, VecWrap & res )
{
    ;
}

void BurgersField::Rhs( VecWrap & u, VecWrap & res )
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

void BurgersField::DumpField( Vec1d & x, VecWrap & u )
{
    for ( int i = 0; i < x.size(); ++ i )
    {
        Global::file_string += std::format( "{:.16f}", x[ i ] );
        for ( int m = 0; m < nequ; ++ m )
        {
            Vec1d & u = this->u.vec( m );
            Global::file_string += std::format( " {:.16f}", u[ i ] );
        }
        Global::file_string += std::format( "\n" );
    }
}
