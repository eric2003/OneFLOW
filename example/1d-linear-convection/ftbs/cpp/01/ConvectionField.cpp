#include "ConvectionField.h"
#include "Weno.h"
#include "hxmath.h"
#include "Grid.h"
#include "Global.h"
#include "cgnslib.h"
#include <iostream>
#include <numbers>
#include <print>

void ConvectionField::InitFieldCommon( Grid * grid )
{
    this->grid = grid;
    this->nequ = Global::nequ;
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
    this->dt = Global::dt;
    this->nt = std::round( Global::total_time / dt );

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

    this->c = 1.0;
}

void ConvectionField::InitFieldAsRestart( Grid * grid )
{
    Vec1d &u = this->u.vec();

    if ( Global::ifinite_volume == 0 )
    {
        //node
        Vec1d & x = grid->x;
        for ( int i = 0; i < ni; ++ i )
        {
            if ( x[ i ] >= 0.5 && x[ i ] <= 1.0 )
            {
                u[ i ] = 2;
            }
            else
            {
                u[ i ] = 1;
            }
        }
    }
    else
    {
        //cell center
        Vec1d & xcc = grid->xcc;

        for ( int i = 0; i < nic; ++ i )
        {
            if ( xcc[ i ] >= 0.5 && xcc[ i ] <= 1.0 )
            {
                u[ i ] = 2;
            }
            else
            {
                u[ i ] = 1;
            }
        }
    }
}

void ConvectionField::ReadFlowField( std::fstream & file, Grid * grid )
{
    if ( Global::ifinite_volume == 1 )
    {
        this->ReadFlowField( file, grid->xcc, u );
    }
    else
    {
        this->ReadFlowField( file, grid->x, u );
    }
}

void ConvectionField::ReadFlowField( std::fstream & file, Vec1d & x, VecWrap & u )
{
    for ( int i = 0; i < x.size(); ++ i )
    {
        std::string line; 
        std::getline( file, line );
        std::stringstream ss( line );
        std::string item;
        std::vector<std::string> row;
        while ( std::getline(ss, item, ' ') )
        {  
            row.push_back( item );
        }
        double um  = std::atof( row[ 1 ].data() );
        this->u[ 0 ][ i ] = um;
    }
    int kkk = 1;
}

void ConvectionField::FTCS( Zone * zone )
{
    this->Rhs( this->u, this->res );
    for ( int m = 0; m < nequ; ++ m )
    {
        Vec1d & u = this->u.vec( m );
        Vec1d & res = this->res.vec( m );
        for ( int i = 0; i < nx; ++ i )
        {
            u[ i ] = u[ i ] + dt * res[ i ];
        }
    }
}

void ConvectionField::CN( Zone * zone )
{
    double rr = 0.5 * this->alpha * dt / ( dx * dx );
    for ( int m = 0; m < nequ; ++ m )
    {
        Vec1d & u = this->u.vec( m );

        std::vector<double> a( ni );
        std::vector<double> b( ni );
        std::vector<double> c( ni );
        std::vector<double> d( ni );

        for ( int i = 0; i < nx; ++ i )
        {
            a[ i ] = - rr;
            b[ i ] = 1.0 + 2.0 * rr;
            c[ i ] = - rr;
        }

        for ( int i = 0; i < nx; ++ i )
        {
            d[ i ] = rr * u[ i - 1 ] + ( 1.0 - 2.0 * rr ) * u[ i ] + rr * u[ i + 1 ];
        }

        double uleft = u[ -1 ] + 2 * rr * ( u[ - 1 ] - 2 * u[ 0 ] + u[ 1 ] );
        double uright = u[ nx ] + 2 * rr * ( u[ nx - 2 ] - 2 * u[ nx - 1 ] + u[ nx ] );

        d[ 0 ] -= a[ 0 ] * uleft;
        d[ nx - 1 ] -= c[ nx - 1 ] * uright;

        std::vector<double> values( d.size() );

        thomas_algorithm( a, b, c, d, values );

        for ( int i = 0; i < nx; ++ i )
        {
            u[ i ] = values[ i ];
        }
    }
}

void ConvectionField::UpdateOldField()
{
    this->un = this->u;
}

void ConvectionField::InviscidResidual( VecWrap & u, VecWrap & res )
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

void ConvectionField::InviscidNonConservative( VecWrap & u, VecWrap & res )
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
    else if ( Global::scheme.reconstruction == to_int( BasicScheme::UpWind1 ) )
    {
        Upwind1L( nx, u, uL );
        Upwind1R( nx, u, uR );
    }

    for ( int m = 0; m < nequ; ++ m )
    {
        Vec1d & u = this->u.vec( m );
        Vec1d & res = this->res.vec( m );
        Vec1d & u_L = uL.vec( m );
        Vec1d & u_R = uR.vec( m );
        double cp = 0.5 * ( c + std::abs( c ) );
        double cn = 0.5 * ( c - std::abs( c ) );
        for ( int i = 0; i < nx; ++ i )
        {
            double fip = cp * u_L[ i + 1 ] + cn * u_R[ i + 1 ];
            double fim = cp * u_L[ i ] + cn * u_R[ i ];
            res[ i ] += - ( fip - fim ) / dx;
        }
    }
}

void ConvectionField::InviscidConservative( VecWrap & u, VecWrap & res )
{
}

void ConvectionField::ViscousResidual( VecWrap & u, VecWrap & res )
{
    double coef = this->alpha / ( dx * dx );
    for ( int m = 0; m < nequ; ++ m )
    {
        Vec1d & u = this->u.vec( m );
        Vec1d & res = this->res.vec( m );
        for ( int i = 0; i < ni; ++ i )
        {
            res[ i ] += coef * ( u[ i + 1 ] - 2.0 * u[ i ] + u[ i - 1 ] );
        }
    }
}
void ConvectionField::Rhs( VecWrap & u, VecWrap & res )
{
    res = 0;
    InviscidResidual( u, res );
    ViscousResidual( u, res );
}

void ConvectionField::DumpField( Grid * grid )
{
    this->DumpField( grid->x, u );
}

void ConvectionField::PostProcess( Grid * grid )
{
    this->DumpField( grid->x, u );
}

void ConvectionField::DumpField( Vec1d & x, VecWrap & u )
{
    for ( int i = 0; i < x.size(); ++ i )
    {
        Global::file_string += std::format( "{:.25f}", x[ i ] );
        for ( int m = 0; m < nequ; ++ m )
        {
            Vec1d & u = this->u.vec( m );
            Global::file_string += std::format( " {:.25f}", u[ i ] );
        }
        Global::file_string += std::format( "\n" );
    }
}