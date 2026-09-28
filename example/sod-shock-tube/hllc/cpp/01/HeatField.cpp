#include "HeatField.h"
#include "hxmath.h"
#include "Grid.h"
#include "Global.h"
#include "cgnslib.h"
#include <iostream>
#include <numbers>

void HeatField::Init( Grid * grid )
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
    std::cout << "ni = " << ni << "\n";

    Vec1d & x = grid->x;
    this->dx = std::abs( x[ 1 ] - x[ 0 ] );
    this->dt = dx / 10.0;
    this->nt = std::round( Global::total_time / dt );

    std::cout << "this->dt   = " << this->dt << "\n";
    std::cout << "this->nt   = " << this->nt << "\n";
    std::cout << "this->ni   = " << this->ni << "\n";
    std::cout << "nt * dt = " << nt * dt << "\n";

    Global::nt = nt;

    this->alpha = 1 / ( std::numbers::pi * std::numbers::pi );
    this->beta = this->alpha * dt / ( dx * dx );

    int ist = 0 - Global::nghost;
    int ied = this->ni - 1 + Global::nghost;

    this->u.Allocate( this->nequ, ist, ied );
    this->un.Allocate( this->nequ, ist, ied );
    this->res.Allocate( this->nequ, 0, this->nx ); //N+1

    Vec1d &u = this->u.vec();

    if ( Global::ifinite_volume == 0 )
    {
        //node
        for ( int i = 0; i < ni; ++ i )
        {
            u[ i ] = - std::sin( std::numbers::pi * x[ i ] );  //initial condition @ t=0
        }
    }
    else
    {
        //cell center
        Vec1d & xcc = grid->xcc;
        for ( int i = 0; i < nic; ++ i )
        {
            u[ i ] = - std::sin( 2.0 * std::numbers::pi * xcc[ i ] );  //initial condition @ t=0
        }

    }
    int kkk = 1;
}

void HeatField::FTCS( Zone * zone )
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

void HeatField::CN( Zone * zone )
{
    double rr = 0.5 * this->alpha * dt / ( dx * dx );
    for ( int m = 0; m < nequ; ++ m )
    {
        Vec1d & u = this->u.vec( m );

        std::vector<double> a( ni );
        std::vector<double> b( ni );
        std::vector<double> c( ni );
        std::vector<double> d( ni );

        for ( int i = 0; i < ni; ++ i )
        {
            a[ i ] = - rr;
            b[ i ] = 1.0 + 2.0 * rr;
            c[ i ] = - rr;
        }

        for ( int i = 0; i < ni; ++ i )
        {
            d[ i ] = rr * u[ i - 1 ] + ( 1.0 - 2.0 * rr ) * u[ i ] + rr * u[ i + 1 ];
        }

        double uleft = u[ -1 ] + 2 * rr * ( u[ - 1 ] - 2 * u[ 0 ] + u[ 1 ] );
        double uright = u[ ni ] + 2 * rr * ( u[ ni - 2 ] - 2 * u[ ni - 1 ] + u[ ni ] );

        d[ 0 ] -= a[ 0 ] * uleft;
        d[ ni - 1 ] -= c[ ni - 1 ] * uright;

        std::vector<double> values( d.size() );

        thomas_algorithm( a, b, c, d, values );

        for ( int i = 0; i < ni; ++ i )
        {
            u[ i ] = values[ i ];
        }
    }
}

void HeatField::ICP( Zone * zone )
{
    double beta = 0.5 * this->alpha * dt / ( dx * dx );
    for ( int m = 0; m < nequ; ++ m )
    {
        Vec1d & u = this->u.vec( m );
        std::vector<double> a( ni );//0:ni-1
        std::vector<double> b( ni );//0:ni-1
        std::vector<double> c( ni );//0:ni-1
        std::vector<double> d( ni );//0:ni-1

        for ( int i = 0; i < ni; ++ i )
        {
            a[ i ] = 1.0 / 12.0 - beta;
            b[ i ] = 10.0 / 12.0 + 2.0 * beta;
            c[ i ] = 1.0 / 12.0 - beta;

            double aa = 1.0 / 12.0 + beta;
            double bb = 10.0 / 12.0 - 2.0 * beta;
            double cc = 1.0 / 12.0 + beta;
            d[ i ] = aa * u[ i - 1 ] + bb * u[ i ] + cc * u[ i + 1 ];
        }


        //double uleft = u[ -1 ] + 2 * beta * ( u[ - 1 ] - 2 * u[ 0 ] + u[ 1 ] );
        //double uright = u[ ni ] + 2 * beta * ( u[ ni - 2 ] - 2 * u[ ni - 1 ] + u[ ni ] );

        //double uleft = u[ -1 ];
        //double uright = u[ ni ];

        double uleft = u[ -1 ] + beta * ( u[ - 1 ] - 2 * u[ 0 ] + u[ 1 ] );
        double uright = u[ ni ] + beta * ( u[ ni - 2 ] - 2 * u[ ni - 1 ] + u[ ni ] );

        d[ 0 ] -= a[ 0 ] * uleft;
        d[ ni - 1 ] -= c[ ni - 1 ] * uright;

        //d[ 0 ] -= a[ 0 ] * u[ -1 ];
        //d[ ni - 1 ] -= c[ ni - 1 ] * u[ ni ];


        std::vector<double> values( d.size() );

        thomas_algorithm( a, b, c, d, values );

        for ( int i = 0; i < ni; ++ i )
        {
            u[ i ] = values[ i ];
        }
    }
}

void HeatField::UpdateOldField()
{
    this->un = this->u;
}

void HeatField::InviscidResidual( VecWrap & u, VecWrap & res )
{
    ;
}

void HeatField::ViscousResidual( VecWrap & u, VecWrap & res )
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
void HeatField::Rhs( VecWrap & u, VecWrap & res )
{
    res = 0;
    InviscidResidual( u, res );
    ViscousResidual( u, res );
}

void HeatField::DumpField( Grid * grid )
{
    this->DumpField( grid->x, u );
}

void HeatField::PostProcess( Grid * grid )
{
    this->DumpField( grid->x, u );
}

void HeatField::DumpField( Vec1d & x, VecWrap & u )
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