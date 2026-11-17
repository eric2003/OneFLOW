#include "HeatField.h"
#include "hxmath.h"
#include "Grid.h"
#include "Global.h"
#include "cgnslib.h"
#include <iostream>
#include <numbers>

void HeatField::Init( Grid * grid )
{
    this->ni = grid->x.size();
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
    u.Allocate( ist, ied, 0 );
    un.Allocate( ist, ied, 0 );
    res.Allocate( 0, this->ni, 0 ); //N+1

    for ( int i = 0; i < ni; ++ i )
    {
        u[ i ] = - std::sin( std::numbers::pi * x[ i ] ); //initial condition @ t=0
    }
    int kkk = 1;
}

void HeatField::FTCS( Zone * zone )
{
    this->Rhs( this->u, this->res );

    for ( int i = 0; i < ni; ++ i )
    {
        u[ i ] = u[ i ] + dt * res[ i ];
    }
}

void HeatField::CN( Zone * zone )
{
    double rr = 0.5 * this->alpha * dt / ( dx * dx );
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

//void HeatField::ICP( Zone * zone )
//{
//    double rr = 0.5 * this->alpha * dt / ( dx * dx );
//    std::vector<double> a( ni );//0:ni-1
//    std::vector<double> b( ni );//0:ni-1
//    std::vector<double> c( ni );//0:ni-1
//    std::vector<double> d( ni );//0:ni-1
//
//    for ( int i = 0; i < ni; ++ i )
//    {
//        a[ i ] = 1.0 / 12.0 - rr;
//        b[ i ] = 10.0 / 12.0 + 2.0 * rr;
//        c[ i ] = 1.0 / 12.0 - rr;
//    }
//
//    a[ 0 ] = 0;
//    b[ 0 ] = 1;
//    c[ 0 ] = 0;
//
//    a[ ni - 1 ] = 0;
//    b[ ni - 1 ] = 1;
//    c[ ni - 1 ] = 0;
//
//    for ( int i = 0; i < ni; ++ i )
//    {
//        double aa = 1.0 / 12.0 + rr;
//        double bb = 10.0 / 12.0 - 2.0 * rr;
//        double cc = 1.0 / 12.0 + rr;
//        d[ i ] = aa * u[ i - 1 ] + bb * u[ i ] + cc * u[ i + 1 ];
//    }
//
//    //d[ 0 ] = 0;
//    //d[ ni - 1 ] = 0;
//
//    std::vector<double> values( d.size() );
//
//    thomas_algorithm( a, b, c, d, values );
//
//    for ( int i = 0; i < ni; ++ i )
//    {
//        u[ i ] = values[ i ];
//    }
//}


void HeatField::ICP( Zone * zone )
{
    double beta = 0.5 * this->alpha * dt / ( dx * dx );
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

void HeatField::UpdateOldField()
{
    this->un = this->u;
}

void HeatField::InviscidResidual( Vec1d & u, Vec1d & res )
{
    ;
}

void HeatField::ViscousResidual( Vec1d & u, Vec1d & res )
{
    double coef = this->alpha / ( dx * dx );

    for ( int i = 0; i < ni; ++ i )
    {
        res[ i ] += coef * ( u[ i + 1 ] - 2.0 * u[ i ] + u[ i - 1 ] );
    }
}
void HeatField::Rhs( Vec1d & u, Vec1d & res )
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

void HeatField::DumpField( Vec1d & x, Vec1d & u )
{
    for ( int i = 0; i < x.size(); ++ i )
    {
        Global::file_string += std::format( "{:.16f} {:.16f}\n", x[ i ], u[ i ] );
    }
}
