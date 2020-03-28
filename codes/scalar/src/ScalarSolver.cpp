/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2020 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "ScalarSolver.h"
#include "ScalarOrder.h"
#include "Numpy.h"
#include "HXMath.h"
#include <iostream>
#include <vector>
using namespace std;

BeginNameSpace( ONEFLOW )

ScalarZone::ScalarZone()
{
}

ScalarZone::~ScalarZone()
{
}

ScalarSolver::ScalarSolver()
{
}

ScalarSolver::~ScalarSolver()
{
    this->FreeScalarZones();
}

void ScalarSolver::FreeScalarZones()
{
    int nSize = this->scalarZones.size();
    for ( int i = 0; i < nSize; ++ i )
    {
        delete this->scalarZones[ i ];
    }
    this->scalarZones.resize( 0 );
}

void ScalarSolver::Init()
{
    this->InitCtrlParameter();

    this->InitGrid();

    this->InitFlowField();
}

//void ScalarSolver::InitCtrlParameter()
//{
//    this->nx = 41;
//    this->len = 2.0;
//    this->dx = len / ( nx - 1.0 );
//    this->nt = 25;
//    this->dt = 0.025;
//    this->c  = 1;
//}

void ScalarSolver::InitCtrlParameter()
{
    this->nx = 41;
    this->len = 2.0;
    this->dx = len / ( nx - 1.0 );
    //this->nt = 25;
    this->nt = 5;
    this->dt = 0.025;
    this->c  = 1;
    double timeN = 0.625;
}

void ScalarSolver::InitCtrlParameterTest( ScalarPara * para )
{
    this->nx = para->nx;
    this->len = para->len;
    this->dx = para->dx;
    this->nt = para->nt;
    this->dt = para->dt;
    this->c  = para->c;
}

void ScalarSolver::InitGrid()
{
    x.resize( nx );
    Numpy::Linspace( x, 0, len );
}

void ScalarSolver::InitFlowField()
{
    u.resize( nx );
    Numpy::Ones( u );

    //int st = 0.5 / dx;
    //int ed = 1 / dx + 1;

    //Numpy::Set( u, st, ed, 2 );

    for ( int i = 0; i < nx; ++ i )
    {
        double xm = this->x[ i ];
        u[ i ] = this->ScalarFun( xm );
    }

    un.resize( nx );
    Numpy::Ones( un );

    u0.resize( nx );
    Numpy::Copy( u, u0 ); 
}

void ScalarSolver::Run()
{
    this->Init();

    this->SolvePartNew();

    this->Visual();
   
}

void ScalarSolver::RunTest( ScalarPara * para )
{
    this->InitCtrlParameterTest( para );

    this->InitGrid();

    this->InitFlowField();

    this->SolvePartNew();

    this->Visual( para );

    this->FreeScalarZones();
}

void ScalarSolver::Solve()
{
    for ( int n = 0; n < nt; ++ n ) //loop for values of n from 0 to nt, so it will run nt times
    {
        Numpy::Copy( u, un ); //copy the existing values of u into un
        for ( int i = 1; i < nx; ++ i ) //you can try commenting this line and...
        {
            u[ i ] = un[ i ] - c * dt / dx * ( un[ i ] - un[ i - 1 ] );
        }
    }
    Numpy::Plot( x, u );
}

void ScalarSolver::SolvePart()
{
    for ( int n = 0; n < nt; ++ n ) //loop for values of n from 0 to nt, so it will run nt times
    {
        Numpy::Copy( u, un ); //copy the existing values of u into un
        int nx1 = 21;
        for ( int i = 1; i < nx1; ++ i ) //you can try commenting this line and...
        {
            u[ i ] = un[ i ] - c * dt / dx * ( un[ i ] - un[ i - 1 ] );
        }

        for ( int i = nx1; i < nx; ++ i ) //you can try commenting this line and...
        {
            u[ i ] = un[ i ] - c * dt / dx * ( un[ i ] - un[ i - 1 ] );
        }
    }
    Numpy::Plot( x, u );
}

void ScalarSolver::SolvePart( int ist, int ied )
{
    for ( int i = ist; i < ied; ++ i )
    {
        u[ i ] = un[ i ] - c * dt / dx * ( un[ i ] - un[ i - 1 ] );
    }
}

void ScalarSolver::SetScalarZone()
{
    int ist = 1;
    int ied = nx/2;

    int jst = nx/2;
    int jed = nx;

    ScalarZone * scalarZone = 0;
    scalarZone = new ScalarZone();
    scalarZone->ist = ist;
    scalarZone->ied = ied;
    this->scalarZones.push_back( scalarZone );

    scalarZone = new ScalarZone();
    scalarZone->ist = jst;
    scalarZone->ied = jed;
    this->scalarZones.push_back( scalarZone );
}

void ScalarSolver::SolveOneStep()
{
    this->SolvePart( this->scalarZones[ 0 ]->ist, this->scalarZones[ 0 ]->ied );
    this->SolvePart( this->scalarZones[ 1 ]->ist, this->scalarZones[ 1 ]->ied );
}

void ScalarSolver::Visual()
{
    Numpy::Plot( x, u );

    double time_final = this->dt * nt;
    this->Theory( time_final );

    Numpy::Plot( "theory.plt", x, utheory );
    Numpy::ToTecplot( "theory_tec.plt", x, u, utheory );
    Numpy::ToTecplot( "du.plt", x, u, du );
    Numpy::ToTecplot( "dua.plt", x, utheory, dua );
}

void ScalarSolver::Visual( ScalarPara * para )
{
    Numpy::Plot( x, u );

    double time_final = this->dt * nt;
    cout << " nt = " << nt << "\n";
    this->Theory( time_final );

    Numpy::Plot( "theory.plt", x, utheory );
    Numpy::ToTecplot( "theory_tec.plt", x, u, utheory );
    Numpy::ToTecplot( "du.plt", x, u, du );
    Numpy::ToTecplot( "dua.plt", x, utheory, dua );
    this->ComputeL1Norm();
    this->ComputeL2Norm();
    para->x = x;
    para->du = dua;
    para->l1Norm = this->l1Norm;
    para->l2Norm = this->l2Norm;
}

void ScalarSolver::SolvePartNew()
{
    this->SetScalarZone();

    for ( int n = 0; n < nt; ++ n )
    {
        Numpy::Copy( u, un );

        this->SolveOneStep();
    }
}

double ScalarSolver::ScalarFun( double xm )
{
    //return SquareFun( xm );
    return MyCurve( xm );
    //return SinFun( xm );
}

double ScalarSolver::SquareFun( double xm )
{
    if ( xm >= 0.5 && xm <= 1.0 )
    {
        return 2.0;
    }
    return 1.0;
}

double ScalarSolver::SinFun( double xm )
{
    if ( xm >= 0.5 && xm <= 1.0 )
    {
        double cccp = 3.14150965;
        double dx = ( xm - 0.5 ) / 0.5 * cccp;
        double ff = 1.0 + sin( dx );
        return ff;
    }
    return 1.0;
}


double ScalarSolver::MyCurve( double xm )
{
    if ( xm < 0.5 ) return 1.0;
    if ( xm > 1.0 ) return 2.0;
    //if ( xm >= 0.5 && xm <= 1.0 )
    double dx = xm - 0.5;
    double ff = 1.0 + dx * 2;
    return ff;
}

void ScalarSolver::Theory( double t )
{
    this->utheory.resize( nx );

    for ( int i = 0; i < nx; ++ i )
    {
        double xm = this->x[ i ];
        double xm_new = xm - c * t;
        utheory[ i ] = this->ScalarFun( xm_new );
    }

    this->du.resize( nx );
    this->dua.resize( nx );

    for ( int i = 0; i < nx; ++ i )
    {
        du[ i ] = utheory[ i ] - u[ i ];
        dua[ i ] = ABS( du[ i ] );
    }

}

void ScalarSolver::ComputeL1Norm()
{
    this->l1Norm = 0.0;

    for ( int i = 0; i < nx; ++ i )
    {
        double um = utheory[ i ];
        double ui = u[ i ];
        double du = ui - um;
        this->l1Norm += ABS( du ) * this->dx;
    }
    cout << " this->l1Norm = " << this->l1Norm << "\n";
}

//void ScalarSolver::ComputeL2Norm()
//{
//    this->l2Norm = 0.0;
//
//    for ( int i = 0; i < nx; ++ i )
//    {
//        double um = utheory[ i ];
//        double ui = u[ i ];
//        double du = ui - um;
//        this->l2Norm += SQR( du ) * this->dx;
//    }
//    this->l2Norm = sqrt( this->l2Norm );
//    cout << " this->l2Norm = " << this->l2Norm << "\n";
//}

void ScalarSolver::ComputeL2Norm()
{
    this->l2Norm = 0.0;

    for ( int i = 0; i < nx; ++ i )
    {
        double um = utheory[ i ];
        double ui = u[ i ];
        double du = ui - um;
        this->l2Norm += SQR( du );
    }
    this->l2Norm = sqrt( this->l2Norm / nx );
    cout << " this->l2Norm = " << this->l2Norm << "\n";
}


EndNameSpace