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

void ScalarZone::Init( int ist, int ied )
{
    this->ist = ist;
    this->ied = ied;
    this->nNode = ied - ist + 1;
    this->u.resize( nNode );
    this->un.resize( nNode );
    this->x.resize( nNode );
}

void ScalarZone::SetBc( int bcL, int bcR )
{
    this->bcPointList.push_back( 0 );
    this->bcTypeList.push_back( bcL );

    this->bcPointList.push_back( nNode - 1 );
    this->bcTypeList.push_back( bcR );
}

void ScalarZone::InitField( vector< double > & uGlobal )
{
    for ( int i = 0; i < nNode; ++ i )
    {
        int j = ( this->ist - 1 ) + i;
        this->u[ i ] = uGlobal[ j ];
    }
}

void ScalarZone::Solve( double coef )
{
    int ist = 1;
    int ied = nNode;
    for ( int i = ist; i < ied; ++ i )
    {
        //u[ i ] = un[ i ] - c * dt / dx * ( un[ i ] - un[ i - 1 ] );
        u[ i ] = un[ i ] - coef * ( un[ i ] - un[ i - 1 ] );
    }
}

void ScalarZone::UpdateUN()
{
    Numpy::Copy( u, un );
}

void ScalarZone::GatherField( vector< double > & ugfield )
{
    for ( int i = 0; i < nNode; ++ i )
    {
        int j = ( this->ist - 1 ) + i;
        ugfield[ j ] = this->u[ i ];
    }
}

void ScalarZone::CompareField( vector< double > & uGlobal )
{
    for ( int i = 0; i < nNode; ++ i )
    {
        int j = ( this->ist - 1 ) + i;
        double ug = uGlobal[ j ];
        double um = this->u[ i ];
        double diff = ABS( ug - um );
        if ( diff > 1.0e-12 )
        {
            cout << " i = " << i << "um = " << um << " ig = " << j << " ug = " << ug << "\n";
        }
    }
}

double ScalarZone::GetRightBcValue()
{
    return this->u[ nNode - 1 ];
}

void ScalarZone::SetLeftBcValue( double lv )
{
    this->u[ 0 ] = lv;
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

void ScalarSolver::InitCtrlParameter()
{
    this->nx = 41;
    this->len = 2.0;
    this->dx = len / ( nx - 1.0 );
    //this->nt = 25;
    this->nt = 5;
    //this->dt = 0.025;
    this->dt = 0.05;
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
    this->SetScalarZone();
}

void ScalarSolver::SetScalarZone()
{
    int nZones = 4;
    int dn = nx / nZones;

    vector< int > idxList( nZones + 1 );
    idxList[ 0 ] = 1;
    idxList[ nZones ] = nx;
    for ( int iZone = 1; iZone < nZones; ++ iZone )
    {
        idxList[ iZone ] = iZone * dn + 1;
    }

    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        ScalarZone * scalarZone = new ScalarZone();
        int ist = idxList[ iZone ];
        int ied = idxList[ iZone + 1 ];

        scalarZone->zoneid = iZone;
        scalarZone->Init( ist, ied );
        int bcL = SCALAR_INTERFACE;
        int bcR = SCALAR_COMPUTE;
        if ( iZone == 0 ) bcL = SCALAR_EXTRAPOLATE;
        scalarZone->SetBc( bcL, bcR );
        
        this->scalarZones.push_back( scalarZone );
        cout << " iZone = " << iZone << " nZones = " << nZones << " ist = " << ist << " ied = " << ied << "\n";
    }
}

void ScalarSolver::InitFlowField()
{
    u.resize( nx );
    Numpy::Ones( u );

    for ( int i = 0; i < nx; ++ i )
    {
        double xm = this->x[ i ];
        u[ i ] = this->ScalarFun( xm );
    }

    un.resize( nx );
    Numpy::Ones( un );

    this->InitZoneFlowField();

}

void ScalarSolver::InitZoneFlowField()
{
    int nZones = this->scalarZones.size();
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        ScalarZone * scalarZone = this->scalarZones[ iZone ];
        scalarZone->InitField( this->u );
    }
}

void ScalarSolver::Run()
{
    this->Init();

    this->SolveFlowField();

    this->Visual();
   
}

void ScalarSolver::RunTest( ScalarPara * para )
{
    this->InitCtrlParameterTest( para );

    this->InitGrid();

    this->InitFlowField();

    this->SolveFlowField();

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


void ScalarSolver::SolvePart( ScalarZone * scalarZone )
{
    int ist = scalarZone->ist;
    int ied = scalarZone->ied;
    for ( int i = ist; i < ied; ++ i )
    {
        u[ i ] = un[ i ] - c * dt / dx * ( un[ i ] - un[ i - 1 ] );
    }

    double coef = c * dt / dx;
    scalarZone->Solve( coef );
}

void ScalarSolver::SolveOneStep()
{
    int nZones = this->scalarZones.size();
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        this->SolvePart( this->scalarZones[ iZone ] );
    }
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
    vector< double > ugfield( u.size() );

    int nZones = this->scalarZones.size();

    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        ScalarZone * scalarZone = this->scalarZones[ iZone ];
        scalarZone->GatherField( ugfield );
    }

    Numpy::ToTecplot( "theory_gather_tec.plt", x, ugfield, utheory );
    Numpy::ToTecplot( "cmp_gather_tec.plt", x, u, ugfield );

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
    this->CalcL1Norm();
    this->CalcL2Norm();
    para->x = x;
    para->du = dua;
    para->l1Norm = this->l1Norm;
    para->l2Norm = this->l2Norm;
}

void ScalarSolver::Boundary()
{
    int nZones = this->scalarZones.size();
    vector< double > bclist;
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        ScalarZone * scalarZone = this->scalarZones[ iZone ];
        double rv = scalarZone->GetRightBcValue();
        bclist.push_back( rv );
    }

    for ( int iZone = 1; iZone < nZones; ++ iZone )
    {
        ScalarZone * scalarZone = this->scalarZones[ iZone ];
        double rv = scalarZone->GetRightBcValue();
        int i = iZone - 1;
        double bcv = bclist[ i ];
        scalarZone->SetLeftBcValue( bcv );
    }
}

void ScalarSolver::SolveFlowField()
{
    for ( int n = 0; n < nt; ++ n )
    {
        this->UpdateUN();

        this->SolveOneStep();

        this->Boundary();

        this->CompareField();
    }
}

void ScalarSolver::CompareField()
{
    int nZones = this->scalarZones.size();
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        ScalarZone * scalarZone = this->scalarZones[ iZone ];
        scalarZone->CompareField( this->u );
    }
}

void ScalarSolver::UpdateUN()
{
    Numpy::Copy( u, un );

    int nZones = this->scalarZones.size();
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        ScalarZone * scalarZone = this->scalarZones[ iZone ];
        scalarZone->UpdateUN();
    }
}


double ScalarSolver::ScalarFun( double xm )
{
    return SquareFun( xm );
    //return MyCurve( xm );
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

void ScalarSolver::CalcL1Norm()
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

void ScalarSolver::CalcL2Norm()
{
    this->l2Norm = 0.0;

    for ( int i = 0; i < nx; ++ i )
    {
        double um = utheory[ i ];
        double ui = u[ i ];
        double du = ui - um;
        this->l2Norm += SQR( du ) * this->dx;
    }
    this->l2Norm = sqrt( this->l2Norm );
    cout << " this->l2Norm = " << this->l2Norm << "\n";
}

EndNameSpace