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

#include "UNsLusgs.h"
#include "UNsSpectrum.h"
#include "UCom.h"
#include "UNsCom.h"
#include "Com.h"
#include "UnsGrid.h"
#include "CellMesh.h"
#include "CellTopo.h"
#include "Zone.h"
#include "NsCtrl.h"
#include "NsIdx.h"
#include "HXMath.h"
#include "Parallel.h"
#include "Iteration.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

UNsLusgs::UNsLusgs()
{
}

UNsLusgs::~UNsLusgs()
{
}

void UNsLusgs::SingleSweep()
{
    this->LowerSweep();
    this->UpperSweep();
}

void UNsLusgs::Init()
{
    UnsGrid * grid = Zone::GetUnsGrid();
    FaceTopo * faceTopo = grid->faceTopo;
    CellTopo * cellTopo = grid->cellMesh->cellTopo;
    cellTopo->CalcC2f( faceTopo );
    ug.Init();
    nslu.Init();
    unsf.Init();
    this->CalcSpectrum();
}

void UNsLusgs::CalcSpectrum()
{
    UNsSpectrum * unsSpectrum = new UNsSpectrum();
    unsSpectrum->CalcImplicitSpectrum();
    delete unsSpectrum;
}

void UNsLusgs::LowerSweep()
{
    this->Init();

    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        ug.cId = cId;
        if ( cId == 9 )
        {
            int kkk = 1;
        }

        gcom.blank = ( * ug.blankf )[ ug.cId ];

        if ( this->IsOversetCell() )
        {
            this->ZeroOversetCell();
        }
        else
        {
            this->PrepareSweep();

            this->ZeroFluxIncrement();
        
            this->SolveLowerCell();

            this->CalcLowerChange();
        }

        this->Update();
    }

    //UploadInterfaceValue( grid, dqField, "dqField",  numberOfTotalEquations );
}

void UNsLusgs::UpperSweep()
{
    this->Init();
    //this->LusgsBoundary();
    //DownloadInterfaceValue( grid, dqField, "dqField",  numberOfTotalEquations );
    
    for ( int cId = ug.nCell - 1; cId >= 0; -- cId )
    {
        ug.cId = cId;
        gcom.blank = ( * ug.blankf )[ ug.cId ];

        if ( this->IsOversetCell() )
        {
            this->ZeroOversetCell();
        }
        else
        {
            this->PrepareSweep();

            this->ZeroFluxIncrement();

            this->SolveUpperCell();

            this->CalcUpperChange();
        }

        this->Update();
    }
}

void UNsLusgs::SolveLowerCell()
{
    int fn = ( * ug.c2f )[ ug.cId ].size();
    for ( int iFace = 0; iFace < fn; ++ iFace )
    {
        int fId = ( * ug.c2f )[ ug.cId ][ iFace ];
        this->SolveLower( fId );
    }
}

void UNsLusgs::SolveUpperCell()
{
    int fn = ( * ug.c2f )[ ug.cId ].size();
    for ( int iFace = 0; iFace < fn; ++ iFace )
    {
        int fId = ( * ug.c2f )[ ug.cId ][ iFace ];
        this->SolveUpper( fId );
    }
}

void UNsLusgs::SolveLower( int fId )
{
    if ( this->CanNotLowerSolve( fId ) ) return;

    this->Solve( fId, - 1 );
}

void UNsLusgs::SolveUpper( int fId )
{
    if ( this->CanNotUpperSolve( fId ) ) return;

    this->Solve( fId, - 1 );
}

bool UNsLusgs::CanNotLowerSolve( int fId )
{
    ug.fId = fId;

    ug.lc = ( * ug.lcf )[ ug.fId ];
    ug.rc = ( * ug.rcf )[ ug.fId ];

    // One of lc  and rc must be cell itself.
    // Now its neighboring cell belongs to lower triangular
    return ( ug.lc > ug.cId || ug.rc > ug.cId );
}

bool UNsLusgs::CanNotUpperSolve( int fId )
{
    ug.fId = fId;

    ug.lc = ( * ug.lcf )[ ug.fId ];
    ug.rc = ( * ug.rcf )[ ug.fId ];

    // One of lc  and rc must be cell itself.
    // Now its neighboring cell belongs to upper triangular
    return ( ug.lc < ug.cId || ug.rc < ug.cId );
}

void UNsLusgs::Solve( int fId, int signValue )
{
    ug.fId = fId;

    if ( fId == 147489 )
    {
        int kkk = 1;
    }

    ug.lc = ( * ug.lcf )[ ug.fId ];
    ug.rc = ( * ug.rcf )[ ug.fId ];

    this->SetMeshGeometry();

    this->PrepareData();

    this->GetStandardFluxIncrement( signValue );

    this->CalcViscousTerm();

    this->AddFluxIncrement();
}

void UNsLusgs::SetMeshGeometry()
{
    gcom.SetGeometry();
}

void UNsLusgs::PrepareData()
{
    for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
    {
        nslu.primj[ iEqu ] = ( * unsf.q )[ iEqu ][ ug.rc ]; //qField存的是原始变量！
    }

    for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
    {
        nslu.dqj[ iEqu ] = ( * unsf.dq )[ iEqu ][ ug.rc ];
    }

    this->PrepareDataFacePrim();

    nslu.gama = ( * unsf.gama )[ 0 ][ ug.rc ];
    nscom.visl = ( * unsf.visl )[ 0 ][ ug.rc ];
    nscom.vist = ( * unsf.vist )[ 0 ][ ug.rc ];

    for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
    {
        nscom.q1[ iEqu ] = ( * unsf.q )[ iEqu ][ ug.lc ];
        nscom.q2[ iEqu ] = ( * unsf.q )[ iEqu ][ ug.rc ];
    }
}

void UNsLusgs::PrepareDataFacePrim()
{
    Real rl = ( * unsf.q )[ IDX::IR ][ ug.lc ];
    Real ul = ( * unsf.q )[ IDX::IU ][ ug.lc ];
    Real vl = ( * unsf.q )[ IDX::IV ][ ug.lc ];
    Real wl = ( * unsf.q )[ IDX::IW ][ ug.lc ];
    Real pl = ( * unsf.q )[ IDX::IP ][ ug.lc ];

    Real rr = ( * unsf.q )[ IDX::IR ][ ug.rc ];
    Real ur = ( * unsf.q )[ IDX::IU ][ ug.rc ];
    Real vr = ( * unsf.q )[ IDX::IV ][ ug.rc ];
    Real wr = ( * unsf.q )[ IDX::IW ][ ug.rc ];
    Real pr = ( * unsf.q )[ IDX::IP ][ ug.rc ];

    Real gl = ( * unsf.gama )[ 0 ][ ug.lc ];
    Real gr = ( * unsf.gama )[ 0 ][ ug.rc ];

    Real hl = gl / ( gl - 1.0 ) * pl / rl + half * SQR( ul, vl, wl );
    Real hr = gr / ( gr - 1.0 ) * pr / rr + half * SQR( ur, vr, wr );

    Real tmp0 = sqrt( rr / rl );
    Real tmp1 = 1.0 / ( 1.0 + tmp0 );

    Real rm = sqrt( rr * rl );
    Real um = ( ul + ur * tmp0 ) * tmp1;
    Real vm = ( vl + vr * tmp0 ) * tmp1;
    Real wm = ( wl + wr * tmp0 ) * tmp1;
    Real hm = ( hl + hr * tmp0 ) * tmp1;
    Real pm = rm * ( hm - half * SQR( um, vm, wm ) ) * ( gl - 1.0 ) / gl;

    if ( pm <= 0.0 ) cout << "pm = " << pm << endl;

    nslu.primF[ IDX::IR ] = rm;
    nslu.primF[ IDX::IU ] = um;
    nslu.primF[ IDX::IV ] = vm;
    nslu.primF[ IDX::IW ] = wm;
    nslu.primF[ IDX::IP ] = pm;

    for ( int iEqu = nslu.nBEqu; iEqu < nslu.nEqu; ++ iEqu ) 
    {
        nslu.primF[ iEqu ] = half * ( ( * unsf.q )[ iEqu ][ ug.lc ] + ( * unsf.q )[ iEqu ][ ug.rc ] ); 
    }
}

void UNsLusgs::CalcViscousTerm()
{
    if ( vis_model.vismodel == 0 ) return;

    if ( nscom.visSRModel == 1 )
    {
        Real density = half * ( nscom.q1[ IDX::IR ] + nscom.q2[ IDX::IR ] );

        Real c1 = 4.0 / 3.0 * ( nscom.visl + nscom.vist );
        Real c2 = nscom.gama * ( nscom.visl * nscom.oprl + nscom.vist * nscom.oprt );
        Real c3 = two * MAX( c1, c2 ) / ( nscom.reynolds * density );
        Real farea2 = SQR( gcom.farea );

        nscom.vissr = farea2 * c3;

        nslu.visrad = nscom.vissr / ( * ug.cvol )[ ug.rc ];
    }
    else
    {
        Real dist = ABS(  gcom.xfn * ( gcom.xcc2 - gcom.xcc1 )
                        + gcom.yfn * ( gcom.ycc2 - gcom.ycc1 )
                        + gcom.zfn * ( gcom.zcc2 - gcom.zcc1 ) );

        Real viscosity = nscom.visl + nscom.vist;
        Real density   = half * ( nscom.q1[ IDX::IR ] + nscom.q2[ IDX::IR ] );

        Real c1  = 2.0 * viscosity / ( density * dist * nscom.reynolds + SMALL );
        nscom.vissr = half * c1 * gcom.farea;
        nslu.visrad = nscom.vissr;
    }

    for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
    {
        nslu.rhs0[ iEqu ] -= nslu.visrad * nslu.dqj[ iEqu ];
    }
}

void UNsLusgs::PrepareSweep()
{
    for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
    {
        gcom.blank = ( * ug.blankf )[ ug.cId ];

        nslu.dqi[ iEqu ] = ( * unsf.dq  )[ iEqu ][ ug.cId ]; //dqField的初值为0（守恒或者原始变量）
        nslu.rhs[ iEqu ] = ( * unsf.rhs )[ iEqu ][ ug.cId ]; //RHS还是存在RHS里面比较好
    }

    if ( nslu.numberOfSweeps > 1 )
    {
        for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
        {
            nslu.dqi0[ iEqu ] = nslu.dqi[ iEqu ];
            nslu.drhs[ iEqu ] = ( * unsf.drhs )[ iEqu ][ ug.cId ];

            nslu.dqi[ iEqu ] = ( * unsf.rhs )[ iEqu ][ ug.cId ] - nslu.drhs[ iEqu ];
            ( * unsf.dq )[ iEqu ][ ug.cId ] = nslu.dqi[ iEqu ];
            nslu.drhs[ iEqu ] = 0.0;
        }
    }

    for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
    {
        nslu.radius[ iEqu ] = ( * unsf.impsr )[ 0 ][ ug.cId ];
    }
}

void UNsLusgs::Update()
{
    if ( nslu.numberOfSweeps > 1 )
    {
        for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
        {
            ( * unsf.dq   )[ iEqu ][ ug.cId ]  = nslu.dqi[ iEqu ];
            ( * unsf.drhs )[ iEqu ][ ug.cId ]  = nslu.drhs[ iEqu ];
        }
    }
    else
    {
        for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
        {
            ( * unsf.dq   )[ iEqu ][ ug.cId ]  = nslu.dqi[ iEqu ];
        }
    }
}

EndNameSpace