/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2024 He Xin and the OneFLOW contributors.
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

#include "UINsLusgs.h"
#include "NsLusgs.h"
#include "UINsSpectrum.h"
#include "UCom.h"
#include "UINsCom.h"
#include "Com.h"
#include "UnsGrid.h"
#include "CellMesh.h"
#include "CellTopo.h"
#include "Zone.h"
#include "INsCtrl.h"
#include "INsIdx.h"
#include "HXMath.h"
#include "Parallel.h"
#include "Iteration.h"
#include <iostream>


BeginNameSpace( ONEFLOW )

UINsLusgs::UINsLusgs()
{
}

UINsLusgs::~UINsLusgs()
{
}

void UINsLusgs::SingleSweep()
{
    this->LowerSweep();
    this->UpperSweep();
}

void UINsLusgs::Init()
{
    UnsGrid * grid = Zone::GetUnsGrid();
    FaceTopo * faceTopo = grid->faceTopo;
    CellTopo * cellTopo = grid->cellMesh->cellTopo;
    cellTopo->CalcC2f( faceTopo );
    ug.Init();
    nslu.Init();
    uinsf.Init();
    this->CalcSpectrum();
}

void UINsLusgs::CalcSpectrum()
{
    UINsSpectrum * unsSpectrum = new UINsSpectrum();
    unsSpectrum->CalcImplicitSpectrum();
    delete unsSpectrum;
}

void UINsLusgs::LowerSweep()
{
    this->Init();

    for ( int cId = 0; cId < ug.nCells; ++ cId )
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

void UINsLusgs::UpperSweep()
{
    this->Init();
    //this->LusgsBoundary();
    //DownloadInterfaceValue( grid, dqField, "dqField",  numberOfTotalEquations );
    
    for ( int cId = ug.nCells - 1; cId >= 0; -- cId )
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

void UINsLusgs::SolveLowerCell()
{
    int fn = ( * ug.c2f )[ ug.cId ].size();
    for ( int iFace = 0; iFace < fn; ++ iFace )
    {
        int fId = ( * ug.c2f )[ ug.cId ][ iFace ];
        this->SolveLower( fId );
    }
}

void UINsLusgs::SolveUpperCell()
{
    int fn = ( * ug.c2f )[ ug.cId ].size();
    for ( int iFace = 0; iFace < fn; ++ iFace )
    {
        int fId = ( * ug.c2f )[ ug.cId ][ iFace ];
        this->SolveUpper( fId );
    }
}

void UINsLusgs::SolveLower( int fId )
{
    if ( this->CanNotLowerSolve( fId ) ) return;

    this->Solve( fId, - 1 );
}

void UINsLusgs::SolveUpper( int fId )
{
    if ( this->CanNotUpperSolve( fId ) ) return;

    this->Solve( fId, - 1 );
}

bool UINsLusgs::CanNotLowerSolve( int fId )
{
    ug.fId = fId;

    ug.lc = ( * ug.lcf )[ ug.fId ];
    ug.rc = ( * ug.rcf )[ ug.fId ];

    // One of lc  and rc must be cell itself.
    // Now its neighboring cell belongs to lower triangular
    return ( ug.lc > ug.cId || ug.rc > ug.cId );
}

bool UINsLusgs::CanNotUpperSolve( int fId )
{
    ug.fId = fId;

    ug.lc = ( * ug.lcf )[ ug.fId ];
    ug.rc = ( * ug.rcf )[ ug.fId ];

    // One of lc  and rc must be cell itself.
    // Now its neighboring cell belongs to upper triangular
    return ( ug.lc < ug.cId || ug.rc < ug.cId );
}

void UINsLusgs::Solve( int fId, int signValue )
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

    this->ComputeViscousTerm();

    this->AddFluxIncrement();
}

void UINsLusgs::SetMeshGeometry()
{
    gcom.SetGeometry();
}

void UINsLusgs::PrepareData()
{
    for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
    {
        nslu.primj[ iEqu ] = ( * uinsf.q )[ iEqu ][ ug.rc ]; //Qfield is the original variable!
    }

    for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
    {
        nslu.dqj[ iEqu ] = ( * uinsf.dq )[ iEqu ][ ug.rc ];
    }

    this->PrepareDataFacePrim();

    nslu.gama = ( * uinsf.gama )[ 0 ][ ug.rc ];
    nscom.visl = ( * uinsf.visl )[ 0 ][ ug.rc ];
    nscom.vist = ( * uinsf.vist )[ 0 ][ ug.rc ];

    for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
    {
        nscom.q1[ iEqu ] = ( * uinsf.q )[ iEqu ][ ug.lc ];
        nscom.q2[ iEqu ] = ( * uinsf.q )[ iEqu ][ ug.rc ];
    }
}

void UINsLusgs::PrepareDataFacePrim()
{
    Real rl = ( * uinsf.q )[ IIDX::IIR ][ ug.lc ];
    Real ul = ( * uinsf.q )[ IIDX::IIU ][ ug.lc ];
    Real vl = ( * uinsf.q )[ IIDX::IIV ][ ug.lc ];
    Real wl = ( * uinsf.q )[ IIDX::IIW ][ ug.lc ];
    Real pl = ( * uinsf.q )[ IIDX::IIP ][ ug.lc ];

    Real rr = ( * uinsf.q )[ IIDX::IIR ][ ug.rc ];
    Real ur = ( * uinsf.q )[ IIDX::IIU ][ ug.rc ];
    Real vr = ( * uinsf.q )[ IIDX::IIV ][ ug.rc ];
    Real wr = ( * uinsf.q )[ IIDX::IIW ][ ug.rc ];
    Real pr = ( * uinsf.q )[ IIDX::IIP ][ ug.rc ];

    Real gl = ( * uinsf.gama )[ 0 ][ ug.lc ];
    Real gr = ( * uinsf.gama )[ 0 ][ ug.rc ];

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

    if ( pm <= 0.0 ) std::cout << "pm = " << pm << std::endl;

    nslu.primF[ IIDX::IIR ] = rm;
    nslu.primF[ IIDX::IIU ] = um;
    nslu.primF[ IIDX::IIV ] = vm;
    nslu.primF[ IIDX::IIW ] = wm;
    nslu.primF[ IIDX::IIP ] = pm;

    for ( int iEqu = nslu.nBEqu; iEqu < nslu.nEqu; ++ iEqu ) 
    {
        nslu.primF[ iEqu ] = half * ( ( * uinsf.q )[ iEqu ][ ug.lc ] + ( * uinsf.q )[ iEqu ][ ug.rc ] ); 
    }
}

void UINsLusgs::ComputeViscousTerm()
{
    if ( vis_model.vismodel == 0 ) return;

    if ( nscom.visSRModel == 1 )
    {
        Real density = half * ( nscom.q1[ IIDX::IIR ] + nscom.q2[ IIDX::IIR ] );

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
        Real density   = half * ( nscom.q1[ IIDX::IIR ] + nscom.q2[ IIDX::IIR ] );

        Real c1  = 2.0 * viscosity / ( density * dist * nscom.reynolds + SMALL );
        nscom.vissr = half * c1 * gcom.farea;
        nslu.visrad = nscom.vissr;
    }

    for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
    {
        nslu.rhs0[ iEqu ] -= nslu.visrad * nslu.dqj[ iEqu ];
    }
}

void UINsLusgs::PrepareSweep()
{
    for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
    {
        gcom.blank = ( * ug.blankf )[ ug.cId ];

        nslu.dqi[ iEqu ] = ( * uinsf.dq  )[ iEqu ][ ug.cId ]; //The initial value of dqfield is 0 (conserved or original variable)
        nslu.rhs[ iEqu ] = ( * uinsf.rhs )[ iEqu ][ ug.cId ]; //It is better to have RHS in RHS
    }

    if ( nslu.numberOfSweeps > 1 )
    {
        for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
        {
            nslu.dqi0[ iEqu ] = nslu.dqi[ iEqu ];
            nslu.drhs[ iEqu ] = ( * uinsf.drhs )[ iEqu ][ ug.cId ];

            nslu.dqi[ iEqu ] = ( * uinsf.rhs )[ iEqu ][ ug.cId ] - nslu.drhs[ iEqu ];
            ( * uinsf.dq )[ iEqu ][ ug.cId ] = nslu.dqi[ iEqu ];
            nslu.drhs[ iEqu ] = 0.0;
        }
    }

    for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
    {
        nslu.radius[ iEqu ] = ( * uinsf.impsr )[ 0 ][ ug.cId ];
    }
}

void UINsLusgs::Update()
{
    if ( nslu.numberOfSweeps > 1 )
    {
        for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
        {
            ( * uinsf.dq   )[ iEqu ][ ug.cId ]  = nslu.dqi[ iEqu ];
            ( * uinsf.drhs )[ iEqu ][ ug.cId ]  = nslu.drhs[ iEqu ];
        }
    }
    else
    {
        for ( int iEqu = 0; iEqu < nslu.nEqu; ++ iEqu )
        {
            ( * uinsf.dq   )[ iEqu ][ ug.cId ]  = nslu.dqi[ iEqu ];
        }
    }
}

EndNameSpace
