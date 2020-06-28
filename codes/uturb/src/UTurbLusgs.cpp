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

#include "UTurbLusgs.h"
#include "UNsSpectrum.h"
#include "UTurbSpectrum.h"
#include "UTurbCom.h"
#include "UCom.h"
#include "TurbCom.h"
#include "TurbRhs.h"
#include "UNsCom.h"
#include "Com.h"
#include "UnsGrid.h"
#include "CellMesh.h"
#include "CellTopo.h"
#include "Zone.h"
#include "Ctrl.h"
#include "NsIdx.h"
#include "HXMath.h"
#include "Parallel.h"
#include "Iteration.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

UTurbLusgs::UTurbLusgs()
{
}

UTurbLusgs::~UTurbLusgs()
{
}

void UTurbLusgs::SingleSweep()
{
    this->LowerSweep();
    this->UpperSweep();
}

void UTurbLusgs::Init()
{
    UnsGrid * grid = Zone::GetUnsGrid();
    FaceTopo * faceTopo = grid->faceTopo;
    CellTopo * cellTopo = grid->cellMesh->cellTopo;
    cellTopo->CalcC2f( faceTopo );
    ug.Init();
    turblu.Init();
    uturbf.Init();
}

void UTurbLusgs::ReadTmp()
{
    static int iii = 0;
    if ( iii ) return;
    iii = 1;
    fstream file;
    file.open( "turbtmpres.dat", ios_base::in | ios_base::binary );

    for ( int iCell = 0; iCell < ug.nTCell; ++ iCell )
    {
        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            file.read( reinterpret_cast< char * >( & ( * uturbf.res )[ iEqu ][ iCell ] ), sizeof( double ) );
        }
    }
    file.close();
    file.clear();
}

void UTurbLusgs::LowerSweep()
{
    this->Init();

    //ReadTmp();
    
    for ( int cId = 0; cId < ug.nCell; ++ cId )
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
        
            this->SolveLowerCell();

            this->CalcLowerChange();
        }

        this->Update();
    }
    
    //UploadInterfaceValue( grid, dqField, "dqField",  numberOfTotalEquations );
}

void UTurbLusgs::UpperSweep()
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

    //CalcTurbulentViscosity();
}

void UTurbLusgs::SolveLowerCell()
{
    int fn = ( * ug.c2f )[ ug.cId ].size();
    for ( int iFace = 0; iFace < fn; ++ iFace )
    {
        int fId = ( * ug.c2f )[ ug.cId ][ iFace ];
        this->SolveLower( fId );
    }
}

void UTurbLusgs::SolveUpperCell()
{
    int fn = ( * ug.c2f )[ ug.cId ].size();
    for ( int iFace = 0; iFace < fn; ++ iFace )
    {
        int fId = ( * ug.c2f )[ ug.cId ][ iFace ];
        this->SolveUpper( fId );
    }
}

void UTurbLusgs::SolveLower( int fId )
{
    if ( this->CanNotLowerSolve( fId ) ) return;

    this->Solve( fId, - 1 );
}

void UTurbLusgs::SolveUpper( int fId )
{
    if ( this->CanNotUpperSolve( fId ) ) return;

    this->Solve( fId, - 1 );
}

bool UTurbLusgs::CanNotLowerSolve( int fId )
{
    ug.fId = fId;

    ug.lc = ( * ug.lcf )[ ug.fId ];
    ug.rc = ( * ug.rcf )[ ug.fId ];

    // One of lc  and rc must be cell itself.
    // Now its neighboring cell belongs to lower triangular
    return ( ug.lc > ug.cId || ug.rc > ug.cId );
}

bool UTurbLusgs::CanNotUpperSolve( int fId )
{
    ug.fId = fId;

    ug.lc = ( * ug.lcf )[ ug.fId ];
    ug.rc = ( * ug.rcf )[ ug.fId ];

    // One of lc  and rc must be cell itself.
    // Now its neighboring cell belongs to upper triangular
    return ( ug.lc < ug.cId || ug.rc < ug.cId );
}

void UTurbLusgs::Solve( int fId, int signValue )
{
    ug.fId = fId;

    ug.lc = ( * ug.lcf )[ ug.fId ];
    ug.rc = ( * ug.rcf )[ ug.fId ];

    this->SetMeshGeometry();

    this->PrepareData();

    this->GetStandardFluxIncrement( signValue );

    this->AddFluxIncrement();
}

void UTurbLusgs::SetMeshGeometry()
{
    gcom.SetGeometry();
}

void UTurbLusgs::PrepareData()
{
    if ( gcom.swapflag )
    {
        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            turblu.matrix[ iEqu ] = ( * uturbf.matrix_l )[ iEqu ][ ug.fId ];
        }
    }
    else
    {
        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            turblu.matrix[ iEqu ] = ( * uturbf.matrix_r )[ iEqu ][ ug.fId ];
        }
    }

    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turblu.dqj[ iEqu ] = ( * uturbf.dq )[ iEqu ][ ug.rc ];
    }

}

void UTurbLusgs::PrepareSweep()
{
    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        gcom.blank = ( * ug.blankf )[ ug.cId ];

        turblu.dqi[ iEqu ] = ( * uturbf.dq  )[ iEqu ][ ug.cId ]; //dqField的初值为0（守恒或者原始变量）
        turblu.rhs[ iEqu ] = ( * uturbf.rhs )[ iEqu ][ ug.cId ]; //RHS还是存在RHS里面比较好
    }

    if ( turblu.numberOfSweeps > 1 )
    {
        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            turblu.dqi0[ iEqu ] = turblu.dqi[ iEqu ];
            turblu.drhs[ iEqu ] = ( * uturbf.drhs )[ iEqu ][ ug.cId ];

            turblu.dqi[ iEqu ] = ( * uturbf.rhs )[ iEqu ][ ug.cId ] - turblu.drhs[ iEqu ];
            ( * unsf.dq )[ iEqu ][ ug.cId ] = turblu.dqi[ iEqu ];
            turblu.drhs[ iEqu ] = 0.0;
        }
    }

    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turblu.radius[ iEqu ] = ( * uturbf.impsr )[ iEqu ][ ug.cId ];
    }
}

void UTurbLusgs::Update()
{
    if ( turblu.numberOfSweeps > 1 )
    {
        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            ( * uturbf.dq   )[ iEqu ][ ug.cId ]  = turblu.dqi[ iEqu ];
            ( * uturbf.drhs )[ iEqu ][ ug.cId ]  = turblu.drhs[ iEqu ];
        }
    }
    else
    {
        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            ( * uturbf.dq   )[ iEqu ][ ug.cId ]  = turblu.dqi[ iEqu ];
        }
    }
}



EndNameSpace