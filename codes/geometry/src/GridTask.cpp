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

#include "GridTask.h"
#include "CmxTask.h"
#include "TaskCom.h"
#include "TaskState.h"
#include "Ctrl.h"
#include "SolverDef.h"
#include "SolverState.h"
#include "Parallel.h"
#include "Zone.h"
#include "ZoneState.h"
#include "UnsGrid.h"
#include "CellMesh.h"
#include "ActionState.h"
#include "DataBook.h"
#include "DataBaseIO.h"
#include "InterFace.h"
#include "WallDist.h"
#include "FaceTopo.h"
#include "SolverRegister.h"

BeginNameSpace( ONEFLOW )

void SetGridFunc()
{
    REGISTER_DATA_CLASS( AllocWallDist );
    REGISTER_DATA_CLASS( ReadWallDist  );
    REGISTER_DATA_CLASS( DumpWallDist  );
    REGISTER_DATA_CLASS( CreateCalcMetricsTask  );
    REGISTER_DATA_CLASS( CalcMetrics  );
    REGISTER_DATA_CLASS( SwapCellCenter  );
    REGISTER_DATA_CLASS( DecodeCellCenter );
    SetWallTask();
}

void AllocWallDist( StringField & data )
{
    UnsGrid * grid = Zone::GetUnsGrid();
    grid->cellMesh->AllocDist();
}

void ReadWallDist( StringField & data )
{
    UnsGrid * grid = Zone::GetUnsGrid();
    grid->cellMesh->ReadDist();
}

void DumpWallDist( StringField & data )
{
    UnsGrid * grid = Zone::GetUnsGrid();
    grid->cellMesh->DumpDist();
}

void CreateCalcMetricsTask( StringField & data )
{
    CalcMetricsTask * task = new CalcMetricsTask();
    TaskState::task = task;
}

void CalcMetrics( StringField & data )
{
    UnsGrid * grid = Zone::GetUnsGrid();
    grid->CalcMetrics();
}

void CalcMetricsTask::Run()
{
    ActionState::dataBook = this->dataBook;
    for ( int zId = 0; zId < ZoneState::nZones; ++ zId )
    {
        ZoneState::zid = zId;

        Client2Server( this, WriteScreen );
    }
}

void SwapCellCenter( StringField & data )
{
    UnsGrid * grid = Zone::GetUnsGrid();
    InterFace * interFace = grid->interFace;
    if ( ! ONEFLOW::IsValid( interFace ) ) return;

    int iNei = interFace->z2n[ ZoneState::rzid ];
    int nIFace  = interFace->interFacePairs[ iNei ]->nIFace;

    IntField & interfaceId = interFace->GetInterfaceId( iNei, GREAT_SEND );

    ActionState::dataBook->MoveToBegin();

    CellMesh * cellMesh = grid->cellMesh;

    RealField & xcc = cellMesh->xcc;
    RealField & ycc = cellMesh->ycc;
    RealField & zcc = cellMesh->zcc;
    RealField & vol = cellMesh->vol;

    for ( int iLocalFace = 0; iLocalFace < nIFace; ++ iLocalFace )
    {
        int s1;
        int iFace = interfaceId[ iLocalFace ];

        grid->faceTopo->GetSId( iFace, 1, s1 );

        HXWrite( ActionState::dataBook, xcc[ s1 ] );
        HXWrite( ActionState::dataBook, ycc[ s1 ] );
        HXWrite( ActionState::dataBook, zcc[ s1 ] );
        HXWrite( ActionState::dataBook, vol[ s1 ] );
    }
}

void DecodeCellCenter( StringField & data )
{
    UnsGrid * grid = Zone::GetUnsGrid();
    InterFace * interFace = grid->interFace;
    if ( ! ONEFLOW::IsValid( interFace ) ) return;

    int iNei = interFace->z2n[ ZoneState::szid ];
    int nIFace  = interFace->interFacePairs[ iNei ]->nIFace;
    IntField & interfaceId = interFace->GetInterfaceId( iNei, GREAT_RECV );

    ActionState::dataBook->MoveToBegin();

    CellMesh * cellMesh = grid->cellMesh;

    RealField & xcc = cellMesh->xcc;
    RealField & ycc = cellMesh->ycc;
    RealField & zcc = cellMesh->zcc;
    RealField & vol = cellMesh->vol;

    for ( int iLocalFace = 0; iLocalFace < nIFace; ++ iLocalFace )
    {
        int iFace = interfaceId[ iLocalFace ];
        int t1;
        grid->faceTopo->GetTId( iFace, 1, t1 );

        HXRead( ActionState::dataBook, xcc[ t1 ] );
        HXRead( ActionState::dataBook, ycc[ t1 ] );
        HXRead( ActionState::dataBook, zcc[ t1 ] );
        HXRead( ActionState::dataBook, vol[ t1 ] );
    }
}

SolverRegData gridReg;

SolverRegData * GetGridReg()
{
    gridReg.sTid = GRID_SOLVER;
    gridReg.func = & SetGridFunc;
    gridReg.solverName = "grid";
    gridReg.baseKind = GRID_BASED;
    gridReg.dataFlag = NO_DATA;

    return & gridReg;
}

REGISTER_REG_DATA( GetGridReg );

EndNameSpace