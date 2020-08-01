/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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

#include "Restart.h"
#include "NsRestart.h"
#include "TurbRestart.h"
#include "ActionState.h"
#include "SolverDef.h"
#include "Zone.h"
#include "Grid.h"
#include "Ctrl.h"
#include "InterFace.h"
#include "DataBook.h"
#include "DataBase.h"
#include "DataBaseIO.h"
#include "Iteration.h"
#include "FieldImp.h"
#include "FieldWrap.h"
#include "FieldAlloc.h"
#include "UsdPara.h"
#include "RegisterUtil.h"
#include "INsRestart.h"

BeginNameSpace( ONEFLOW )

Restart * CreateRestart( int sTid )
{
    if ( sTid == NS_SOLVER )
    {
        return CreateNsRestart();
    }
    else if ( sTid == INC_NS_SOLVER )
    {
        return CreateINsRestart();
    }
    else if ( sTid == TURB_SOLVER )
    {
        return CreateTurbRestart();
    }

    return 0;
}

Restart::Restart()
{
    ;
}

Restart::~Restart()
{
    ;
}

void Restart::ReadUnsteady( int sTid )
{
    FieldManager * fieldManager = FieldFactory::GetFieldManager( sTid );

    UsdPara * usdPara = fieldManager->usdPara;

    Grid * grid = Zone::GetGrid();

    MRField * q  = GetFieldPointer< MRField > ( grid, usdPara->flow[ 0 ] );
    MRField * q1 = GetFieldPointer< MRField > ( grid, usdPara->flow[ 1 ] );
    MRField * q2 = GetFieldPointer< MRField > ( grid, usdPara->flow[ 2 ] );

    HXRead( ActionState::dataBook, q1 );
    HXRead( ActionState::dataBook, q2 );
    SetField( q, q1 );

    MRField * res  = GetFieldPointer< MRField > ( grid, usdPara->residual[ 0 ] );
    MRField * res1 = GetFieldPointer< MRField > ( grid, usdPara->residual[ 1 ] );
    MRField * res2 = GetFieldPointer< MRField > ( grid, usdPara->residual[ 2 ] );

    HXRead( ActionState::dataBook, res1 );
    HXRead( ActionState::dataBook, res2 );
    SetField( res, res1 );
}

void Restart::DumpUnsteady( int sTid )
{
    FieldManager * fieldManager = FieldFactory::GetFieldManager( sTid );

    UsdPara * usdPara = fieldManager->usdPara;

    Grid * grid = Zone::GetGrid();

    MRField * q  = GetFieldPointer< MRField > ( grid, usdPara->flow[ 0 ] );
    MRField * q1 = GetFieldPointer< MRField > ( grid, usdPara->flow[ 1 ] );
    MRField * q2 = GetFieldPointer< MRField > ( grid, usdPara->flow[ 2 ] );

    HXWrite( ActionState::dataBook, q1 );
    HXWrite( ActionState::dataBook, q2 );

    MRField * res  = GetFieldPointer< MRField > ( grid, usdPara->residual[ 0 ] );
    MRField * res1 = GetFieldPointer< MRField > ( grid, usdPara->residual[ 1 ] );
    MRField * res2 = GetFieldPointer< MRField > ( grid, usdPara->residual[ 2 ] );

    HXWrite( ActionState::dataBook, res1 );
    HXWrite( ActionState::dataBook, res2 );
}

void Restart::InitUnsteady( int sTid )
{
    FieldManager * fieldManager = FieldFactory::GetFieldManager( sTid );
    UsdPara * usdPara = fieldManager->usdPara;
    Grid * grid = Zone::GetGrid();

    MRField * q  = GetFieldPointer< MRField > ( grid, usdPara->flow[ 0 ] );
    MRField * q1 = GetFieldPointer< MRField > ( grid, usdPara->flow[ 1 ] );
    MRField * q2 = GetFieldPointer< MRField > ( grid, usdPara->flow[ 2 ] );

    SetField( q1, q );
    SetField( q2, q );

    MRField * res  = GetFieldPointer< MRField > ( grid, usdPara->residual[ 0 ] );
    MRField * res1 = GetFieldPointer< MRField > ( grid, usdPara->residual[ 1 ] );
    MRField * res2 = GetFieldPointer< MRField > ( grid, usdPara->residual[ 2 ] );

    SetField( res , 0.0 );
    SetField( res1, res );
    SetField( res2, res );
}

void Restart::Read( int sTid )
{
    ActionState::dataBook->MoveToBegin();

    ReadRestartHeader();

    this->ReadUnsteady( sTid );

    RwInterface( sTid, GREAT_READ );
}


void Restart::Dump( int sTid )
{
    ActionState::dataBook->MoveToBegin();

    DumpRestartHeader();

    this->DumpUnsteady( sTid );

    RwInterface( sTid, GREAT_WRITE );
}

void Restart::InitRestart( int sTid )
{
    Iteration::outerSteps = 0;
    ctrl.currTime = 0.0;
}

void Restart::InitinsRestart( int sTid )
{
	Iteration::outerSteps = 0;
	ctrl.currTime = 0.0;
}

void ReadRestartHeader()
{
    HXRead( ActionState::dataBook, Iteration::outerSteps );
    HXRead( ActionState::dataBook, ctrl.currTime );
}

void ReadinsRestartHeader()
{
	HXRead(ActionState::dataBook, Iteration::outerSteps);
	HXRead(ActionState::dataBook, ctrl.currTime);
}

void DumpRestartHeader()
{
    HXWrite( ActionState::dataBook, Iteration::outerSteps );
    HXWrite( ActionState::dataBook, ctrl.currTime );
}

void RwInterface( int sTid, int readOrWrite )
{
    Grid * grid = Zone::GetGrid();
    InterFace * interFace = grid->interFace;

    if ( ! IsValid( interFace ) ) return;

     VarNameSolver * varNameSolver = VarNameFactory::GetVarNameSolver( sTid, INTERFACE_GRADIENT_DATA );
    StringField fieldNameList = varNameSolver->data;

    for ( int ghostId = MAX_GHOST_LEVELS - 1; ghostId >= 0; -- ghostId )
    {
        RwInterfaceRecord( interFace->dataSend[ ghostId ], fieldNameList, readOrWrite );
    }

    for ( int ghostId = MAX_GHOST_LEVELS - 1; ghostId >= 0; -- ghostId )
    {
        RwInterfaceRecord( interFace->dataRecv[ ghostId ], fieldNameList, readOrWrite );
    }
}

void RwInterfaceRecord( DataStorage * storage, StringField & fieldNameList, int readOrWrite )
{
    if ( readOrWrite == GREAT_READ )
    {
        ReadFieldRecord( storage, fieldNameList );
    }
    else if ( readOrWrite == GREAT_WRITE )
    {
        WriteFieldRecord( storage, fieldNameList );
    }
    else if ( readOrWrite == GREAT_ZERO )
    {
        ZeroFieldRecord( storage, fieldNameList );
    }
}

void ReadFieldRecord( DataStorage * storage, StringField & fieldNameList )
{
    for ( int iField = 0; iField < fieldNameList.size(); ++ iField )
    {
        string & filedName = fieldNameList[ iField ];
        MRField * field = ONEFLOW::GetFieldPointer< MRField >( storage, filedName );

        HXRead( ActionState::dataBook, field );
    }
}

void WriteFieldRecord( DataStorage * storage, StringField & fieldNameList )
{
    for ( int iField = 0; iField < fieldNameList.size(); ++ iField )
    {
        string & filedName = fieldNameList[ iField ];
        MRField * field = ONEFLOW::GetFieldPointer< MRField >( storage, filedName );

        HXWrite( ActionState::dataBook, field );
    }
}

void ZeroFieldRecord( DataStorage * storage, StringField & fieldNameList )
{
    for ( int iField = 0; iField < fieldNameList.size(); ++ iField )
    {
        string & filedName = fieldNameList[ iField ];
        MRField * field = ONEFLOW::GetFieldPointer< MRField >( storage, filedName );

        SetField( field, 0.0 );
    }
}


EndNameSpace