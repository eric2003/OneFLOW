/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
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
#include <memory>
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
#include "DataStorage.h"
#include "Iteration.h"
#include "FieldManager.h"
#include "FieldWrap.h"
#include "UnsteadyFieldView.h"
#include "Fatal.h"
#include "RegisterUtils.h"
#include "INsRestart.h"

BeginNameSpace( ONEFLOW )

namespace
{
    void BindUnsteadyFields(
        UnsteadyFieldView & fieldView,
        int solverType )
    {
        FieldManager * fieldManager =
            FieldManagerRegistry::GetFieldManager(
                solverType );

        if ( fieldManager == nullptr )
        {
            Fatal(
                "FieldManager is not registered for solverType" );
        }

        UnsGrid * grid = Zone::GetUnsGrid();

        fieldView.BindFields(
            grid,
            fieldManager->GetUnsteadyFieldNames() );
    }
}

std::unique_ptr<Restart> CreateRestart( int solverType )
{
    if ( solverType == NS_SOLVER )
    {
        return CreateNsRestart();
    }
    else if ( solverType == INC_NS_SOLVER )
    {
        return CreateINsRestart();
    }
    else if ( solverType == TURB_SOLVER )
    {
        return CreateTurbRestart();
    }

    return nullptr;
}

Restart::Restart()
{
    ;
}

Restart::~Restart()
{
    ;
}

void Restart::ReadUnsteady( int solverType )
{
    UnsteadyFieldView fieldView;
    BindUnsteadyFields( fieldView, solverType );

    // The current level is reconstructed from the first stored history level.
    HXRead(
        ActionState::dataBook,
        fieldView.GetFlow( 1 ) );

    SetField(
        fieldView.GetFlow( 0 ),
        fieldView.GetFlow( 1 ) );

    for ( std::size_t level = 2;
        level < fieldView.GetFlowCount();
        ++ level )
    {
        HXRead(
            ActionState::dataBook,
            fieldView.GetFlow( level ) );
    }

    // Residual history follows the same restart layout as flow history.
    HXRead(
        ActionState::dataBook,
        fieldView.GetResidual( 1 ) );

    SetField(
        fieldView.GetResidual( 0 ),
        fieldView.GetResidual( 1 ) );

    for ( std::size_t level = 2;
        level < fieldView.GetResidualCount();
        ++ level )
    {
        HXRead(
            ActionState::dataBook,
            fieldView.GetResidual( level ) );
    }
}

void Restart::DumpUnsteady( int solverType )
{
    UnsteadyFieldView fieldView;
    BindUnsteadyFields( fieldView, solverType );

    // Keep the current level out of the restart stream.
    // It is reconstructed from the first stored history level on read.
    for ( std::size_t level = 1;
        level < fieldView.GetFlowCount();
        ++ level )
    {
        HXWrite(
            ActionState::dataBook,
            fieldView.GetFlow( level ) );
    }

    for ( std::size_t level = 1;
        level < fieldView.GetResidualCount();
        ++ level )
    {
        HXWrite(
            ActionState::dataBook,
            fieldView.GetResidual( level ) );
    }
}

void Restart::InitUnsteady( int solverType )
{
    UnsteadyFieldView fieldView;
    BindUnsteadyFields( fieldView, solverType );

    // Initialize every configured history level from the current field.
    for ( std::size_t level = 1;
        level < fieldView.GetFlowCount();
        ++ level )
    {
        SetField(
            fieldView.GetFlow( level ),
            fieldView.GetFlow( 0 ) );
    }

    SetField(
        fieldView.GetResidual( 0 ),
        0.0 );

    for ( std::size_t level = 1;
        level < fieldView.GetResidualCount();
        ++ level )
    {
        SetField(
            fieldView.GetResidual( level ),
            fieldView.GetResidual( 0 ) );
    }
}

void Restart::Read( int solverType )
{
    ActionState::dataBook->MoveToBegin();

    ReadRestartHeader();

    this->ReadUnsteady( solverType );

    RwInterface( solverType, GREAT_READ );
}

void Restart::Dump( int solverType )
{
    ActionState::dataBook->MoveToBegin();

    DumpRestartHeader();

    this->DumpUnsteady( solverType );

    RwInterface( solverType, GREAT_WRITE );
}

void Restart::InitRestart( int solverType )
{
    Iteration::outerSteps = 0;
    ctrl.currTime = 0.0;
}

void Restart::InitinsRestart( int solverType )
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

void RwInterface( int solverType, int readOrWrite )
{
    Grid & grid = Zone::GetGridReference();
    InterFace * interFace = grid.interFace.get();

    if ( ! IsValid( interFace ) ) return;

    VarNameSolver * varNameSolver = VarNameFactory::GetVarNameSolver( solverType, INTERFACE_GRADIENT_DATA );
    StringField fieldNameList = varNameSolver->data;

    for ( int ghostId = MAX_GHOST_LEVELS - 1; ghostId >= 0; -- ghostId )
    {
        RwInterfaceRecord( &interFace->GetSendStorage( ghostId ), fieldNameList, readOrWrite );
    }

    for ( int ghostId = MAX_GHOST_LEVELS - 1; ghostId >= 0; -- ghostId )
    {
        RwInterfaceRecord( &interFace->GetRecvStorage( ghostId ), fieldNameList, readOrWrite );
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
        std::string & filedName = fieldNameList[ iField ];
        MRField * field = ONEFLOW::GetFieldPointer< MRField >( storage, filedName );

        HXRead( ActionState::dataBook, field );
    }
}

void WriteFieldRecord( DataStorage * storage, StringField & fieldNameList )
{
    for ( int iField = 0; iField < fieldNameList.size(); ++ iField )
    {
        std::string & filedName = fieldNameList[ iField ];
        MRField * field = ONEFLOW::GetFieldPointer< MRField >( storage, filedName );

        HXWrite( ActionState::dataBook, field );
    }
}

void ZeroFieldRecord( DataStorage * storage, StringField & fieldNameList )
{
    for ( int iField = 0; iField < fieldNameList.size(); ++ iField )
    {
        std::string & filedName = fieldNameList[ iField ];
        MRField * field = ONEFLOW::GetFieldPointer< MRField >( storage, filedName );

        SetField( field, 0.0 );
    }
}


EndNameSpace
