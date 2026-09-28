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

#include "InterfaceTaskReg.h"
#include "InterField.h"
#include "Fatal.h"
#include "ActionState.h"
#include "HXMath.h"
#include "SolverDef.h"
#include "SolverInfo.h"
#include "Restart.h"
#include "Unsteady.h"
#include "UnsteadyImp.h"
#include "Update.h"
#include "FieldWrap.h"
#include "CmxTask.h"
#include "DataBase.h"
#include "DataBook.h"
#include "Lusgs.h"
#include "Lhs.h"
#include "FieldManager.h"
#include "SolverState.h"
#include "Zone.h"
#include "Grid.h"
#include "GridState.h"
#include "UnsGrid.h"
#include "InterFace.h"
#include "RegisterUtils.h"
#include "FieldRecord.h"
#include "UResidual.h"
#include "UNsCom.h"
#include "TaskRegister.h"

BeginNameSpace( ONEFLOW )

REGISTER_TASK( RegisterInterfaceTask )

void RegisterInterfaceTask()
{
    REGISTER_DATA_CLASS( CalcInterfaceGrad );
    REGISTER_DATA_CLASS( UploadInterfaceData );
    REGISTER_DATA_CLASS( DownloadInterfaceData );
    REGISTER_DATA_CLASS( PrepareInterfaceField );
    REGISTER_DATA_CLASS( PrepareOversetInterfaceField );
    REGISTER_DATA_CLASS( ExchangeInterfaceDq );
}


void CalcInterfaceGrad( StringField & data )
{
    ;
}

//void UploadInterfaceData( StringField & data )
//{
//    int solverType = SolverState::solverType;
//
//    FieldManager * fieldManager = FieldManagerRegistry::GetFieldManager( solverType );
//
//    fieldManager->GetInterfaceFieldProperty().UploadInterfaceValue();
//
//}

//void DownloadInterfaceData( StringField & data )
//{
//    int solverType = SolverState::solverType;
//
//    FieldManager * fieldManager = FieldManagerRegistry::GetFieldManager( solverType );
//    fieldManager->GetInterfaceFieldProperty().DownloadInterfaceValue();
//}

void UploadInterfaceData( StringField & data )
{
    Grid * gridIn = Zone::GetGrid();
    if ( ! ONEFLOW::IsUnsGrid( gridIn->type ) )
    {
        return;
    }

    UnsGrid * grid = ONEFLOW::UnsGridCast( gridIn );

    int solverType = SolverState::solverType;
    FieldManager * fieldManager =
        FieldManagerRegistry::GetFieldManager( solverType );

    const auto & fieldData =
        fieldManager->GetInterfaceFieldProperty().GetData();

    for ( const auto & [ fieldName, nEqu ] : fieldData )
    {
        MRField * field =
            ONEFLOW::GetFieldPointer< MRField >(
                grid,
                fieldName );

        if ( field == nullptr )
        {
            Fatal(
                "Grid field is not allocated for interface upload: "
                + fieldName );
        }

        ONEFLOW::UploadInterfaceValue(
            grid,
            field,
            fieldName,
            nEqu );
    }
}

void DownloadInterfaceData( StringField & data )
{
    Grid * gridIn = Zone::GetGrid();
    if ( ! ONEFLOW::IsUnsGrid( gridIn->type ) )
    {
        return;
    }

    UnsGrid * grid = ONEFLOW::UnsGridCast( gridIn );

    int solverType = SolverState::solverType;
    FieldManager * fieldManager =
        FieldManagerRegistry::GetFieldManager( solverType );

    const auto & fieldData =
        fieldManager->GetInterfaceFieldProperty().GetData();

    for ( const auto & [ fieldName, nEqu ] : fieldData )
    {
        MRField * field =
            ONEFLOW::GetFieldPointer< MRField >(
                grid,
                fieldName );

        if ( field == nullptr )
        {
            Fatal(
                "Grid field is not allocated for interface download: "
                + fieldName );
        }

        ONEFLOW::DownloadInterfaceValue(
            grid,
            field,
            fieldName,
            nEqu );
    }
}

void PrepareInterfaceField( StringField & data )
{
    Grid * grid = Zone::GetGrid();
    InterFace * interFace = grid->interFace;
    if ( ! ONEFLOW::IsValid( interFace ) ) return;

    int solverType = SolverState::solverType;
    int iFk  = ( * interfaceMap )[ data[ 0 ] ];
    int iSr  = ( * sendRecvMap )[ data[ 1 ] ];

    FieldRecord * fieldRecord = new FieldRecord();

    PrepareInterfaceFieldRecord( solverType, iFk, iSr, fieldRecord );

    SetInterfaceFieldData( iSr, fieldRecord );

    delete fieldRecord;
}

void PrepareOversetInterfaceField( StringField & data )
{
    ;
}

void ExchangeInterfaceDq( StringField & data )
{
    ;
}

EndNameSpace
