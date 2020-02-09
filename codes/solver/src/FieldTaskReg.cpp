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

#include "FieldTaskReg.h"
#include "InterField.h"
#include "ActionState.h"
#include "HXMath.h"
#include "SolverDef.h"
#include "SolverInfo.h"
#include "Restart.h"
#include "Unsteady.h"
#include "UnsteadyImp.h"
#include "Update.h"
#include "FieldWrap.h"
#include "FieldAlloc.h"
#include "CmxTask.h"
#include "DataBase.h"
#include "DataBook.h"
#include "Lusgs.h"
#include "Lhs.h"
#include "FieldImp.h"
#include "FieldWrap.h"
#include "SolverState.h"
#include "Zone.h"
#include "Grid.h"
#include "UnsGrid.h"
#include "InterFace.h"
#include "RegisterUtil.h"
#include "FieldRecord.h"
#include "UVisualize.h"
#include "UResidual.h"
#include "UNsCom.h"
#include "TaskRegister.h"
#include <map>
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

REGISTER_TASK( RegisterFieldTask )

void RegisterFieldTask()
{
    REGISTER_DATA_CLASS( LoadResiduals );
    REGISTER_DATA_CLASS( ZeroResiduals );
    REGISTER_DATA_CLASS( LoadQ );
    REGISTER_DATA_CLASS( StoreRHS );
    REGISTER_DATA_CLASS( SetField );
}

void LoadResiduals( StringField & data )
{
    int sTid = SolverState::tid;
    int fieldId = ( * fieldIdMap )[ data[ 0 ] ];

    UnsGrid * grid = Zone::GetUnsGrid();

    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( sTid );
    string & residualName = solverInfo->residualName;

    MRField * resField = ONEFLOW::GetFieldPointer< MRField >( grid, residualName );
    MRField * rhsField = FieldHome::GetUnsField( fieldId );

    NegField( resField, rhsField );

}

void ZeroResiduals( StringField & data )
{
    int sTid = SolverState::tid;

    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( sTid );
    string & fieldName = solverInfo->residualName;

    FieldHome::SetField( fieldName, zero );
}

void LoadQ( StringField & data )
{
    string & fName = data[ 0 ];

    int sTid = SolverState::tid;
    int fieldId = ( * fieldIdMap )[ fName ];

    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( sTid );

    string & fieldName = solverInfo->gradString[ 0 ];

    FieldHome::SetField( fieldId, fieldName, FLOW_RHS_ORDER );
}

void SetField( StringField & data )
{
    int sTid = SolverState::tid;
    string & fieldName   = data[ 0 ];
    string & valueString = data[ 1 ];

    Real value = StringToDigit< Real >( valueString );
    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( sTid );

    FieldHome::SetField( fieldName, value );
}

void StoreRHS( StringField & data )
{
    string & fName = data[ 0 ];

    int sTid = SolverState::tid;
    int fieldId = ( * fieldIdMap )[ fName ];

    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( sTid );
    string & residualName = solverInfo->residualName;

    FieldHome::SetField( fieldId, residualName, FLOW_RHS_ORDER );
}

EndNameSpace