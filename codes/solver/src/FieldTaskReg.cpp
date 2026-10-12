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

#include "FieldTaskReg.h"
#include "NodeField.h"
#include "Fatal.h"
#include "DataBase.h"
#include "Grid.h"
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
#include "CmxTask.h"
#include "DataBase.h"
#include "DataBook.h"
#include "Lusgs.h"
#include "Lhs.h"
#include "FieldManager.h"
#include "SolverState.h"
#include "Zone.h"
#include "Grid.h"
#include "UnsGrid.h"
#include "InterFace.h"
#include "RegisterUtils.h"
#include "FieldRecord.h"
#include "UResidual.h"
#include "UNsCom.h"
#include "TaskRegister.h"
#include <map>
#include <iostream>


BeginNameSpace( ONEFLOW )

REGISTER_TASK( RegisterFieldTask )

void RegisterFieldTask()
{
    REGISTER_DATA_CLASS( LoadResiduals );
    REGISTER_DATA_CLASS( ZeroResiduals );
    REGISTER_DATA_CLASS( LoadQ );
    REGISTER_DATA_CLASS( StoreRHS );
    REGISTER_DATA_CLASS( SetField );
    REGISTER_DATA_CLASS( ConvertNodeFieldToCell );
}

void LoadResiduals( StringField & data )
{
    int solverType = SolverState::solverType;
    int fieldId = ( * fieldIdMap )[ data[ 0 ] ];

    UnsGrid * grid = Zone::GetUnsGrid();

    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( solverType );
    std::string & residualName = solverInfo->residualName;

    MRField * resField = FieldHome::GetResidualField( grid );
    MRField * rhsField = FieldHome::GetUnsField( fieldId );

    NegField( resField, rhsField );

}

void ZeroResiduals( StringField & data )
{
    int solverType = SolverState::solverType;

    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( solverType );
    std::string & fieldName = solverInfo->residualName;

    FieldHome::SetField( fieldName, zero );
}

void LoadQ( StringField & data )
{
    std::string & fName = data[ 0 ];

    int solverType = SolverState::solverType;
    int fieldId = ( * fieldIdMap )[ fName ];

    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( solverType );

    std::string & fieldName = solverInfo->gradString[ 0 ];

    FieldHome::SetField( fieldId, fieldName, FLOW_RHS_ORDER );
}

void SetField( StringField & data )
{
    int solverType = SolverState::solverType;
    std::string & fieldName   = data[ 0 ];
    std::string & valueString = data[ 1 ];

    Real value = StringToDigit< Real >( valueString );
    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( solverType );

    FieldHome::SetField( fieldName, value );
}

void StoreRHS( StringField & data )
{
    std::string & fName = data[ 0 ];

    int solverType = SolverState::solverType;
    int fieldId = ( * fieldIdMap )[ fName ];

    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( solverType );
    std::string & residualName = solverInfo->residualName;

    FieldHome::SetField( fieldId, residualName, FLOW_RHS_ORDER );
}


void ConvertNodeFieldToCell( StringField & data )
{
    // data[0]: source cell-centered field name on current zone grid.
    // data[1]: optional target cell field name (default: same as source — in-place update).
    // Real coupling: cell -> node -> cell, then write back into the named target field.
    if ( data.empty() )
    {
        return;
    }
    const std::string & srcName = data[ 0 ];
    const std::string & dstName = ( data.size() > 1 ) ? data[ 1 ] : data[ 0 ];

    auto nodeField = InterpolateCellToNode( srcName );

    Grid * grid = Zone::GetGrid();
    MRField * target = ONEFLOW::GetFieldPointer< MRField >( grid, dstName );
    if ( target == nullptr )
    {
        Fatal( "ConvertNodeFieldToCell: target field not allocated: " + dstName );
    }

    const int nEqu = nodeField->GetNEqu();
    if ( target->GetNEqu() < nEqu )
    {
        Fatal( "ConvertNodeFieldToCell: target field nEqu smaller than source" );
    }
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        InterpolateNodeToCellForComponent( ( * target )[ iEqu ], ( * nodeField )[ iEqu ] );
    }
}

EndNameSpace
