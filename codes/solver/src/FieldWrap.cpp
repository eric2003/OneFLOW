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

#include "FieldWrap.h"
#include "SolverMap.h"
#include "BgField.h"
#include "Solver.h"
#include "SolverInfo.h"
#include "SolverState.h"
#include "GridState.h"
#include "BgGrid.h"
#include "Zone.h"
#include "ZoneState.h"
#include "Mesh.h"
#include "Grid.h"
#include "UnsGrid.h"
#include "DataBase.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

FieldWrap::FieldWrap()
{
    unsField   = 0;
    deleteFlag = false;
}

FieldWrap::~FieldWrap()
{
    if ( deleteFlag )
    {
        delete unsField;
    }
}

MRField * FieldWrap::GetUnsField()
{
    return unsField;
}

void FieldWrap::SetUnsField( MRField * unsField, bool deleteFlag )
{
    this->unsField = unsField;
    this->deleteFlag = deleteFlag;
}

FieldHome::FieldHome()
{
    ;
}

FieldHome::~FieldHome()
{
    ;
}

FieldWrap * FieldHome::CreateField()
{
    SolverState::SetTidById( SolverState::id );
    return FieldHome::CreateField( SolverState::tid, GridState::gridLevel );
}

FieldWrap * FieldHome::CreateField( int sTid )
{
    return FieldHome::CreateField( sTid, GridState::gridLevel );
}

FieldWrap * FieldHome::CreateField( int sTid, int level )
{
    //Solver * solver = SolverMap::GetSolver( id, level );

    SolverInfo * info = SolverInfoFactory::GetSolverInfo( sTid );
    //SolverInfo * info = solver->info;

    Grid * grid = Zone::GetGrid();

    int nTCell = grid->nCell + grid->nBFace;

    MRField * field = new MRField( info->nTEqu, nTCell );

    FieldWrap * fieldWrap = new FieldWrap();

    fieldWrap->SetUnsField( field, true );

    return fieldWrap;
}

FieldWrap * FieldHome::GetFieldWrap( const string & fieldName )
{
    Grid * grid = Zone::GetGrid();

    MRField * field = ONEFLOW::GetFieldPointer< MRField >( grid, fieldName );

    FieldWrap * fieldWrap = new FieldWrap();

    fieldWrap->SetUnsField( field );

    return fieldWrap;
}

void FieldHome::SetField( const string & fieldName, Real value )
{
    Grid * gridIn = Zone::GetGrid();

    UnsGrid * grid = ONEFLOW::UnsGridCast( gridIn );

    FieldWrap * fieldWrap = FieldHome::GetFieldWrap( fieldName );

    ONEFLOW::SetField( fieldWrap, value );
}

void FieldHome::SetField( int fieldId, const string & fieldName, int orderFlag )
{
    FieldHome::SetUnsField( fieldId, fieldName, orderFlag );
}

void FieldHome::SetUnsField( int fieldId, const string & fieldName, int orderFlag )
{
    Grid * gridIn = Zone::GetGrid();
    UnsGrid * grid = ONEFLOW::UnsGridCast( gridIn );

    MRField * sField, * tField;
    FieldHome::GetSourceTargetField( grid, fieldId, fieldName, sField, tField, orderFlag );

    ONEFLOW::SetField( tField, sField );
}

FieldWrap * FieldHome::GetFieldWrap( int fieldId )
{
    FieldWrap * fieldWrap = BgField::GetFieldWrap( ZoneState::zid, SolverState::id, fieldId, GridState::gridLevel );
    return fieldWrap;
}

MRField * FieldHome::GetUnsField( int fieldId )
{
    FieldWrap * fieldWrap = FieldHome::GetFieldWrap( fieldId );
    MRField * field = fieldWrap->GetUnsField();
    return field;
}

MRField * FieldHome::GetUnsField( Grid * grid, const string & fieldName )
{
    MRField * field = ONEFLOW::GetFieldPointer< MRField >( grid, fieldName );
    return field;
}

void FieldHome::GetSourceTargetField( Grid * grid, int fieldId, const string & fieldName, MRField *& sField, MRField *& tField, int orderFlag )
{
    if ( orderFlag == FLOW_RHS_ORDER )
    {
        sField = FieldHome::GetUnsField( grid, fieldName );
        tField = FieldHome::GetUnsField( fieldId );
    }
    else
    {
        tField = FieldHome::GetUnsField( grid, fieldName );
        sField = FieldHome::GetUnsField( fieldId );
    }
}

void SetField( FieldWrap * fieldWrap, Real value )
{
    MRField * field = fieldWrap->GetUnsField();
    SetField( field, value );
}

void SetField( MRField * field, Real value )
{
    int nTEqu = field->GetNEqu();
    for ( int iEqu = 0; iEqu < nTEqu; ++ iEqu )
    {
        ( * field )[ iEqu ] = value;
    }
}

void SetField( RealField & field, Real value )
{
    field = value;
}

void SetField( MRField * field1, MRField * field2 )
{
    int nTEqu = field1->GetNEqu();
    for ( int iEqu = 0; iEqu < nTEqu; ++ iEqu )
    {
        int nElements = ( * field1 )[ iEqu ].size();
        for ( int iElement = 0; iElement < nElements; ++ iElement )
        {
            ( * field1 )[ iEqu ][ iElement ] = ( * field2 )[ iEqu ][ iElement ];
        }
    }
}

void NegField( MRField * field1, MRField * field2 )
{
    int nTEqu = field1->GetNEqu();
    for ( int iEqu = 0; iEqu < nTEqu; ++ iEqu )
    {
        int nElements = ( * field1 )[ iEqu ].size();

        for ( int iElement = 0; iElement < nElements; ++ iElement )
        {
            ( * field1 )[ iEqu ][ iElement ] = - ( * field2 )[ iEqu ][ iElement ];
        }
    }
}


EndNameSpace