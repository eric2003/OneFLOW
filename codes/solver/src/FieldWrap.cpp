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

#include "FieldWrap.h"
#include <memory>
#include "Fatal.h"
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


BeginNameSpace( ONEFLOW )

MRField * FieldWrap::GetUnsField()
{
    return view;
}

const MRField * FieldWrap::GetConstUnsField() const
{
    return view;
}

void FieldWrap::SetUnsField( MRField * unsField, bool deleteFlag )
{
    if ( deleteFlag )
    {
        owned.reset( unsField );
        view = owned.get();
    }
    else
    {
        owned.reset();
        view = unsField;
    }
}

void FieldWrap::SetOwnedField( std::unique_ptr<MRField> field )
{
    owned = std::move( field );
    view = owned.get();
}

std::unique_ptr<FieldWrap> FieldHome::CreateField()
{
    SolverState::SetSolverTypeBySolverIndex( SolverState::solverIndex );
    return FieldHome::CreateField( SolverState::solverType, GridState::gridLevel );
}

std::unique_ptr<FieldWrap> FieldHome::CreateField( int solverType )
{
    return FieldHome::CreateField( solverType, GridState::gridLevel );
}

std::unique_ptr<FieldWrap> FieldHome::CreateField( int solverType, int level )
{
    SolverInfo * info = SolverInfoFactory::GetSolverInfo( solverType );

    Grid & grid = Zone::GetGridReference();

    int nTCell = grid.nCells + grid.nBFaces;

    auto field = std::make_unique<MRField>( info->nTEqu, nTCell );

    auto fieldWrap = std::make_unique<FieldWrap>();

    fieldWrap->SetOwnedField( std::move( field ) );

    return fieldWrap;
}

std::unique_ptr<FieldWrap> FieldHome::GetFieldWrap( const std::string & fieldName )
{
    Grid & grid = Zone::GetGridReference();

    MRField * field = ONEFLOW::GetFieldPointer< MRField >( &grid, fieldName );

    auto fieldWrap = std::make_unique<FieldWrap>();

    fieldWrap->SetUnsField( field );

    return fieldWrap;
}

void FieldHome::SetField( const std::string & fieldName, Real value )
{
    // Direct path: no temporary FieldWrap for constant init from alloc/init.txt.
    Grid & grid = Zone::GetGridReference();
    MRField * field =
        ONEFLOW::GetFieldPointer< MRField >( &grid, fieldName );

    if ( field == nullptr )
    {
        Fatal(
            "Field is not allocated: " + fieldName );
    }

    ONEFLOW::SetField( field, value );
}

void FieldHome::SetField( int fieldId, const std::string & fieldName, int orderFlag )
{
    FieldHome::SetUnsField( fieldId, fieldName, orderFlag );
}

void FieldHome::SetUnsField( int fieldId, const std::string & fieldName, int orderFlag )
{
    Grid & grid = Zone::GetGridReference();
    UnsGrid * unsGrid = ONEFLOW::UnsGridCast( &grid );

    MRField * sField, * tField;
    FieldHome::GetSourceTargetField( unsGrid, fieldId, fieldName, sField, tField, orderFlag );

    ONEFLOW::SetField( tField, sField );
}

FieldWrap * FieldHome::GetFieldWrap( int fieldId )
{
    FieldWrap * fieldWrap = BgField::GetFieldWrap( ZoneState::zid, SolverState::solverIndex, fieldId, GridState::gridLevel );
    return fieldWrap;
}

MRField * FieldHome::GetUnsField( int fieldId )
{
    FieldWrap * fieldWrap = FieldHome::GetFieldWrap( fieldId );
    MRField * field = fieldWrap->GetUnsField();
    return field;
}


const FieldWrap * FieldHome::GetConstFieldWrap( int fieldId )
{
    return BgField::GetConstFieldWrap( ZoneState::zid, SolverState::solverIndex, fieldId, GridState::gridLevel );
}

const MRField * FieldHome::GetConstUnsField( int fieldId )
{
    const FieldWrap * fieldWrap = FieldHome::GetConstFieldWrap( fieldId );
    return fieldWrap->GetConstUnsField();
}

const MRField * FieldHome::GetConstUnsField( const Grid * grid, const std::string & fieldName )
{
    return ONEFLOW::GetFieldPointer< MRField >( grid, fieldName );
}

const std::string & FieldHome::GetResidualName()
{
    return FieldHome::GetResidualName( SolverState::solverType );
}

const std::string & FieldHome::GetResidualName( int solverType )
{
    return SolverInfoFactory::GetSolverInfo( solverType )->residualName;
}

MRField * FieldHome::GetResidualField( Grid * grid )
{
    return FieldHome::GetResidualField( grid, SolverState::solverType );
}

MRField * FieldHome::GetResidualField( Grid * grid, int solverType )
{
    return FieldHome::GetUnsField( grid, FieldHome::GetResidualName( solverType ) );
}

const MRField * FieldHome::GetConstResidualField( const Grid * grid )
{
    return FieldHome::GetConstResidualField( grid, SolverState::solverType );
}

const MRField * FieldHome::GetConstResidualField( const Grid * grid, int solverType )
{
    return FieldHome::GetConstUnsField( grid, FieldHome::GetResidualName( solverType ) );
}

const std::string & FieldHome::GetSolutionName()
{
    return FieldHome::GetSolutionName( SolverState::solverType );
}

const std::string & FieldHome::GetSolutionName( int solverType )
{
    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( solverType );
    if ( ! solverInfo->solutionGradNames.empty() )
    {
        return solverInfo->solutionGradNames[ 0 ];
    }
    static const std::string kDefaultQ = "q";
    return kDefaultQ;
}

MRField * FieldHome::GetSolutionField( Grid * grid )
{
    return FieldHome::GetSolutionField( grid, SolverState::solverType );
}

MRField * FieldHome::GetSolutionField( Grid * grid, int solverType )
{
    return FieldHome::GetUnsField( grid, FieldHome::GetSolutionName( solverType ) );
}

const MRField * FieldHome::GetConstSolutionField( const Grid * grid )
{
    return FieldHome::GetConstUnsField( grid, FieldHome::GetSolutionName() );
}

const std::string & FieldHome::GetGradName( int index )
{
    return FieldHome::GetGradName( SolverState::solverType, index );
}

const std::string & FieldHome::GetGradName( int solverType, int index )
{
    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( solverType );
    if ( index >= 0 && index < static_cast< int >( solverInfo->solutionGradNames.size() ) )
    {
        return solverInfo->solutionGradNames[ index ];
    }
    static const std::string kEmpty;
    return kEmpty;
}

MRField * FieldHome::GetGradField( Grid * grid, int index )
{
    return FieldHome::GetGradField( grid, SolverState::solverType, index );
}

MRField * FieldHome::GetGradField( Grid * grid, int solverType, int index )
{
    const std::string & name = FieldHome::GetGradName( solverType, index );
    if ( name.empty() ) return nullptr;
    return FieldHome::GetUnsField( grid, name );
}

const MRField * FieldHome::GetConstGradField( const Grid * grid, int index )
{
    const std::string & name = FieldHome::GetGradName( index );
    if ( name.empty() ) return nullptr;
    return FieldHome::GetConstUnsField( grid, name );
}

const std::string & FieldHome::GetUpdateDeltaName()
{
    return FieldHome::GetUpdateDeltaName( SolverState::solverType );
}

const std::string & FieldHome::GetUpdateDeltaName( int solverType )
{
    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( solverType );
    if ( solverInfo->implicitFieldNames.size() >= 2 )
    {
        return solverInfo->implicitFieldNames[ 1 ];
    }
    static const std::string kDefaultDq = "dq";
    return kDefaultDq;
}

MRField * FieldHome::GetUpdateDeltaField( Grid * grid )
{
    return FieldHome::GetUpdateDeltaField( grid, SolverState::solverType );
}

MRField * FieldHome::GetUpdateDeltaField( Grid * grid, int solverType )
{
    return FieldHome::GetUnsField( grid, FieldHome::GetUpdateDeltaName( solverType ) );
}

const std::string & FieldHome::GetBcGradName( int index )
{
    return FieldHome::GetBcGradName( SolverState::solverType, index );
}

const std::string & FieldHome::GetBcGradName( int solverType, int index )
{
    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( solverType );
    if ( index >= 0 && index < static_cast< int >( solverInfo->bcGradFieldNames.size() ) )
    {
        return solverInfo->bcGradFieldNames[ index ];
    }
    static const std::string kEmpty;
    return kEmpty;
}

MRField * FieldHome::GetBcGradField( Grid * grid, int index )
{
    return FieldHome::GetBcGradField( grid, SolverState::solverType, index );
}

MRField * FieldHome::GetBcGradField( Grid * grid, int solverType, int index )
{
    const std::string & name = FieldHome::GetBcGradName( solverType, index );
    if ( name.empty() ) return nullptr;
    return FieldHome::GetUnsField( grid, name );
}

const std::string & FieldHome::GetStageSolutionName( int index )
{
    return FieldHome::GetStageSolutionName( SolverState::solverType, index );
}

const std::string & FieldHome::GetStageSolutionName( int solverType, int index )
{
    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( solverType );
    if ( index >= 0 && index < static_cast< int >( solverInfo->stageSolutionNames.size() ) )
    {
        return solverInfo->stageSolutionNames[ index ];
    }
    static const std::string kEmpty;
    return kEmpty;
}

MRField * FieldHome::GetStageSolutionField( Grid * grid, int index )
{
    return FieldHome::GetStageSolutionField( grid, SolverState::solverType, index );
}

MRField * FieldHome::GetStageSolutionField( Grid * grid, int solverType, int index )
{
    const std::string & name = FieldHome::GetStageSolutionName( solverType, index );
    if ( name.empty() )
    {
        return nullptr;
    }
    return FieldHome::GetUnsField( grid, name );
}

const std::string & FieldHome::GetStageResidualName( int index )
{
    return FieldHome::GetStageResidualName( SolverState::solverType, index );
}

const std::string & FieldHome::GetStageResidualName( int solverType, int index )
{
    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( solverType );
    if ( index >= 0 && index < static_cast< int >( solverInfo->stageResidualNames.size() ) )
    {
        return solverInfo->stageResidualNames[ index ];
    }
    static const std::string kEmpty;
    return kEmpty;
}

MRField * FieldHome::GetStageResidualField( Grid * grid, int index )
{
    return FieldHome::GetStageResidualField( grid, SolverState::solverType, index );
}

MRField * FieldHome::GetStageResidualField( Grid * grid, int solverType, int index )
{
    const std::string & name = FieldHome::GetStageResidualName( solverType, index );
    if ( name.empty() )
    {
        return nullptr;
    }
    return FieldHome::GetUnsField( grid, name );
}

MRField * FieldHome::GetUnsField( Grid * grid, const std::string & fieldName )
{
    MRField * field = ONEFLOW::GetFieldPointer< MRField >( grid, fieldName );
    return field;
}

void FieldHome::GetSourceTargetField( Grid * grid, int fieldId, const std::string & fieldName, MRField *& sField, MRField *& tField, int orderFlag )
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
