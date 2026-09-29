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

    OneFLOW is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\---------------------------------------------------------------------------*/

#include "UnsteadyFieldView.h"
#include "FieldManager.h"
#include "FieldWrap.h"
#include "DataBase.h"
#include "Zone.h"
#include "Fatal.h"
#include "UnsGrid.h"

BeginNameSpace( ONEFLOW )


UnsteadyFieldView::UnsteadyFieldView()
{
}

UnsteadyFieldView::~UnsteadyFieldView()
{
}

void UnsteadyFieldView::BindFields( int solverType )
{
    FieldManager * fieldManager =
        FieldManagerRegistry::GetFieldManager(
            solverType );

    if ( fieldManager == nullptr )
    {
        Fatal(
            "FieldManager is not registered for solverType" );
    }

    this->BindFields(
        fieldManager->GetUnsteadyFieldNames() );
}

void UnsteadyFieldView::BindFields(
    const UnsteadyFieldNames & fieldNames )
{
    // Build unsteady access view from registered field names.
    // This class does not allocate field storage.
    UnsGrid * grid = Zone::GetUnsGrid();

    this->flow.resize(
        fieldNames.flow.size() );

    for ( std::size_t i = 0;
        i < fieldNames.flow.size();
        ++ i )
    {
        this->flow[ i ] =
            GetFieldPointer< MRField >(
                grid,
                fieldNames.flow[ i ] );
    }

    this->residual.resize(
        fieldNames.residual.size() );

    if ( this->flow.size() < 3 )
    {
        Fatal(
            "Unsteady flow fields require at least 3 time levels." );
    }

    if ( this->residual.size() < 3 )
    {
        Fatal(
            "Unsteady residual fields require at least 3 time levels." );
    }

    for ( std::size_t i = 0;
        i < fieldNames.residual.size();
        ++ i )
    {
        this->residual[ i ] =
            GetFieldPointer< MRField >(
                grid,
                fieldNames.residual[ i ] );
    }
}

MRField * UnsteadyFieldView::GetFlow( std::size_t level )
{
    return this->flow[ level ];
}

MRField * UnsteadyFieldView::GetResidual( std::size_t level )
{
    return this->residual[ level ];
}

MRField * UnsteadyFieldView::GetFlow( HistoryLevel level )
{
    return this->GetFlow(
        static_cast< std::size_t >( level ) );
}

MRField * UnsteadyFieldView::GetResidual( HistoryLevel level )
{
    return this->GetResidual(
        static_cast< std::size_t >( level ) );
}


EndNameSpace
