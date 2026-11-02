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
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
    License for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "Unsteady.h"
#include "FieldWrap.h"
#include "Fatal.h"

BeginNameSpace( ONEFLOW )

void Unsteady::BindFields(
    UnsGrid * grid,
    const UnsteadyFieldNames & fieldNames )
{
    const std::size_t requiredHistoryLevels =
        GetHistoryIndex( HistoryLevel::Old ) + 1;

    if ( fieldNames.flow.size() < requiredHistoryLevels )
    {
        Fatal(
            "Unsteady requires Current, Previous, and Old flow time levels." );
    }

    if ( fieldNames.residual.size() < requiredHistoryLevels )
    {
        Fatal(
            "Unsteady requires Current, Previous, and Old residual time levels." );
    }

    this->fieldView.BindFields( grid, fieldNames );
}

std::size_t Unsteady::GetHistoryIndex( Unsteady::HistoryLevel level )
{
    switch ( level )
    {
    case Unsteady::HistoryLevel::Current:
        return 0;

    case Unsteady::HistoryLevel::Previous:
        return 1;

    case Unsteady::HistoryLevel::Old:
        return 2;
    }

    Fatal( "Invalid unsteady history level" );
    return 0;
}

MRField * Unsteady::GetFlow( Unsteady::HistoryLevel level )
{
    return this->fieldView.GetFlow(
        GetHistoryIndex( level ) );
}

MRField * Unsteady::GetResidual( Unsteady::HistoryLevel level )
{
    return this->fieldView.GetResidual(
        GetHistoryIndex( level ) );
}

void Unsteady::UpdateUnsteady()
{
    // Reuse the persistent field view instead of rebuilding a second view.
    UnsteadyFieldView & fieldView = this->fieldView;

    // Shift from the oldest configured level toward the current level.
    // Reverse order prevents overwriting a history level before it is copied.
    for ( std::size_t level = fieldView.GetFlowCount();
        level > 1;
        -- level )
    {
        SetField(
            fieldView.GetFlow( level - 1 ),
            fieldView.GetFlow( level - 2 ) );
    }
}

EndNameSpace
