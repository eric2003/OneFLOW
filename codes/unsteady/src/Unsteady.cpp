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
    if ( fieldNames.flow.size() < 3 )
    {
        Fatal(
            "Unsteady requires at least 3 flow time levels." );
    }

    if ( fieldNames.residual.size() < 3 )
    {
        Fatal(
            "Unsteady requires at least 3 residual time levels." );
    }

    this->field.BindFields( grid, fieldNames );
}

MRField * Unsteady::GetFlow( Unsteady::HistoryLevel level )
{
    return this->field.GetFlow(
        static_cast< std::size_t >( level ) );
}

MRField * Unsteady::GetResidual( Unsteady::HistoryLevel level )
{
    return this->field.GetResidual(
        static_cast< std::size_t >( level ) );
}

void Unsteady::UpdateUnsteady()
{
    // Reuse the persistent field view instead of rebuilding a second view.
    UnsteadyFieldView & fieldView = this->field;

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
