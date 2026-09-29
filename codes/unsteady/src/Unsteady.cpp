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
#include "Fatal.h"
#include "FieldManager.h"
#include "FieldWrap.h"
#include "Zone.h"

BeginNameSpace( ONEFLOW )

Unsteady::Unsteady()
{
    solverType = -1;
}

Unsteady::~Unsteady()
{
}

void Unsteady::BindFields()
{
    FieldManager * fieldManager =
        FieldManagerRegistry::GetFieldManager(
            this->solverType );

    if ( fieldManager == nullptr )
    {
        Fatal(
            "FieldManager is not registered for solverType" );
    }

    UnsGrid * grid = Zone::GetUnsGrid();

    this->field.BindFields(
        grid,
        fieldManager->GetUnsteadyFieldNames() );
}

void Unsteady::UpdateUnsteady()
{
    // The concrete unsteady object initializes this view in its constructor.
    // Reuse the persistent view instead of rebuilding a second field view.
    UnsteadyFieldView & fieldView = this->field;

    // Shift from the oldest configured level toward the current level.
    // Reverse order prevents overwriting a history level before it is copied.
    for ( std::size_t level = fieldView.flow.size();
        level > 1;
        -- level )
    {
        SetField(
            fieldView.GetFlow( level - 1 ),
            fieldView.GetFlow( level - 2 ) );
    }
}

EndNameSpace
