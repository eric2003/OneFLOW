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
#include "FieldDefinitionTable.h"
#include "Fatal.h"

BeginNameSpace( ONEFLOW )

void FieldDefinitionTable::AddField(
    const std::string & fieldName,
    int nEqu )
{
    // Policy:
    // - first registration: insert
    // - same name, same nEqu: idempotent (no-op)
    // - same name, different nEqu: configuration conflict
    auto iter = this->data.find( fieldName );

    if ( iter == this->data.end() )
    {
        this->data[ fieldName ] = nEqu;
        return;
    }

    if ( iter->second != nEqu )
    {
        Fatal(
            "Conflicting field definition: "
            + fieldName
            + " existing nEqu="
            + std::to_string( iter->second )
            + " new nEqu="
            + std::to_string( nEqu ) );
    }
}

bool FieldDefinitionTable::HasField(
    const std::string & fieldName ) const
{
    return this->data.find( fieldName ) != this->data.end();
}

int FieldDefinitionTable::GetNEqu(
    const std::string & fieldName ) const
{
    auto iter = this->data.find( fieldName );

    if ( iter == this->data.end() )
    {
        Fatal(
            "Field is not defined: "
            + fieldName );
    }

    return iter->second;
}

bool FieldDefinitionTable::Empty() const
{
    return this->data.empty();
}

const FieldDefinitionTable::Data & FieldDefinitionTable::GetData() const
{
    return this->data;
}

void FieldDefinitionTable::Dump(
    std::ostream & output ) const
{
    for ( FieldDefinitionTable::Data::const_iterator iter =
        this->data.begin();
        iter != this->data.end();
        ++ iter )
    {
        output
            << "    "
            << iter->first
            << "  nEqu="
            << iter->second
            << '\n';
    }
}

FieldDefinitionTable & FieldDefinitionSet::GetFieldDefinition(
    FieldLocation location )
{
    switch ( location )
    {
    case FieldLocation::Inner:
        return innerField;

    case FieldLocation::Face:
        return faceField;

    case FieldLocation::Boundary:
        return bcField;
    }

    Fatal( "Invalid field location" );
    return innerField;
}

const FieldDefinitionTable & FieldDefinitionSet::GetFieldDefinition(
    FieldLocation location ) const
{
    switch ( location )
    {
    case FieldLocation::Inner:
        return innerField;

    case FieldLocation::Face:
        return faceField;

    case FieldLocation::Boundary:
        return bcField;
    }

    Fatal( "Invalid field location" );
    return innerField;
}

EndNameSpace