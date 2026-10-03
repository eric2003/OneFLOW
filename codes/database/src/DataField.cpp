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

#include "DataField.h"
#include "DataPointer.h"
#include <memory>

BeginNameSpace( ONEFLOW )

FieldEntry::FieldEntry()
{
    this->name = "";
}

FieldEntry::FieldEntry( const std::string & name, std::unique_ptr<PointerWrap> data )
{
    this->name = name;
    this->data = std::move( data );
}

FieldEntry::~FieldEntry()
{
}

DataField::DataField()
{
}

DataField::~DataField()
{
    Clear();
}

void DataField::Clear()
{
    dataMap.clear();
}

void DataField::UpdateFieldEntry( std::unique_ptr<FieldEntry> fieldEntry )
{
    if ( fieldEntry == nullptr ) return;

    const std::string name = fieldEntry->name;
    auto it = dataMap.find( name );
    if ( it == dataMap.end() )
    {
        dataMap[ name ] = std::move( fieldEntry );
    }
    // else: already exists - discard the new entry (unique_ptr destroys it)
}

FieldEntry * DataField::GetFieldEntry( const std::string & name )
{
    auto it = dataMap.find( name );
    if ( it != dataMap.end() )
    {
        return it->second.get();
    }
    return nullptr;
}

void DataField::DeleteFieldEntry( const std::string & name )
{
    dataMap.erase( name );
}

EndNameSpace
