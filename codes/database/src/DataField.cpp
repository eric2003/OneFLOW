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

BeginNameSpace( ONEFLOW )

DataF::DataF()
{
    this->name = "";
    this->data = nullptr;
}

DataF::DataF( const std::string & name, PointerWrap * data )
{
    this->name = name;
    this->data = data;
}

DataF::~DataF()
{
    // data is owned and deleted by DataField
}

DataField::DataField()
{
    dataMap = new DataMap;
}

DataField::~DataField()
{
    for ( auto & pair : *dataMap )
    {
        delete pair.second->data;   // delete PointerWrap
        delete pair.second;         // delete DataF
    }
    dataMap->clear();
    delete dataMap;
}

void DataField::UpdateDataF( DataF * dataf )
{
    auto it = dataMap->find( dataf->name );
    if ( it == dataMap->end() )
    {
        // Not exist ¡ú take ownership
        ( *dataMap )[ dataf->name ] = dataf;
    }
    else
    {
        // Already exist ¡ú discard the new one
        if ( it->second != dataf )
        {
            delete dataf;
        }
    }
}

DataF * DataField::GetDataF( const std::string & name )
{
    auto it = dataMap->find( name );
    if ( it != dataMap->end() )
    {
        return it->second;
    }
    return nullptr;
}

void DataField::DeleteDataF( const std::string & name )
{
    auto it = dataMap->find( name );
    if ( it != dataMap->end() )
    {
        delete it->second->data;
        delete it->second;
        dataMap->erase( it );
    }
}

EndNameSpace