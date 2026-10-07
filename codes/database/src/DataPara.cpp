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

#include "DataPara.h"
#include <memory>
#include "DataObject.h"
#include "DataBaseType.h"
#include <iostream>

BeginNameSpace( ONEFLOW )

DataEntry::DataEntry()
{
    this->name = "";
    this->data = nullptr;
}

DataEntry::DataEntry( const std::string & name, int type, int size, std::unique_ptr<DataObject> data )
{
    this->name = name;
    this->type = type;
    this->size = size;
    this->data = std::move( data );
}

DataEntry::~DataEntry()
{
}

void DataEntry::Copy( DataEntry * inputData )
{
    if ( inputData == nullptr )
    {
        throw std::invalid_argument(
            "DataEntry::Copy: inputData is null" );
    }

    if ( this->data == nullptr || inputData->data == nullptr )
    {
        throw std::runtime_error(
            "DataEntry::Copy: data pointer is null" );
    }

    if ( this->type != inputData->type )
    {
        throw std::runtime_error(
            "DataEntry::Copy: data type mismatch for entry '" +
            this->name + "'" );
    }

    if ( this->size != inputData->size )
    {
        throw std::runtime_error(
            "DataEntry::Copy: data size mismatch for entry '" +
            this->name + "'" );
    }

    // Copy only the data value.
    // Name, type, and size belong to the existing DataEntry.
    this->data->Copy( inputData->data.get() );
}

void DataEntry::Dump( std::fstream & file )
{
    file << name << " , " << DataBaseType::GetName( type ) << " : ";
    this->data->Dump( file );
    file << "\n";
}

DataPara::DataPara()
{
}

DataPara::~DataPara()
{
    Clear();
}

void DataPara::UpdateDataPointer( std::unique_ptr<DataEntry> data )
{
    if ( data == nullptr )
    {
        return;
    }

    auto it = dataMap.find( data->name );

    if ( it == dataMap.end() )
    {
        // No entry with the same name exists.
        // DataPara takes ownership of the new DataEntry.
        const std::string name = data->name;
        dataMap[ name ] = std::move( data );
        return;
    }

    // Copy() validates type and size before updating the value.
    // Temporary DataEntry is destroyed automatically when unique_ptr goes out of scope.
    it->second->Copy( data.get() );
}

DataEntry * DataPara::GetDataPointer( const std::string & name )
{
    auto it = dataMap.find( name );
    if ( it != dataMap.end() )
    {
        return it->second.get();
    }
    return nullptr;
}

const DataEntry * DataPara::GetDataPointer( const std::string & name ) const
{
    auto it = dataMap.find( name );
    if ( it != dataMap.end() )
    {
        return it->second.get();
    }
    return nullptr;
}

void DataPara::DeleteDataPointer( const std::string & name )
{
    dataMap.erase( name );
}

void DataPara::Clear()
{
    dataMap.clear();
}

void DataPara::DumpData( std::fstream & file )
{
    std::cout << " Dumping database:\n";
    int count = 0;
    for ( auto & pair : dataMap )
    {
        file << ++ count << ": ";
        pair.second->Dump( file );
    }
}

EndNameSpace
