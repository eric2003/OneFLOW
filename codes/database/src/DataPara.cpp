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
#include "DataObject.h"
#include "DataBaseType.h"
#include <iostream>

BeginNameSpace( ONEFLOW )

DataV::DataV()
{
    this->name = "";
    this->data = nullptr;
}

DataV::DataV( const std::string & name, int type, int size, DataObject * data )
{
    this->name = name;
    this->type = type;
    this->size = size;
    this->data = data;
}

DataV::~DataV()
{
    delete data;
}

void DataV::Copy( DataV * inputData )
{
    this->data->Copy( inputData->data );
}

void DataV::Dump( std::fstream & file )
{
    file << name << " , " << DataBaseType::GetName( type ) << " : ";
    this->data->Dump( file );
    file << "\n";
}

DataPara::DataPara()
{
    dataMap = new DataMap;
}

DataPara::~DataPara()
{
    for ( auto & pair : *dataMap )
    {
        delete pair.second;     // DataV destructor deletes the DataObject
    }
    dataMap->clear();
    delete dataMap;
}

void DataPara::UpdateDataPointer( DataV * data )
{
    auto it = dataMap->find( data->name );
    if ( it != dataMap->end() )
    {
        // Already exists ¡ú copy content and discard the new object
        it->second->Copy( data );
        delete data;
        return;
    }
    // New key ¡ú take ownership
    ( *dataMap )[ data->name ] = data;
}

DataV * DataPara::GetDataPointer( const std::string & name )
{
    auto it = dataMap->find( name );
    if ( it != dataMap->end() )
    {
        return it->second;
    }
    return nullptr;
}

void DataPara::DeleteDataPointer( const std::string & name )
{
    auto it = dataMap->find( name );
    if ( it != dataMap->end() )
    {
        delete it->second;
        dataMap->erase( it );
    }
}

void DataPara::DumpData( std::fstream & file )
{
    std::cout << " Dumping database:\n";
    int count = 0;
    for ( auto & pair : *dataMap )
    {
        file << ++count << ": ";
        pair.second->Dump( file );
    }
}

EndNameSpace
