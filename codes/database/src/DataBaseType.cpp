/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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

#include "DataBaseType.h"

BeginNameSpace( ONEFLOW )

map< int, std::string > DataBaseType::nameMap;
map< std::string, int > DataBaseType::indexMap;
bool DataBaseType::init_flag = false;

DataBaseType::DataBaseType()
{
}

DataBaseType::~DataBaseType()
{
    ;
}

void DataBaseType::Init()
{
    if ( DataBaseType::init_flag ) return;
    DataBaseType::init_flag = true;
    DataBaseType::AddItem( "int", HX_INT );
    DataBaseType::AddItem( "float", HX_FLOAT );
    DataBaseType::AddItem( "double", HX_DOUBLE );
    DataBaseType::AddItem( "Real", HX_REAL );
    DataBaseType::AddItem( "string", HX_STRING );
    DataBaseType::AddItem( "bool", HX_BOOL );
}

void DataBaseType::AddItem( const std::string &name, int index )
{
    DataBaseType::indexMap.insert( pair< std::string, int >( name, index ) );
    DataBaseType::nameMap.insert( pair< int, std::string >( index, name ) );
}

int DataBaseType::GetIndex( const std::string & name )
{
    return DataBaseType::indexMap[ name ];
}

string & DataBaseType::GetName( int index )
{
    return DataBaseType::nameMap[ index ];
}

EndNameSpace
