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
#include "ScalarFieldRecord.h"
#include "DataStorage.h"
#include "DataBase.h"
#include <stdexcept>
#include <vector>

BeginNameSpace( ONEFLOW )

std::map< std::string, int > GFieldDim::data;

GFieldDim::GFieldDim()
{
}

GFieldDim::~GFieldDim()
{
}

void GFieldDim::AddField( const std::string & fileName, int nEqu )
{
    GFieldDim::data[ fileName ] = nEqu;
}

int GFieldDim::GetNEqu( const std::string & fileName )
{
    std::map< std::string, int >::iterator iter;
    iter = GFieldDim::data.find( fileName );
    if ( iter != GFieldDim::data.end() )
    {
        return iter->second;
    }
    return -1;
}


ScalarFieldRecord::ScalarFieldRecord()
{
}

ScalarFieldRecord::~ScalarFieldRecord()
{
}

void ScalarFieldRecord::AddField( MRField * field, int nEqu )
{
    if ( field == nullptr )
    {
        throw std::invalid_argument( "ScalarFieldRecord::AddField: field must not be null" );
    }

    this->nEquList.push_back( nEqu );
    this->fields.push_back( field );
}

MRField * ScalarFieldRecord::GetField( int id )
{
    if ( id < 0 || static_cast< size_t >( id ) >= this->fields.size() )
    {
        throw std::out_of_range( "ScalarFieldRecord::GetField: field index is out of range" );
    }

    return this->fields[ id ];
}

void ScalarFieldRecord::AddFieldRecord( DataStorage * dataStorage, StringField & fieldNameList )
{
    if ( dataStorage == nullptr )
    {
        throw std::invalid_argument( "ScalarFieldRecord::AddFieldRecord: data storage must not be null" );
    }

    // Resolve all fields before mutating the record, so a missing field cannot
    // leave a partially populated collection of non-owning pointers.
    std::vector< MRField * > resolvedFields;
    resolvedFields.reserve( fieldNameList.size() );

    for ( const std::string & fieldName : fieldNameList )
    {
        MRField * field = ONEFLOW::GetFieldPointer< MRField >( dataStorage, fieldName );
        if ( field == nullptr )
        {
            throw std::runtime_error(
                "ScalarFieldRecord::AddFieldRecord: field '" + fieldName + "' was not found in data storage" );
        }
        resolvedFields.push_back( field );
    }

    for ( size_t iField = 0; iField < fieldNameList.size(); ++ iField )
    {
        const int nEqu = GFieldDim::GetNEqu( fieldNameList[ iField ] );
        this->AddField( resolvedFields[ iField ], nEqu );
    }
}


EndNameSpace
