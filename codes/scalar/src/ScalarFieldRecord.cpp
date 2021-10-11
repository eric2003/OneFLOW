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
#include "ScalarFieldRecord.h"
#include "DataStorage.h"
#include "DataBase.h"

BeginNameSpace( ONEFLOW )

std::map< string, int > GFieldDim::data;

GFieldDim::GFieldDim()
{
}

GFieldDim::~GFieldDim()
{
}

void GFieldDim::AddField( const string & fileName, int nEqu )
{
    GFieldDim::data[ fileName ] = nEqu;
}

int GFieldDim::GetNEqu( const string & fileName )
{
    std::map< string, int >::iterator iter;
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
    this->nEquList.push_back( nEqu );
    this->fields.push_back( field );
}

MRField * ScalarFieldRecord::GetField( int id )
{
    return this->fields[ id ];
}

void ScalarFieldRecord::AddFieldRecord( DataStorage * dataStorage, StringField & fieldNameList )
{
    for ( int iField = 0; iField < fieldNameList.size(); ++ iField )
    {
        string & fieldName = fieldNameList[ iField ];
        MRField * field = ONEFLOW::GetFieldPointer< MRField >( dataStorage, fieldName );
        int nEqu = GFieldDim::GetNEqu( fieldName );
        this->AddField( field, nEqu );
    }
}


EndNameSpace
