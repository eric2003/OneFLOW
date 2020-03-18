/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2020 He Xin and the OneFLOW contributors.
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

BeginNameSpace( ONEFLOW )

DataV::DataV()
{
    this->name = "";
    this->data = 0;
}

DataV::DataV( const string & name, int type, int size, DataObject * data )
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

DataPara::DataPara()
{
    dataSet = new DataSET;
}

DataPara::~DataPara()
{
    DataSET::iterator iter;
    for ( iter = dataSet->begin(); iter != dataSet->end(); ++ iter )
    {
        DataObject * dataObject = reinterpret_cast< DataObject * > ( ( * iter )->data );
        delete dataObject;
    }

    dataSet->clear();

    delete dataSet;
}

void DataPara::UpdateDataPointer( DataV * data )
{
    DataV * findData = this->GetDataPointer( data->name );
    if ( findData )
    {
        findData->Copy( data );
        delete data;
        return;
    }

    dataSet->insert( data );
}

DataV * DataPara::GetDataPointer( const string & name )
{
    DataV * data = new DataV( name, 0, 0, 0 );
    DataSET::iterator iter = dataSet->find( data );
    delete data;
    if ( iter != dataSet->end() )
    {
        return ( * iter );
    }
    else
    {
        return 0;
    }
}

void DataPara::DeleteDataPointer( const string & name )
{
    DataV * data = new DataV( name, 0, 0, 0 );
    DataSET::iterator iter = dataSet->find( data );
    if ( iter != dataSet->end() )
    {
        delete ( * iter );
        dataSet->erase( iter );
    }
    delete data;
}

EndNameSpace
