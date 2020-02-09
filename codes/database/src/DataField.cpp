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


#include "DataField.h"
#include "DataObject.h"
#include "DataPointer.h"

BeginNameSpace( ONEFLOW )

DataF::DataF()
{
    this->name = "";
    this->data = 0;
}

DataF::DataF( const string & name, PointerWrap * data )
{
    this->name = name;
    this->data = data;
}

DataF::~DataF()
{
    delete data;
}

DataField::DataField()
{
    dataSet = new DataSET;
}

DataField::~DataField()
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

void DataField::UpdateDataF( DataF * dataf )
{
    DataF * findData = this->GetDataF( dataf->name );
    if ( ! findData )
    {
        dataSet->insert( dataf );
    }
    else
    {
        if ( findData != dataf )
        {
            delete dataf;
        }
    }
}

DataF * DataField::GetDataF( const string & name )
{
    DataF * data = new DataF( name, 0 );
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

void DataField::DeleteDataF( const string & name )
{
    DataF * data = new DataF( name, 0 );
    DataSET::iterator iter = dataSet->find( data );
    if ( iter != dataSet->end() )
    {
        delete ( * iter );
        dataSet->erase( iter );
    }
    delete data;
}

EndNameSpace
