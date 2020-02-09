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

#include "DataBase.h"
#include "Stop.h"
#include "DataPara.h"
#include "DataObject.h"
#include "DataField.h"

BeginNameSpace( ONEFLOW )

DataBase * globalDataBase = 0;

DataBase * GetGlobalDataBase()
{
    return globalDataBase;
}

class HXInitGlobalDataBase
{
public:
    HXInitGlobalDataBase()
    {
        globalDataBase = new DataBase();
    };
    ~HXInitGlobalDataBase()
    {
        delete globalDataBase;
    }
};

HXInitGlobalDataBase initGlobalDataBase;


DataBase::DataBase()
{
    dataPara = new DataPara();
    dataField = new DataField();
}

DataBase::~DataBase()
{
    delete dataPara;
    delete dataField;
}

void HXWriteVoid( DataBook * dataBook, DataV * datav )
{
    datav->data->Write( dataBook );
}

void HXWriteDataV( DataBook * dataBook, DataV * datav )
{
    ONEFLOW::HXWrite( dataBook, datav->name );
    ONEFLOW::HXWrite( dataBook, datav->type );
    ONEFLOW::HXWrite( dataBook, datav->size );
    ONEFLOW::HXWriteVoid( dataBook, datav );
}

void HXReadDataV( DataBook * dataBook, DataV * datav )
{
    ONEFLOW::HXRead( dataBook, datav->name );
    ONEFLOW::HXRead( dataBook, datav->type );
    ONEFLOW::HXRead( dataBook, datav->size );
    ONEFLOW::HXReadVoid( dataBook, datav );
}

void HXReadVoid( DataBook * dataBook, DataV * datav )
{
    DataObject * o = CreateDataObject( datav->type, datav->size );
    datav->data = o;
    datav->data->Read( dataBook, datav->size );
}

void ProcessData( const string & name, string * value, int type, int size )
{
    DataV * datav = new DataV();
    datav->name = name;
    datav->type = type;
    datav->size = size;
    if ( type == ONEFLOW::HX_STRING )
    {
        TDataObject< string > * stringObject = new TDataObject< string >( size );
        stringObject->CopyValue( value );
        datav->data = stringObject;
    }
    else if ( type == HX_INT )
    {
        TDataObject< int > * intObject = new TDataObject< int >( size );
        intObject->CopyValue( value );
        datav->data = intObject;
    }
    else if ( type == HX_REAL )
    {
        TDataObject< Real > * realObject = new TDataObject< Real >( size );
        realObject->CopyValue( value );
        datav->data = realObject;
    }
    else
    {
        Stop( " Parameter Type Error \n" );
    }
    DataBase * dataBase = ONEFLOW::GetGlobalDataBase();
    dataBase->dataPara->UpdateDataPointer( datav );
}

DataObject * CreateDataObject( int type, int size )
{
    if ( type == ONEFLOW::HX_STRING )
    {
        TDataObject< string > * stringObject = new TDataObject< string >( size );
        return stringObject;
    }
    else if ( type == HX_INT )
    {
        TDataObject< int > * intObject = new TDataObject< int >( size );
        return intObject;
    }
    else if ( type == HX_REAL )
    {
        TDataObject< Real > * realObject = new TDataObject< Real >( size );
        return realObject;
    }
    else
    {
        Stop( " Parameter Type Error In CreateDataObject\n" );
        return 0;
    }
}

void SetDataInt( const std::string & varName, int & value )
{
    SetData( varName, & value, HX_INT, 1 );
}

void SetDataReal( const std::string & varName, Real & value )
{
    SetData( varName, & value, HX_REAL, 1 );
}

void SetDataString( const std::string & varName, Real & value )
{
    SetData( varName, & value, HX_STRING, 1 );
}

PointerWrap * GetPointerWrap( DataField * dataField, const string & dataObjectName )
{
    DataF * dataf = dataField->GetDataF( dataObjectName );
    return dataf->GetPointerWrap();
}

void CreateFieldPointer( DataBase * database, PointerWrap * pointerWrap, const string & dataObjectName )
{
    DataF * dataf = new DataF( dataObjectName, pointerWrap );
    database->dataField->UpdateDataF( dataf );
}

void * GetFieldPointerVoid( DataBase * database, const string & dataObjectName )
{
    PointerWrap * pointerWrap = GetPointerWrap( database->dataField, dataObjectName );
    if ( pointerWrap )
    {
        void * p = pointerWrap->GetPointer();
        return p;
    }
    return 0;
}


EndNameSpace
