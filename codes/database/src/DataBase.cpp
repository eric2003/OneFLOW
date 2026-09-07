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

#include "DataBase.h"
#include "Fatal.h"
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

void HXWriteVoid( DataBook * dataBook, DataEntry * dataEntry )
{
    dataEntry->data->Write( dataBook );
}

void HXWriteDataEntry( DataBook * dataBook, DataEntry * dataEntry )
{
    ONEFLOW::HXWrite( dataBook, dataEntry->name );
    ONEFLOW::HXWrite( dataBook, dataEntry->type );
    ONEFLOW::HXWrite( dataBook, dataEntry->size );
    ONEFLOW::HXWriteVoid( dataBook, dataEntry );
}

void HXReadDataEntry( DataBook * dataBook, DataEntry * dataEntry )
{
    ONEFLOW::HXRead( dataBook, dataEntry->name );
    ONEFLOW::HXRead( dataBook, dataEntry->type );
    ONEFLOW::HXRead( dataBook, dataEntry->size );
    ONEFLOW::HXReadVoid( dataBook, dataEntry );
}

void HXReadVoid( DataBook * dataBook, DataEntry * dataEntry )
{
    DataObject * o = CreateDataObject( dataEntry->type, dataEntry->size );
    dataEntry->data = o;
    dataEntry->data->Read( dataBook, dataEntry->size );
}

void ProcessData( const std::string & name, std::string * value, int type, int size )
{
    DataEntry * dataEntry = new DataEntry();
    dataEntry->name = name;
    dataEntry->type = type;
    dataEntry->size = size;
    if ( type == ONEFLOW::HX_STRING )
    {
        TDataObject< std::string > * stringObject = new TDataObject< std::string >( size );
        stringObject->CopyValue( value );
        dataEntry->data = stringObject;
    }
    else if ( type == HX_INT )
    {
        TDataObject< int > * intObject = new TDataObject< int >( size );
        intObject->CopyValue( value );
        dataEntry->data = intObject;
    }
    else if ( type == HX_REAL )
    {
        TDataObject< Real > * realObject = new TDataObject< Real >( size );
        realObject->CopyValue( value );
        dataEntry->data = realObject;
    }
    else
    {
        Fatal( " Parameter Type Error \n" );
    }
    DataBase * dataBase = ONEFLOW::GetGlobalDataBase();
    dataBase->dataPara->UpdateDataPointer( dataEntry );
}

DataObject * CreateDataObject( int type, int size )
{
    if ( type == ONEFLOW::HX_STRING )
    {
        TDataObject< std::string > * stringObject = new TDataObject< std::string >( size );
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
        Fatal( "Parameter Type Error In CreateDataObject" );
    }
}

void SetDataInt( const std::string & varName, const int & value )
{
    // Make a local copy so we can take its address safely
    int tmp = value;
    SetData( varName, & tmp, HX_INT, 1 );
}

void SetDataReal( const std::string & varName, const Real & value )
{
    // Make a local copy so we can take its address safely
    Real tmp = value;
    SetData( varName, & tmp, HX_REAL, 1 );
}

void SetDataString( const std::string & varName, const std::string & value )
{
    // Make a local copy so we can take its address safely
    std::string tmp = value;
    SetData( varName, & tmp, HX_STRING, 1 );
}

PointerWrap * GetPointerWrap( DataField * dataField, const std::string & dataObjectName )
{
    FieldEntry * fieldEntry = dataField->GetFieldEntry( dataObjectName );
    if ( fieldEntry == nullptr )
    {
        return nullptr;
    }
    return fieldEntry->GetPointerWrap();
}

void CreateFieldPointer( DataBase * database, PointerWrap * pointerWrap, const std::string & dataObjectName )
{
    FieldEntry * fieldEntry = new FieldEntry( dataObjectName, pointerWrap );
    database->dataField->UpdateFieldEntry( fieldEntry );
}

void * GetFieldPointerVoid( DataBase * database, const std::string & dataObjectName )
{
    PointerWrap * pointerWrap = GetPointerWrap( database->dataField, dataObjectName );
    if ( pointerWrap )
    {
        return pointerWrap->GetPointer();
    }
    return nullptr;
}

void DumpDataBase( std::fstream & file )
{
    DataBase * dataBase = ONEFLOW::GetGlobalDataBase();
    dataBase->dataPara->DumpData( file );
}

EndNameSpace
