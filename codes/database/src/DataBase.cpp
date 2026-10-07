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
#include <memory>
#include "Fatal.h"
#include "DataPara.h"
#include "DataObject.h"
#include "DataField.h"

BeginNameSpace( ONEFLOW )

std::unique_ptr<DataBase> globalDataBase;

DataBase * GetGlobalDataBase()
{
    return globalDataBase.get();
}

class HXInitGlobalDataBase
{
public:
    HXInitGlobalDataBase()
    {
        globalDataBase = std::make_unique<DataBase>();
    };
    ~HXInitGlobalDataBase()
    {
        globalDataBase.reset();
    }
};

HXInitGlobalDataBase initGlobalDataBase;


DataBase::DataBase()
    : dataPara( std::make_unique<DataPara>() )
    , dataField( std::make_unique<DataField>() )
{
}

DataBase::~DataBase()
{
}

void HXWriteVoid( DataBook * dataBook, const DataEntry * dataEntry )
{
    dataEntry->data->Write( dataBook );
}

void HXWriteDataEntry( DataBook * dataBook, const DataEntry * dataEntry )
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
    dataEntry->data = CreateDataObject( dataEntry->type, dataEntry->size );
    dataEntry->data->Read( dataBook, dataEntry->size );
}

void ProcessData( const std::string & name, const std::string * value, int type, int size )
{
    auto dataEntry = std::make_unique<DataEntry>();
    dataEntry->name = name;
    dataEntry->type = type;
    dataEntry->size = size;
    if ( type == ONEFLOW::HX_STRING )
    {
        auto stringObject = std::make_unique<TDataObject< std::string > >( size );
        stringObject->CopyValue( value, size );
        dataEntry->data = std::move( stringObject );
    }
    else if ( type == HX_INT )
    {
        auto intObject = std::make_unique<TDataObject< int > >( size );
        intObject->AssignFromString( value, size );
        dataEntry->data = std::move( intObject );
    }
    else if ( type == HX_REAL )
    {
        auto realObject = std::make_unique<TDataObject< Real > >( size );
        realObject->AssignFromString( value, size );
        dataEntry->data = std::move( realObject );
    }
    else
    {
        Fatal( " Parameter Type Error \n" );
    }
    DataBase * dataBase = ONEFLOW::GetGlobalDataBase();
    dataBase->GetDataPara()->UpdateDataPointer( std::move( dataEntry ) );
}

std::unique_ptr<DataObject> CreateDataObject( int type, int size )
{
    if ( type == ONEFLOW::HX_STRING )
    {
        return std::make_unique<TDataObject< std::string > >( size );
    }
    else if ( type == HX_INT )
    {
        return std::make_unique<TDataObject< int > >( size );
    }
    else if ( type == HX_REAL )
    {
        return std::make_unique<TDataObject< Real > >( size );
    }
    else
    {
        Fatal( "Parameter Type Error In CreateDataObject" );
    }
    return nullptr;
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

const PointerWrap * GetPointerWrap( const DataField * dataField, const std::string & dataObjectName )
{
    const FieldEntry * fieldEntry = dataField->GetFieldEntry( dataObjectName );
    if ( fieldEntry == nullptr )
    {
        return nullptr;
    }
    return fieldEntry->GetPointerWrap();
}

void CreateFieldPointer( DataBase * database, std::unique_ptr<PointerWrap> pointerWrap, const std::string & dataObjectName )
{
    if ( database == nullptr )
    {
        throw std::runtime_error( "DataBase: database is not initialized" );
    }

    auto fieldEntry = std::make_unique<FieldEntry>(
        dataObjectName, std::move( pointerWrap ) );
    database->GetDataField()->UpdateFieldEntry( std::move( fieldEntry ) );
}

void * GetFieldPointerVoid( DataBase * database, const std::string & dataObjectName )
{
    if ( database == nullptr )
    {
        throw std::runtime_error( "DataBase: database is not initialized" );
    }

    PointerWrap * pointerWrap = GetPointerWrap( database->dataField.get(), dataObjectName );
    if ( pointerWrap )
    {
        return pointerWrap->GetPointer();
    }
    return nullptr;
}

const void * GetFieldPointerVoid( const DataBase * database, const std::string & dataObjectName )
{
    if ( database == nullptr )
    {
        throw std::runtime_error( "DataBase: database is not initialized" );
    }

    const PointerWrap * pointerWrap = GetPointerWrap( database->dataField.get(), dataObjectName );
    if ( pointerWrap )
    {
        return pointerWrap->GetPointer();
    }
    return nullptr;
}

void DumpDataBase( std::fstream & file )
{
    DataBase * dataBase = ONEFLOW::GetGlobalDataBase();
    dataBase->GetDataPara()->DumpData( file );
}

EndNameSpace
