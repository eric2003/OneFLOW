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

#pragma once
#include "NamespaceMacros.h"
#include "DataBook.h"
#include "DataPara.h"
#include "DataField.h"
#include "DataObject.h"
#include "DataPointer.h"
#include "DataBaseType.h"
#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <map>

BeginNameSpace( ONEFLOW )

class DataObject;
class DataEntry;
class DataField;

class DataBase
{
public:
    DataBase();
    ~DataBase();
public:
    DataPara *dataPara;
    DataField *dataField;
};
void HXReadDataEntry( DataBook * dataBook, DataEntry * dataEntry );
void HXWriteDataEntry( DataBook * dataBook, DataEntry * dataEntry );
void HXWriteVoid( DataBook * dataBook, DataEntry * dataEntry );
void HXReadVoid( DataBook * dataBook, DataEntry * dataEntry );

DataBase * GetGlobalDataBase();
void ProcessData( const std::string & name, std::string * value, int type, int size );
DataObject * CreateDataObject( int type, int size );

class DataBase;
template < typename T >
void SetData( const std::string & name, T * value, int type, int size );

template < typename T >
T GetDataValue( const std::string & varName, DataBase * database = ONEFLOW::GetGlobalDataBase() );

//Read the value of the variable with parameter type T and name Varname from the database
template < typename T >
T GetDataValue( const std::string & varName, DataBase * database )
{
    DataEntry * dataEntry = database->dataPara->GetDataPointer( varName );

    if (dataEntry != nullptr )
    {
        DataObject * data = dataEntry->data;
        return GetDataValue< T >(data);
    }
    else
    {
        // Short-term: throw instead of exit, so unit tests can catch it
        throw std::runtime_error( "DataBase: cannot find variable \"" + varName + "\"" );
    }   
}

template < typename T >
void SetData( const std::string & name, T * value, int type, int size )
{
    DataEntry * dataEntry = new DataEntry();
    dataEntry->name = name;
    dataEntry->type = type;
    dataEntry->size = size;
    TDataObject< T > * o = new TDataObject< T >( size );
    o->CopyValue( value );
    dataEntry->data = o;

    DataBase * dataBase = ONEFLOW::GetGlobalDataBase();
    dataBase->dataPara->UpdateDataPointer( dataEntry );
}

void SetDataInt( const std::string & varName, const int & value );
void SetDataReal( const std::string & varName, const Real & value );
void SetDataString( const std::string & varName, const std::string & value );

template < typename T >
T * GetDataPointer( const std::string & varName )
{
    DataBase * database = ONEFLOW::GetGlobalDataBase();
    DataEntry * dataEntry = database->dataPara->GetDataPointer( varName );
    DataObject * data = dataEntry->data;
    return static_cast< T * >( data->GetVoidPointer() );
}

class PointerWrap;
PointerWrap * GetPointerWrap( DataField * dataField, const std::string & dataObjectName );

void * GetFieldPointerVoid( DataBase * database, const std::string & dataObjectName );

template < typename T >
T * GetFieldPointer( DataBase * database, const std::string & dataObjectName );
template < typename T, typename TStorage >
T * GetFieldPointer( TStorage * storage, const std::string & dataObjectName );
template < typename T >
T & GetFieldReference( DataBase * database, const std::string & dataObjectName );
template < typename T, typename TStorage >
T & GetFieldReference( TStorage * storage, const std::string & dataObjectName );

void CreateFieldPointer( DataBase * database, PointerWrap * pointerWrap, const std::string & dataObjectName );
template < typename TStorage >
void CreateFieldPointer( TStorage * storage, PointerWrap * pointerWrap, const std::string & dataObjectName );

template < typename T >
T * GetFieldPointer( DataBase * database, const std::string & dataObjectName )
{
    void * p = GetFieldPointerVoid( database, dataObjectName );
    if ( p )
    {
        T * pointer = reinterpret_cast< T * >( p );
        return pointer;
    }
    return 0;
}

template < typename T, typename TStorage >
T * GetFieldPointer( TStorage * storage, const std::string & dataObjectName )
{
    DataBase * database = storage->GetDataBase();
    T * pointer = ONEFLOW::GetFieldPointer< T >( database, dataObjectName );
    return pointer;
}

template < typename T >
T & GetFieldReference( DataBase * database, const std::string & dataObjectName )
{
    return * ONEFLOW::GetFieldPointer< T >( database, dataObjectName );
}

template < typename T, typename TStorage >
T & GetFieldReference( TStorage * storage, const std::string & dataObjectName )
{
    return * ONEFLOW::GetFieldPointer< T, TStorage >( storage, dataObjectName );
}

template < typename TStorage >
void CreateFieldPointer( TStorage * storage, PointerWrap * pointerWrap, const std::string & dataObjectName )
{
    DataBase * database = storage->GetDataBase();
    ONEFLOW::CreateFieldPointer( database, pointerWrap, dataObjectName );
}

void DumpDataBase( std::fstream & file );

EndNameSpace
