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
#include <memory>
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
#include <stdexcept>

BeginNameSpace( ONEFLOW )

class DataObject;
class DataEntry;
class DataField;

class DataBase
{
public:
    DataBase();
    ~DataBase();
private:
    std::unique_ptr<DataPara> dataPara;
    std::unique_ptr<DataField> dataField;

public:
    // The database owns these stores for its entire lifetime.
    // Use pointer access only for compatibility with existing nullable APIs.
    DataPara * GetDataPara() { return dataPara.get(); }
    const DataPara * GetDataPara() const { return dataPara.get(); }

    DataPara & RequireDataPara() { return *dataPara; }
    const DataPara & RequireDataPara() const { return *dataPara; }

    DataField * GetDataField() { return dataField.get(); }
    const DataField * GetDataField() const { return dataField.get(); }

    DataField & RequireDataField() { return *dataField; }
    const DataField & RequireDataField() const { return *dataField; }
};
std::unique_ptr<DataEntry> HXReadDataEntry( DataBook * dataBook );
void HXWriteDataEntry( DataBook * dataBook, const DataEntry * dataEntry );
void HXWriteVoid( DataBook * dataBook, const DataEntry * dataEntry );

DataBase * GetGlobalDataBase();
void ProcessData( const std::string & name, const std::string * value, int type, int size );
std::unique_ptr<DataObject> CreateDataObject( int type, int size );

class DataBase;
template < typename T >
void SetData( const std::string & name, T * value, int type, int size );

template < typename T >
T GetDataValue( const std::string & varName, DataBase * database = ONEFLOW::GetGlobalDataBase() );
template < typename T >
T GetDataValue( const std::string & varName, const DataBase * database );

//Read the value of the variable with parameter type T and name Varname from the database
template < typename T >
T GetDataValue( const std::string & varName, DataBase * database )
{
    if ( database == nullptr )
    {
        throw std::runtime_error( "DataBase: database is not initialized" );
    }

    DataEntry * dataEntry = database->RequireDataPara().FindDataEntry( varName );

    if (dataEntry != nullptr )
    {
        DataObject & data = dataEntry->GetDataObject();
        return GetDataValue< T >( &data );
    }
    else
    {
        // Short-term: throw instead of exit, so unit tests can catch it
        throw std::runtime_error( "DataBase: cannot find variable \"" + varName + "\"" );
    }
}

template < typename T >
T GetDataValue( const std::string & varName, const DataBase * database )
{
    if ( database == nullptr )
    {
        throw std::runtime_error( "DataBase: database is not initialized" );
    }

    const DataEntry * dataEntry = database->RequireDataPara().FindDataEntry( varName );

    if ( dataEntry != nullptr )
    {
        const DataObject & data = dataEntry->GetDataObject();
        return GetDataValue< T >( &data );
    }

    throw std::runtime_error( "DataBase: cannot find variable \"" + varName + "\"" );
}

template < typename T >
void SetData( const std::string & name, T * value, int type, int size )
{
    auto o = std::make_unique<TDataObject< T > >( size );
    o->CopyValue( value, size );
    auto dataEntry = std::make_unique<DataEntry>(
        name, type, size, std::move( o ) );

    DataBase * dataBase = ONEFLOW::GetGlobalDataBase();
    dataBase->RequireDataPara().SetDataEntry( std::move( dataEntry ) );
}

void SetDataInt( const std::string & varName, const int & value );
void SetDataReal( const std::string & varName, const Real & value );
void SetDataString( const std::string & varName, const std::string & value );

template < typename T >
T * GetDataPointer( const std::string & varName )
{
    DataBase * database = ONEFLOW::GetGlobalDataBase();
    DataEntry * dataEntry = database->RequireDataPara().FindDataEntry( varName );

    // Required lookup: match GetDataValue -- missing name must not
    // dereference a null DataEntry.
    if ( dataEntry == nullptr )
    {
        // Short-term: throw instead of exit, so unit tests can catch it
        throw std::runtime_error(
            "DataBase: cannot find variable \"" + varName + "\"" );
    }

    DataObject & data = dataEntry->GetDataObject();
    return static_cast< T * >( data.GetVoidPointer() );
}

class PointerWrap;
PointerWrap * GetPointerWrap( DataField * dataField, const std::string & dataObjectName );
const PointerWrap * GetPointerWrap( const DataField * dataField, const std::string & dataObjectName );

// Field storage lookup (optional): returns nullptr if the named field
// is not registered. The DataBase itself is required and must be initialized.
void * GetFieldPointerVoid( DataBase * database, const std::string & dataObjectName );
const void * GetFieldPointerVoid( const DataBase * database, const std::string & dataObjectName );

template < typename T >
T * GetFieldPointer( DataBase * database, const std::string & dataObjectName );
template < typename T >
const T * GetFieldPointer( const DataBase * database, const std::string & dataObjectName );
template < typename T, typename TStorage >
T * GetFieldPointer( TStorage * storage, const std::string & dataObjectName );
template < typename T, typename TStorage >
const T * GetFieldPointer( const TStorage * storage, const std::string & dataObjectName );

// Required field access: throws when the named field is not registered.
// The DataBase itself is required and must be initialized.
template < typename T >
T & GetFieldReference( DataBase * database, const std::string & dataObjectName );
template < typename T >
const T & GetFieldReference( const DataBase * database, const std::string & dataObjectName );
template < typename T, typename TStorage >
T & GetFieldReference( TStorage * storage, const std::string & dataObjectName );
template < typename T, typename TStorage >
const T & GetFieldReference( const TStorage * storage, const std::string & dataObjectName );

void CreateFieldPointer( DataBase * database, std::unique_ptr<PointerWrap> pointerWrap, const std::string & dataObjectName );
template < typename TStorage >
void CreateFieldPointer( TStorage * storage, std::unique_ptr<PointerWrap> pointerWrap, const std::string & dataObjectName );

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
    DataBase & database = storage->RequireDataBase();
    T * pointer = ONEFLOW::GetFieldPointer< T >( &database, dataObjectName );
    return pointer;
}

template < typename T >
const T * GetFieldPointer( const DataBase * database, const std::string & dataObjectName )
{
    const void * p = GetFieldPointerVoid( database, dataObjectName );
    if ( p )
    {
        return reinterpret_cast< const T * >( p );
    }
    return nullptr;
}

template < typename T, typename TStorage >
const T * GetFieldPointer( const TStorage * storage, const std::string & dataObjectName )
{
    const DataBase * database = storage->GetDataBase();
    return ONEFLOW::GetFieldPointer< T >( database, dataObjectName );
}

template < typename T >
T & GetFieldReference( DataBase * database, const std::string & dataObjectName )
{
    T * pointer = ONEFLOW::GetFieldPointer< T >( database, dataObjectName );
    if ( pointer == nullptr )
    {
        throw std::runtime_error(
            "DataBase: cannot find field \"" + dataObjectName + "\"" );
    }
    return * pointer;
}


template < typename T >
const T & GetFieldReference( const DataBase * database, const std::string & dataObjectName )
{
    const T * pointer = ONEFLOW::GetFieldPointer< T >( database, dataObjectName );
    if ( pointer == nullptr )
    {
        throw std::runtime_error(
            "DataBase: cannot find field \"" + dataObjectName + "\"" );
    }
    return * pointer;
}

template < typename T, typename TStorage >
T & GetFieldReference( TStorage * storage, const std::string & dataObjectName )
{
    return ONEFLOW::GetFieldReference< T >( &storage->RequireDataBase(), dataObjectName );
}

template < typename T, typename TStorage >
const T & GetFieldReference( const TStorage * storage, const std::string & dataObjectName )
{
    return ONEFLOW::GetFieldReference< T >( &storage->RequireDataBase(), dataObjectName );
}

template < typename TStorage >
void CreateFieldPointer( TStorage * storage, std::unique_ptr<PointerWrap> pointerWrap, const std::string & dataObjectName )
{
    DataBase & database = storage->RequireDataBase();
    ONEFLOW::CreateFieldPointer( &database, std::move( pointerWrap ), dataObjectName );
}

void DumpDataBase( std::fstream & file );

EndNameSpace
