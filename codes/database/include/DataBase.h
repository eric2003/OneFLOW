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

#pragma once
#include "Configure.h"
#include "DataBook.h"
#include "DataPara.h"
#include "DataField.h"
#include "DataObject.h"
#include "DataPointer.h"
#include <iostream>
#include <string>
#include <set>
using namespace std;

BeginNameSpace( ONEFLOW )

const int HX_INT    = 1;
const int HX_FLOAT  = 2;
const int HX_DOUBLE = 3;
const int HX_REAL   = 4;
const int HX_STRING = 5;
const int HX_BOOL   = 6;

class DataObject;
class DataV;
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
void HXReadDataV( DataBook * dataBook, DataV * datav );
void HXWriteDataV( DataBook * dataBook, DataV * datav );
void HXWriteVoid( DataBook * dataBook, DataV * datav );
void HXReadVoid( DataBook * dataBook, DataV * datav );

DataBase * GetGlobalDataBase();
void ProcessData( const string & name, string * value, int type, int size );
DataObject * CreateDataObject( int type, int size );

class DataBase;
template < typename T >
void SetData( const string & name, T * value, int type, int size );

template < typename T >
T GetDataValue( const std::string & varName, DataBase * database = ONEFLOW::GetGlobalDataBase() );

//从数据库中读取参数类型为T，名称为varName的变量的值
template < typename T >
T GetDataValue( const std::string & varName, DataBase * database )
{
    DataV * datav = database->dataPara->GetDataPointer( varName );

    if (datav != NULL)
    {
        DataObject * data = datav->data;
        return GetDataValue< T >(data);
    }
    else
    {
        cerr << "can't find:" << varName << " in database!!" << endl;
        exit(EXIT_FAILURE);
    }   
}

template < typename T >
void SetData( const string & name, T * value, int type, int size )
{
    DataV * datav = new DataV();
    datav->name = name;
    datav->type = type;
    datav->size = size;
    TDataObject< T > * o = new TDataObject< T >( size );
    o->CopyValue( value );
    datav->data = o;

    DataBase * dataBase = ONEFLOW::GetGlobalDataBase();
    dataBase->dataPara->UpdateDataPointer( datav );
}

void SetDataInt( const std::string & varName, int & value );
void SetDataReal( const std::string & varName, Real & value );
void SetDataString( const std::string & varName, Real & value );

template < typename T >
T * GetDataPointer( const std::string & varName )
{
    DataBase * database = ONEFLOW::GetGlobalDataBase();
    DataV * datav = database->dataPara->GetDataPointer( varName );
    DataObject * data = datav->data;
    return static_cast< T * >( data->GetVoidPointer() );
}

class PointerWrap;
PointerWrap * GetPointerWrap( DataField * dataField, const string & dataObjectName );

void * GetFieldPointerVoid( DataBase * database, const string & dataObjectName );

template < typename T >
T * GetFieldPointer( DataBase * database, const string & dataObjectName );
template < typename T, typename TStorage >
T * GetFieldPointer( TStorage * storage, const string & dataObjectName );
template < typename T >
T & GetFieldReference( DataBase * database, const string & dataObjectName );
template < typename T, typename TStorage >
T & GetFieldReference( TStorage * storage, const string & dataObjectName );

void CreateFieldPointer( DataBase * database, PointerWrap * pointerWrap, const string & dataObjectName );
template < typename TStorage >
void CreateFieldPointer( TStorage * storage, PointerWrap * pointerWrap, const string & dataObjectName );

template < typename T >
T * GetFieldPointer( DataBase * database, const string & dataObjectName )
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
T * GetFieldPointer( TStorage * storage, const string & dataObjectName )
{
    DataBase * database = storage->GetDataBase();
    T * pointer = ONEFLOW::GetFieldPointer< T >( database, dataObjectName );
    return pointer;
}

template < typename T >
T & GetFieldReference( DataBase * database, const string & dataObjectName )
{
    return * ONEFLOW::GetFieldPointer< T >( database, dataObjectName );
}

template < typename T, typename TStorage >
T & GetFieldReference( TStorage * storage, const string & dataObjectName )
{
    return * ONEFLOW::GetFieldPointer< T, TStorage >( storage, dataObjectName );
}

template < typename TStorage >
void CreateFieldPointer( TStorage * storage, PointerWrap * pointerWrap, const string & dataObjectName )
{
    DataBase * database = storage->GetDataBase();
    ONEFLOW::CreateFieldPointer( database, pointerWrap, dataObjectName );
}

EndNameSpace
