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
#include "Word.h"
#include "DataBook.h"
#include "DataBaseIO.h"
#include "LogFile.h"
#include <string>
#include <set>
using namespace std;

BeginNameSpace( ONEFLOW )

class DataObject
{
public:
    DataObject() {};
    virtual ~DataObject() {};
public:
    virtual void * GetVoidPointer() { return 0; };
    virtual void Write( DataBook * dataBook ) {};
    virtual void Read( DataBook * dataBook, int numberOfElements ) {};
    virtual void Copy( DataObject * dataObject ) {};
};

template < typename T >
T GetDataValue( DataObject * dataObject, int iElement = 0 );

template < typename T >
T GetDataValue( DataObject * dataObject, int iElement )
{
    T * data = static_cast< T *>( dataObject->GetVoidPointer() );
    return data[ iElement ];
}

template < typename T >
class TDataObject : public DataObject
{
public:
    TDataObject( int nSize )
    {
        this->data.resize( nSize );
    }
    ~TDataObject(){}
public:
    vector< T > data;
public:
    void * GetVoidPointer() { return & data[ 0 ]; };
    void CopyValue( string * valueIn )
    {
        UInt nSize = this->data.size();
        for ( UInt i = 0; i < nSize; ++ i )
        {
            data[ i ] = StringToDigit< T >( valueIn[ i ], std::dec );
        }
    }

    void CopyValue( T * valueIn )
    {
        UInt size = this->data.size();
        for ( UInt i = 0; i < size; ++ i )
        {
            data[ i ] = valueIn[ i ];
        }
    }

    void Write( DataBook * dataBook )
    {
        UInt numberOfElements = this->data.size();
        for ( UInt iElement = 0; iElement < numberOfElements; ++ iElement )
        {
            T & value = this->data[ iElement ];
            ONEFLOW::HXWrite( dataBook, value );
        }
    }

    void Read( DataBook * dataBook, int numberOfElements )
    {
        this->data.resize( numberOfElements );
        for ( int iElement = 0; iElement < numberOfElements; ++ iElement )
        {
            ONEFLOW::HXRead( dataBook, this->data[ iElement ] );
        }
    }

    void Copy( DataObject * dataObject )
    {
        UInt numberOfElements = this->data.size();
        for ( UInt iElement = 0; iElement < numberOfElements; ++ iElement )
        {
            TDataObject< T > * tDataObject = static_cast<TDataObject< T > *>( dataObject );
            data[ iElement ] = tDataObject->data[ iElement ];
        }
    }
};


template <>
class TDataObject< string > : public DataObject
{
public:
    TDataObject( int nSize )
    {
        this->data.resize( nSize );
    }
    ~TDataObject() {}
public:
    vector< string > data;
public:
    void * GetVoidPointer() { return & data[ 0 ]; };
    void CopyValue( string * valueIn )
    {
        UInt nSize = this->data.size();
        for ( UInt i = 0; i < nSize; ++ i )
        {
            data[ i ] = valueIn[ i ];
        }
    }

    void Write( DataBook * dataBook )
    {
        UInt numberOfElements = this->data.size();
        for ( UInt iElement = 0; iElement < numberOfElements; ++ iElement )
        {
            string & value = this->data[ iElement ];
            ONEFLOW::HXWrite( dataBook, value );
        }
    }

    void Read( DataBook * dataBook, int numberOfElements )
    {
        this->data.resize( numberOfElements );
        for ( int iElement = 0; iElement < numberOfElements; ++ iElement )
        {
            ONEFLOW::HXRead( dataBook, this->data[ iElement ] );
        }
    }

    void Copy( DataObject * dataObject )
    {
        UInt numberOfElements = this->data.size();
        for ( UInt iElement = 0; iElement < numberOfElements; ++ iElement )
        {
            TDataObject< string > * tDataObject = static_cast<TDataObject< string > *>( dataObject );
            data[ iElement ] = tDataObject->data[ iElement ];
        }
    }
};

EndNameSpace
