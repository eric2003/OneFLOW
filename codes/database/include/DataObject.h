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
#include "Word.h"
#include "DataBook.h"
#include "DataBaseIO.h"
#include <string>
#include <set>
#include <algorithm>

BeginNameSpace( ONEFLOW )

class DataObject
{
public:
    DataObject() {};
    virtual ~DataObject() {};
public:
    virtual void * GetVoidPointer() { return 0; };
    virtual const void * GetVoidPointer() const { return nullptr; };
    virtual void Write( DataBook * dataBook ) {};
    virtual void Read( DataBook * dataBook, int numberOfElements ) {};
    virtual void Copy( DataObject * dataObject ) {};
    virtual void Dump( std::fstream & file ) {};
};

template < typename T >
T GetDataValue( DataObject * dataObject, int iElement = 0 );
template < typename T >
T GetDataValue( const DataObject * dataObject, int iElement );

template < typename T >
T GetDataValue( DataObject * dataObject, int iElement )
{
    T * data = static_cast< T *>( dataObject->GetVoidPointer() );
    return data[ iElement ];
}

template < typename T >
T GetDataValue( const DataObject * dataObject, int iElement )
{
    const T * data = static_cast< const T * >( dataObject->GetVoidPointer() );
    return data[ iElement ];
}

template < typename T >
void TDataObjectDump( std::fstream &file, std::vector< T > data )
{
    if ( data.size() == 0 ) return;
    file << data[ 0 ];
    for ( int i = 1; i < static_cast<int>(data.size()); ++ i )
    {
        file << " , ";
        file << data[ i ];
    }
}

template < typename T >
class TDataObject : public DataObject
{
public:
    explicit TDataObject( int nSize )
    {
        this->data.resize( static_cast<HXSize_t>(nSize) );
    }

    // Virtual destructor for safe polymorphic deletion via base‑class pointer
    virtual ~TDataObject() override = default;

    // Delete copy constructor to avoid object slicing in inheritance hierarchy
    TDataObject(const TDataObject&) = delete;
    TDataObject& operator=(const TDataObject&) = delete;

    // Enable move semantics
    TDataObject(TDataObject&&) noexcept = default;
    TDataObject& operator=(TDataObject&&) noexcept = default;

public:
    std::vector< T > data;

public:
    // Return raw pointer to underlying storage; return nullptr if container is empty
    void* GetVoidPointer() override
    {
        if (data.empty())
            return nullptr;
        return &data[0];
    };

    const void* GetVoidPointer() const override
    {
        if (data.empty())
            return nullptr;
        return data.data();
    };

    // Copy values from typed‑array, copy at most nCopyElements items
    void AssignFromString( const std::string * valueIn, HXSize_t nCopyElements )
    {
        const HXSize_t nSize = this->data.size();
        const HXSize_t nActual = std::min(nSize, nCopyElements);
        for ( HXSize_t i = 0; i < nActual; ++ i )
        {
            data[ i ] = StringToDigit< T >( valueIn[ i ], std::dec );
        }
    }

    // Copy values from typed‑array, copy at most nCopyElements items
    void CopyValue( const T* valueIn, HXSize_t nCopyElements )
    {
        const HXSize_t size = this->data.size();
        const HXSize_t nActual = std::min(size, nCopyElements);
        for ( HXSize_t i = 0; i < nActual; ++ i )
        {
            data[ i ] = valueIn[ i ];
        }
    }

    // Serialize internal data to DataBook
    void Write( DataBook* dataBook ) override
    {
        const HXSize_t numberOfElements = this->data.size();
        for ( HXSize_t iElement = 0; iElement < numberOfElements; ++ iElement )
        {
            T& value = this->data[ iElement ];
            ONEFLOW::HXWrite( dataBook, value );
        }
    }

    // Deserialize numberOfElements items from DataBook into internal storage
    void Read( DataBook* dataBook, int numberOfElements ) override
    {
        this->data.resize( static_cast<HXSize_t>(numberOfElements) );
        for ( int iElement = 0; iElement < numberOfElements; ++ iElement )
        {
            ONEFLOW::HXRead( dataBook, this->data[ iElement ] );
        }
    }

    // Copy content from another DataObject instance
    // Use dynamic_cast for runtime type checking to prevent undefined behaviour
    void Copy( DataObject* dataObject ) override
    {
        if (dataObject == nullptr)
            return;

        TDataObject< T >* tDataObject = dynamic_cast<TDataObject< T >*>(dataObject);
        if (tDataObject == nullptr)
            return;

        const HXSize_t srcSize = tDataObject->data.size();
        const HXSize_t dstSize = this->data.size();
        const HXSize_t nCopy = std::min(srcSize, dstSize);

        for ( HXSize_t iElement = 0; iElement < nCopy; ++ iElement )
        {
            data[ iElement ] = tDataObject->data[ iElement ];
        }
    }

    // Dump data content to output file stream
    void Dump( std::fstream& file ) override
    {
        TDataObjectDump( file, data );
    }
};

EndNameSpace
