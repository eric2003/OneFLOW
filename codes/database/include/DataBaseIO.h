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

#include "HXDefine.h"
#include "HXArray.h"
#include "DataBook.h"

#include <cstddef>
#include <fstream>
#include <string>
#include <vector>

BeginNameSpace( ONEFLOW )

class DataBook;

// ===========================================================================
// Forward declarations
//
// Only the "helper" templates that are called by other templates defined
// later in this header need a forward declaration. Leaf overloads do not,
// because a template definition is itself a declaration.
// ===========================================================================

template < typename TIO, typename T >
void HXReadVector( TIO * tio, std::vector< T > & field );

template < typename T >
void HXReadVector( std::fstream * file, std::vector< T > & field );

template < typename T >
void HXReadVector( std::fstream * file, HXVector< T > & field );

template < typename TIO, typename T >
void HXWriteVector( TIO * tio, const std::vector< T > & field );

template < typename T >
void HXWriteVector( std::fstream * file, const std::vector< T > & field );

template < typename vect2D >
void HXReadVector2D( DataBook * dataBook, vect2D & field2D );

template < typename vect2D >
void HXWriteVector2D( DataBook * dataBook, const vect2D & field2D );

template < typename vect2D >
void HXAppendVector2D( DataBook * dataBook, vect2D & field2D );

// ===========================================================================
// Non-template overloads
// ===========================================================================

void HXRead( DataBook * dataBook, std::string & cs );
void HXWrite( DataBook * dataBook, const std::string & cs );

void HXRead( DataBook * dataBook, MRField * field );
void HXWrite( DataBook * dataBook, MRField * field );

// ===========================================================================
// HXRead
// ===========================================================================

// Read a single value from a generic IO object.
template < typename TIO, typename T >
void HXRead( TIO * tio, T & value )
{
    tio->Read( reinterpret_cast< char * >( & value ), sizeof( T ) );
}

// Read a single value from a binary file stream.
template < typename T >
void HXRead( std::fstream * file, T & value )
{
    file->read( reinterpret_cast< char * >( & value ), sizeof( T ) );
}

// Read an array of nElement values from a generic IO object.
template < typename TIO, typename T >
void HXRead( TIO * tio, T * field, int nElement )
{
    if ( nElement <= 0 ) return;
    tio->Read( field, nElement * sizeof( T ) );
}

// Read an array of nElement values from a binary file stream.
template < typename T >
void HXRead( std::fstream * file, T * field, int nElement )
{
    if ( nElement <= 0 ) return;
    file->read( reinterpret_cast< char * >( field ), nElement * sizeof( T ) );
}

// Read the content of a std::vector from a generic IO object.
// The vector must already have the correct size.
template < typename TIO, typename T >
void HXReadVector( TIO * tio, std::vector< T > & field )
{
    std::size_t nElement = field.size();
    if ( nElement == 0 ) return;
    tio->Read( field.data(), nElement * sizeof( T ) );
}

// Read the content of a std::vector from a generic IO object.
template < typename TIO, typename T >
void HXRead( TIO * tio, std::vector< T > & field )
{
    HXReadVector( tio, field );
}

// Read the content of a HXVector from a generic IO object.
template < typename TIO, typename T >
void HXRead( TIO * tio, HXVector< T > & field )
{
    HXReadVector( tio, field );
}

// Read the content of a std::vector from a binary file stream.
template < typename T >
void HXReadVector( std::fstream * file, std::vector< T > & field )
{
    int nElement = static_cast< int >( field.size() );
    HXRead( file, field.data(), nElement );
}

// Read the content of a std::vector from a binary file stream.
template < typename T >
void HXRead( std::fstream * file, std::vector< T > & field )
{
    HXReadVector( file, field );
}

// Read the content of a HXVector from a binary file stream.
// NOTE: HXVector must expose data() and size() compatible with std::vector.
template < typename T >
void HXReadVector( std::fstream * file, HXVector< T > & field )
{
    int nElement = static_cast< int >( field.size() );
    HXRead( file, field.data(), nElement );
}

// Read the content of a HXVector from a binary file stream.
template < typename T >
void HXRead( std::fstream * file, HXVector< T > & field )
{
    HXReadVector( file, field );
}

// ===========================================================================
// HXWrite
// ===========================================================================

// Write a single value to a generic IO object.
template < typename TIO, typename T >
void HXWrite( TIO * tio, const T & value )
{
    tio->Write( reinterpret_cast< const char * >( & value ), sizeof( T ) );
}

// Write a single value to a binary file stream.
template < typename T >
void HXWrite( std::fstream * file, const T & value )
{
    file->write( reinterpret_cast< const char * >( & value ), sizeof( T ) );
}

// Write an array of nElement values to a generic IO object.
template < typename TIO, typename T >
void HXWrite( TIO * tio, const T * field, int nElement )
{
    if ( nElement <= 0 ) return;
    tio->Write( field, nElement * sizeof( T ) );
}

// Write an array of nElement values to a binary file stream.
template < typename T >
void HXWrite( std::fstream * file, const T * field, int nElement )
{
    if ( nElement <= 0 ) return;
    file->write( reinterpret_cast< const char * >( field ), nElement * sizeof( T ) );
}

// Write the content of a std::vector to a generic IO object.
template < typename TIO, typename T >
void HXWriteVector( TIO * tio, const std::vector< T > & field )
{
    int nElement = static_cast< int >( field.size() );
    if ( nElement <= 0 ) return;
    tio->Write( field.data(), nElement * sizeof( T ) );
}

// Write the content of a std::vector to a generic IO object.
template < typename TIO, typename T >
void HXWrite( TIO * tio, const std::vector< T > & field )
{
    HXWriteVector( tio, field );
}

// Write the content of a HXVector to a generic IO object.
template < typename TIO, typename T >
void HXWrite( TIO * tio, const HXVector< T > & field )
{
    HXWriteVector( tio, field );
}

// Write the content of a std::vector to a binary file stream.
template < typename T >
void HXWriteVector( std::fstream * file, const std::vector< T > & field )
{
    int nElement = static_cast< int >( field.size() );
    HXWrite( file, field.data(), nElement );
}

// Write the content of a std::vector to a binary file stream.
template < typename T >
void HXWrite( std::fstream * file, const std::vector< T > & field )
{
    HXWriteVector( file, field );
}

// Write the content of a HXVector to a binary file stream.
template < typename T >
void HXWrite( std::fstream * file, const HXVector< T > & field )
{
    HXWriteVector( file, field );
}

// ===========================================================================
// 2D read / write (vector of vectors, HXVector of HXVector)
// ===========================================================================

// Read a 2D field: for each outer element, read its size, resize, then read.
template < typename vect2D >
void HXReadVector2D( DataBook * dataBook, vect2D & field2D )
{
    HXSize_t nElem = field2D.size();
    if ( nElem == 0 ) return;

    for ( HXSize_t iElem = 0; iElem < nElem; ++ iElem )
    {
        auto & field = field2D[ iElem ];

        int nSubElem = 0;
        HXRead( dataBook, nSubElem );

        field.resize( nSubElem );
        HXRead( dataBook, field );
    }
}

// Read a std::vector< std::vector< T > >.
template < typename T >
void HXRead( DataBook * dataBook, std::vector< std::vector< T > > & field2D )
{
    HXReadVector2D( dataBook, field2D );
}

// Read a HXVector< HXVector< T > >.
template < typename T >
void HXRead( DataBook * dataBook, HXVector< HXVector< T > > & field2D )
{
    HXReadVector2D( dataBook, field2D );
}

// Write a 2D field: for each outer element, write its size, then the data.
template < typename vect2D >
void HXWriteVector2D( DataBook * dataBook, const vect2D & field2D )
{
    HXSize_t nElem = field2D.size();
    if ( nElem == 0 ) return;

    for ( HXSize_t iElem = 0; iElem < nElem; ++ iElem )
    {
        const auto & field = field2D[ iElem ];

        int nSubElem = static_cast< int >( field.size() );
        HXWrite( dataBook, nSubElem );

        HXWrite( dataBook, field );
    }
}

// Write a std::vector< std::vector< T > >.
template < typename T >
void HXWrite( DataBook * dataBook, const std::vector< std::vector< T > > & field2D )
{
    HXWriteVector2D( dataBook, field2D );
}

// Write a HXVector< HXVector< T > >.
template < typename T >
void HXWrite( DataBook * dataBook, const HXVector< HXVector< T > > & field2D )
{
    HXWriteVector2D( dataBook, field2D );
}

// ===========================================================================
// HXAppend
// ===========================================================================

// Append a single value to a DataBook.
template < typename T >
void HXAppend( DataBook * dataBook, const T & value )
{
    dataBook->Append( & value, sizeof( T ) );
}

// Append an array of nElement values to a DataBook.
template < typename T >
void HXAppend( DataBook * dataBook, const T * field, int nElement )
{
    if ( nElement <= 0 ) return;
    dataBook->Append( field, nElement * sizeof( T ) );
}

// Append the content of a std::vector to a DataBook.
template < typename T >
void HXAppendVector( DataBook * dataBook, const std::vector< T > & field )
{
    HXSize_t nElement = field.size();
    if ( nElement == 0 ) return;
    dataBook->Append( field.data(), nElement * sizeof( T ) );
}

// Append the content of a std::vector to a DataBook.
template < typename T >
void HXAppend( DataBook * dataBook, const std::vector< T > & field )
{
    HXAppendVector( dataBook, field );
}

// Append the content of a HXVector to a DataBook.
template < typename T >
void HXAppend( DataBook * dataBook, const HXVector< T > & field )
{
    HXAppendVector( dataBook, field );
}

// Append a 2D field: for each outer element, append its size, then the data.
template < typename vect2D >
void HXAppendVector2D( DataBook * dataBook, const vect2D & field2D )
{
    HXSize_t nElem = field2D.size();
    if ( nElem == 0 ) return;

    for ( HXSize_t iElem = 0; iElem < nElem; ++ iElem )
    {
        const auto & field = field2D[ iElem ];

        HXSize_t nSubElem = field.size();
        HXAppend( dataBook, nSubElem );

        HXAppend( dataBook, field );
    }
}

// Append a std::vector< std::vector< T > >.
template < typename T >
void HXAppend( DataBook * dataBook, const std::vector< std::vector< T > > & field2D )
{
    HXAppendVector2D( dataBook, field2D );
}

// Append a HXVector< HXVector< T > >.
template < typename T >
void HXAppend( DataBook * dataBook, const HXVector< HXVector< T > > & field2D )
{
    HXAppendVector2D( dataBook, field2D );
}

EndNameSpace
