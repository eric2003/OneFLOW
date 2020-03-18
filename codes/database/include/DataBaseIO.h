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
#include "HXDefine.h"
#include "HXArray.h"
#include "DataBook.h"
#include <string>
#include <fstream>
using namespace std;

BeginNameSpace( ONEFLOW )

class DataBook;

template < typename TIO, typename T >
void HXRead( TIO * tio, T & value );

template < typename T >
void HXRead( fstream * file, T & value );

template < typename TIO, typename T >
void HXRead( TIO * tio, T * field, int nElement );

template < typename T >
void HXRead( fstream * file, T * field, int nElement );

template < typename TIO, typename T >
void HXRead( TIO * tio, HXVector< T > & field );

template < typename T >
void HXRead( fstream * file, HXVector< T > & field );

template < typename TIO, typename T >
void HXWrite( TIO * tio, T & value );

template < typename T >
void HXWrite( fstream * file, T & value );

template < typename TIO, typename T >
void HXWrite( TIO * tio, T * field, int nElement );

template < typename T >
void HXWrite( fstream * file, T * field, int nElement );

template < typename TIO, typename T >
void HXWrite( TIO * tio, HXVector< T > & field );

template < typename T >
void HXWrite( fstream * file, HXVector< T > & field );

template < typename T >
void HXAppend( DataBook * dataBook, T & value );

template < typename T >
void HXAppend( DataBook * dataBook, T * field, int nElement );

template < typename T >
void HXAppend( DataBook * dataBook, HXVector< T > & field );

void HXRead( DataBook * dataBook, string & cs );
void HXWrite( DataBook * dataBook, string & cs );

void HXRead( DataBook * dataBook, MRField * field );
void HXWrite( DataBook * dataBook, MRField * field );

template < typename T >
void HXRead( DataBook * dataBook, HXVector< HXVector< T > > & field2D );

template < typename T >
void HXWrite( DataBook * dataBook, HXVector< HXVector< T > > & field2D );

template < typename T >
void HXAppend( DataBook * dataBook, HXVector< HXVector< T > > & field2D );

template < typename TIO, typename T >
void HXRead( TIO * tio, T & value )
{
    tio->Read( reinterpret_cast< char * >( & value ), sizeof( T ) );
}

template < typename T >
void HXRead( fstream * file, T & value )
{
    file->read( reinterpret_cast< char * >( & value ), sizeof( T ) );
}


template < typename TIO, typename T >
void HXRead( TIO * tio, T * field, int nElement )
{
    if ( nElement <= 0 ) return;
    tio->Read( field, nElement * sizeof( T ) );
}

template < typename T >
void HXRead( fstream * file, T * field, int nElement )
{
    if ( nElement <= 0 ) return;
    file->read( reinterpret_cast< char * >( field ), nElement * sizeof( T ) );
}

template < typename TIO, typename T >
void HXRead( TIO * tio, HXVector< T > & field )
{
    size_t nElement = field.size();
    if ( nElement <= 0 ) return;
    tio->Read( & field[ 0 ], nElement * sizeof( T ) );
}

template < typename T >
void HXRead( fstream * file, HXVector< T > & field )
{
    int nElement = field.size();
    HXRead( file, & field[ 0 ], nElement );
}


template < typename TIO, typename T >
void HXWrite( TIO * tio, T & value )
{
    tio->Write( reinterpret_cast< char * >( & value ), sizeof( T ) );
}

template < typename T >
void HXWrite( fstream * file, T & value )
{
    file->write( reinterpret_cast< char * >( & value ), sizeof( T ) );
}

template < typename TIO, typename T >
void HXWrite( TIO * tio, T * field, int nElement )
{
    if ( nElement <= 0 ) return;
    tio->Write( field, nElement * sizeof( T ) );
}

template < typename T >
void HXWrite( fstream * file, T * field, int nElement )
{
    if ( nElement <= 0 ) return;
    file->write( reinterpret_cast< char * >( field ), nElement * sizeof( T ) );
}

template < typename TIO, typename T >
void HXWrite( TIO * tio, HXVector< T > & field )
{
    int nElement = field.size();
    if ( nElement <= 0 ) return;
    tio->Write( & field[ 0 ], nElement * sizeof( T ) );
}

template < typename T >
void HXWrite( fstream * file, HXVector< T > & field )
{
    UInt nElement = static_cast<int> (field.size());
    HXWrite( file, & field[ 0 ], nElement );
}


template < typename T >
void HXRead( DataBook * dataBook, HXVector< HXVector< T > > & field2D )
{
    UInt nElem = field2D.size();
    if ( nElem == 0 ) return;
    for ( UInt iElem = 0; iElem < nElem; ++ iElem )
    {
        HXVector< T > & field = field2D[ iElem ];

        int nSubElem = 0;
        HXRead( dataBook, nSubElem );

        field.resize( nSubElem );
        HXRead( dataBook, field );
    }
}

template < typename T >
void HXWrite( DataBook * dataBook, HXVector< HXVector< T > > & field2D )
{
    UInt nElem = field2D.size();
    if ( nElem == 0 ) return;
    for ( UInt iElem = 0; iElem < nElem; ++ iElem )
    {
        HXVector< T > & field = field2D[ iElem ];

        int nSubElem = field.size();
        HXWrite( dataBook, nSubElem );

        HXWrite( dataBook, field );
    }
}

template < typename T >
void HXAppend( DataBook * dataBook, T & value )
{
    dataBook->Append( & value, sizeof( T ) );
}

template < typename T >
void HXAppend( DataBook * dataBook, T * field, int nElement )
{
    if ( nElement <= 0 ) return;
    dataBook->Append( field, nElement * sizeof( T ) );
}

template < typename T >
void HXAppend( DataBook * dataBook, HXVector< T > & field )
{
    UInt nElement = field.size();
    if ( nElement <= 0 ) return;
    dataBook->Append( & field[ 0 ], nElement * sizeof( T ) );
}

template < typename T >
void HXAppend( DataBook * dataBook, HXVector< HXVector< T > > & field2D )
{
    UInt nElem = field2D.size();
    if ( nElem == 0 ) return;
    for ( UInt iElem = 0; iElem < nElem; ++ iElem )
    {
        HXVector< T > & field = field2D[ iElem ];

        UInt nSubElem = field.size();
        HXAppend( dataBook, nSubElem );

        HXAppend( dataBook, field );
    }
}

EndNameSpace
