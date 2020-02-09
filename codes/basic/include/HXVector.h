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
#include "HXType.h"
#include <vector>
using namespace std;

BeginNameSpace( ONEFLOW )

template < typename T >
class HXVector : public vector< T >
{
public:
    HXVector(){};
    ~HXVector(){};
    HXVector( const UInt count )
        : vector< T >( count )
    {
        ;
    }
    HXVector( const UInt count, const T& value )
        : vector< T >( count, value )
    {
        ;
    }
    HXVector( T * first, T * last ) :
        vector< T >( first, last )
    {
        ;
    }
public:
    HXVector< T >& operator =( const T& value )
    {
        for ( UInt i = 0; i < this->size(); ++ i )
        {
            ( *this )[ i ] = value;
        }
        return *this;
    }
};

template < typename T >
void AllocateVector( HXVector< HXVector< T > > & data, int ni, int nj )
{
    if ( nj <= 0 ) return;
    data.resize( ni );
    for ( int i = 0; i < ni; ++ i )
    {
        data[ i ].resize( nj );
    }
}

template < typename T >
void AllocateVector( HXVector< HXVector< HXVector< T > > > & data, int ni, int nj, int nk )
{
    data.resize( ni );
    for ( int i = 0; i < ni; ++ i )
    {
        data[ i ].resize( nj );
        for ( int j = 0; j < nj; ++ j )
        {
            data[ i ][ j ].resize( nk );
        }
    }
}

EndNameSpace