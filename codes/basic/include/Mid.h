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
#include "HXVector.h"

BeginNameSpace( ONEFLOW )

template< typename T >
class Mid
{
public:
    int size, id;
    HXVector< T > data;
public:
    Mid();
    Mid( int size, int id = 0 );
    Mid( const Mid & rhs );
    Mid & operator = ( const Mid & rhs );
    ~Mid( void  );
public:
    bool operator < ( const Mid & rhs ) const;
};

template < typename T >
Mid<T>::Mid()
{
    this->size = 0;
    this->id   = 0;
}

template < typename T >
Mid<T>::Mid( int size, int id )
{
    this->size = size;
    this->id   = id;
    this->data.resize( size );
}

template < typename T >
Mid<T>::Mid( const Mid<T> & rhs )
{
    this->size = rhs.size;
    this->id   = rhs.id;
    this->data.resize( size );

    for ( int i = 0; i < size; ++ i )
    {
        this->data[ i ] = rhs.data[ i ];
    }
}

template < typename T >
Mid<T> & Mid<T>::operator = ( const Mid<T> & rhs )
{
    if ( this == & rhs ) return * this;

    if ( this->size != rhs.size )
    {
        this->size = rhs.size;
        this->data.resize( this->size );
    }

    this->id   = rhs.id;
    for ( int i = 0; i < size; ++ i )
    {
        this->data[ i ] = rhs.data[ i ];
    }

    return * this;
}

template < typename T >
Mid<T>::~Mid()
{
}

template < typename T >
bool Mid<T>::operator < ( const Mid<T> & rhs ) const
{
    if ( this->size != rhs.size ) return this->size < rhs.size;
    for ( int i = 0; i < this->size; ++ i )
    {
        if ( this->data[ i ] != rhs.data[ i ] )
        {
            return this->data[ i ] < rhs.data[ i ];
        }
    }
    return false;
}

EndNameSpace