/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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
class HXMid
{
public:
    std::size_t size, id;
    HXVector< T > data;
public:
    HXMid();
    HXMid( std::size_t size, std::size_t id = 0 );
    HXMid( const HXMid & rhs );
    HXMid & operator = ( const HXMid & rhs );
    ~HXMid( void  );
public:
    bool operator < ( const HXMid & rhs ) const;
};

template < typename T >
HXMid<T>::HXMid()
{
    this->size = 0;
    this->id   = 0;
}

template < typename T >
HXMid<T>::HXMid( std::size_t size, std::size_t id )
{
    this->size = size;
    this->id   = id;
    this->data.resize( size );
}

template < typename T >
HXMid<T>::HXMid( const HXMid<T> & rhs )
{
    this->size = rhs.size;
    this->id   = rhs.id;
    this->data.resize( size );

    for ( std::size_t i = 0; i < size; ++ i )
    {
        this->data[ i ] = rhs.data[ i ];
    }
}

template < typename T >
HXMid<T> & HXMid<T>::operator = ( const HXMid<T> & rhs )
{
    if ( this == & rhs ) return * this;

    if ( this->size != rhs.size )
    {
        this->size = rhs.size;
        this->data.resize( this->size );
    }

    this->id   = rhs.id;
    for ( std::size_t i = 0; i < size; ++ i )
    {
        this->data[ i ] = rhs.data[ i ];
    }

    return * this;
}

template < typename T >
HXMid<T>::~HXMid()
{
}

template < typename T >
bool HXMid<T>::operator < ( const HXMid<T> & rhs ) const
{
    if ( this->size != rhs.size ) return this->size < rhs.size;
    for ( std::size_t i = 0; i < this->size; ++ i )
    {
        if ( this->data[ i ] != rhs.data[ i ] )
        {
            return this->data[ i ] < rhs.data[ i ];
        }
    }
    return false;
}

EndNameSpace
