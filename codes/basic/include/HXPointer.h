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

template < typename T >
class HXPointer
{
public:
    HXPointer();
    HXPointer( UInt nSize );
    ~HXPointer();
protected:
    HXVector< T * > data;
    bool del_flag;
public:
    void SetDeleteFlag( bool del_flag );
    T *& operator[] ( UInt i );
    HXPointer & operator= ( HXPointer &rhs );
public:
    void resize( UInt nSize );
    UInt size();
    void push_back( T * value );
};

template < typename T >
HXPointer<T>::HXPointer()
{
    del_flag = false;
}

template < typename T >
HXPointer<T>::HXPointer( UInt nSize ) :
    data( nSize )
{
    del_flag = false;
}

template < typename T >
HXPointer<T>::~HXPointer()
{
    if ( del_flag )
    {
        DeletePointer( data );
    }
}

template < typename T >
void HXPointer<T>::SetDeleteFlag( bool del_flag )
{
    this->del_flag = del_flag;
}

template < typename T >
T *& HXPointer<T>::operator[] ( UInt i )
{
    return this->data[ i ];
}

template < typename T >
HXPointer<T> & HXPointer<T>::operator= ( HXPointer<T> &rhs )
{
    if ( this == & rhs ) return *this;
    this->data = rhs.data;
    return *this;
}

template < typename T >
void HXPointer<T>::resize( UInt nSize )
{
    this->data.resize( nSize );
}

template < typename T >
UInt HXPointer<T>::size()
{
    return this->data.size();
}

template < typename T >
void HXPointer<T>::push_back( T * value )
{
    this->data.push_back( value );
}

template < typename T >
void CreatePointer( HXVector< T * > & pointer, int nSize );
template < typename T >
void DeletePointer( HXVector< T * > & pointer );

template < typename T >
void CreatePointer( HXVector< T * > & pointer, int nSize )
{
    pointer.resize( nSize );
    for ( int i = 0; i < nSize; ++ i )
    {
        pointer[ i ] = new T();
    }
}

template < typename T >
void DeletePointer( HXVector< T * > & pointer )
{
    for ( UInt i = 0; i < pointer.size(); ++ i )
    {
        delete pointer[ i ];
    }
    pointer.resize( 0 );
}

EndNameSpace