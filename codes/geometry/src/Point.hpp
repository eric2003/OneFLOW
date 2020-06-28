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
#include "HXMath.h"

BeginNameSpace( ONEFLOW )

template < typename T >
Point<T>::Point()
{
    this->x  = static_cast< T > ( 0 );
    this->y  = static_cast< T > ( 0 );
    this->z  = static_cast< T > ( 0 );
    this->id = 0;
}

template < typename T >
Point<T>::Point( const T & x, const T & y, const T & z, int id )
{
    this->x  = x;
    this->y  = y;
    this->z  = z;
    this->id = id;
}

template < typename T >
Point<T>::Point( const Point<T> & rhs )
{
    this->x  = rhs.x;
    this->y  = rhs.y;
    this->z  = rhs.z;
    this->id = rhs.id;
}

template < typename T >
Point<T>::~Point()
{
}

template < typename T >
Point<T> & Point<T>::operator = ( const Point<T> & rhs )
{
    if ( this == & rhs ) return * this;

    this->x  = rhs.x;
    this->y  = rhs.y;
    this->z  = rhs.z;
    this->id = rhs.id;

    return * this;
}

template < typename T >
void Point<T>::SetPoint( const T & x )
{
    this->x  = x;
}

template < typename T >
void Point<T>::SetPoint( const T & x, const T & y )
{
    this->x  = x;
    this->y  = y;
}

template < typename T >
void Point<T>::SetPoint( const T & x, const T & y, const T & z )
{
    this->x  = x;
    this->y  = y;
    this->z  = z;
}

template < typename T >
void Point<T>::SetPoint( const T & x, const T & y, const T & z, int id )
{
    this->x  = x;
    this->y  = y;
    this->z  = z;
    this->id = id;
}

template < typename T >
void Point<T>::MoveRelatively( const T & dx, const T & dy, const T & dz )
{
    x += dx;
    y += dy;
    z += dz;
}

template < typename T >
void Point<T>::MoveTo( const T & newX, const T & newY, const T & newZ )
{
    this->x  = newX;
    this->y  = newY;
    this->z  = newZ;
}

template < typename T >
bool Point<T>::operator < ( const Point<T> & rhs ) const
{
    T dx = x - rhs.x;
    T dy = y - rhs.y;
    T dz = z - rhs.z;

    T diff = 1.0e-10;

    if ( ABS( dx ) > diff ) return x < rhs.x;
    if ( ABS( dy ) > diff ) return y < rhs.y;
    if ( ABS( dz ) > diff ) return z < rhs.z;

    return false;
}

template < typename T >
bool Point<T>::Compare( const Point<T> & rhs, const T & tolerance ) const
{
    T dx = x - rhs.x;
    T dy = y - rhs.y;
    T dz = z - rhs.z;

    if ( ABS( dx ) > tolerance ) return x < rhs.x;
    if ( ABS( dy ) > tolerance ) return y < rhs.y;
    if ( ABS( dz ) > tolerance ) return z < rhs.z;

    return false;
}

EndNameSpace