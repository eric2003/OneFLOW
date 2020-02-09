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
#include "AdtTree.h"

BeginNameSpace( ONEFLOW )

template < typename T >
class Point
{
public:
    typedef Point< T > point_type;
public:
    T x, y, z;
    int id;
public:
    Point();
    Point( const T & x, const T & y, const T & z, int id = 0 );
    Point( const point_type & rhs );
    Point & operator = ( const point_type & rhs );
    ~Point();
public:
    void SetPoint( const T & x );
    void SetPoint( const T & x, const T & y );
    void SetPoint( const T & x, const T & y, const T & z );
    void SetPoint( const T & x, const T & y, const T & z, int id );
    void MoveRelatively( const T & dx, const T & dy, const T & dz );
    void MoveTo( const T & newX, const T & newY, const T & newZ );
public:
    bool operator < ( const point_type & rhs ) const;
    bool Compare( const point_type & rhs, const T & tolerance ) const;
};

EndNameSpace

#include "Point.hpp"