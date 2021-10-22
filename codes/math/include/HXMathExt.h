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
#include "HXDefine.h"
#include "HXMath.h"
#include <vector>
using namespace std;

BeginNameSpace( ONEFLOW )

template < typename T >
T SUM( vector< T > & a );

template < typename T >
T MaxField( vector< T > & field );

template < typename T >
T MinField( vector< T > & field );

template < typename T >
inline T SUM( vector< T > & a )
{
    T sum = 0;
    int nElements = a.size();
    for ( UInt iElement = 0; iElement < nElements; ++ iElement )
    {
        sum += a[ iElement ];
    }
    return sum;
}

template < typename T >
inline T MaxField( vector< T > & field )
{
    T maxValue = field[ 0 ];

    UInt nElements = field.size();
    for ( UInt iElement = 1; iElement < nElements; ++ iElement )
    {
        maxValue = ONEFLOW::MAX( maxValue, field[ iElement ] );
    }
    return maxValue;
}

template < typename T >
inline T MinField( vector< T > & field )
{
    T minValue = field[ 0 ];

    UInt nElements = field.size();
    for ( UInt iElement = 1; iElement < nElements; ++ iElement )
    {
        minValue = ONEFLOW::MIN( minValue, field[ iElement ] );
    }
    return minValue;
}

EndNameSpace
