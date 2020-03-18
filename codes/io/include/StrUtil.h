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
#include "Configure.h"
#include <string>
#include <sstream>
using namespace std;

BeginNameSpace( ONEFLOW )

template < typename T >
T & GetReference( T * x )
{
    return * x;
}

template < typename T1, typename T2 >
string AddString( const T1 & v1, const T2 & v2 )
{
    ostringstream oss;
    oss << v1 << v2;
    return oss.str();
}

template < typename T1, typename T2, typename T3 >
string AddString( const T1 & v1, const T2 & v2, const T3 & v3 )
{
    ostringstream oss;
    oss << v1 << v2 << v3;
    return oss.str();
}

template < typename T1, typename T2, typename T3, typename T4 >
string AddString( const T1 & v1, const T2 & v2, const T3 & v3, const T4 & v4 )
{
    ostringstream oss;
    oss << v1 << v2 << v3 << v4;
    return oss.str();
}

template < typename T1, typename T2, typename T3, typename T4, typename T5>
string AddString( const T1 & v1, const T2 & v2, const T3 & v3, const T4 & v4, const T5 & v5 )
{
    ostringstream oss;
    oss << v1 << v2 << v3 << v4 << v5;
    return oss.str();
}

template < typename T1, typename T2, typename T3, typename T4, typename T5, typename T6 >
string AddString( const T1 & v1, const T2 & v2, const T3 & v3, const T4 & v4, const T5 & v5, const T6 & v6 )
{
    ostringstream oss;
    oss << v1 << v2 << v3 << v4 << v5 << v6;
    return oss.str();
}


EndNameSpace
