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
#include "Constant.h"
#include <vector>
#include <cmath>
using namespace std;

BeginNameSpace( ONEFLOW )

template < typename T >
bool NotANumber( const T & value );

template < typename T >
inline T MIN( const T & value1, const T & value2 );

inline double MIN( const double & value1, const float & value2 );

inline double MIN( const float & value1, const double & value2 );

template < typename T >
inline T MIN( const T & value1, const T & value2, const T & value3 );

template < typename T >
inline T MAX( const T & a, const T & b );

inline double MAX( const double & a, const float & b );

inline double MAX( const float & a, const double & b );

template < typename T >
inline T MAX( const T & a, const T & b, const T & c );

template < typename T >
inline T MAX( const T & v1, const T & v2, const T & v3, const T & v4 );

template < typename T >
inline T MAX( const T & v1, const T & v2, const T & v3, const T & v4, const T & v5 );

template < typename T >
inline T MAX( const T & v1, const T & v2, const T & v3, const T & v4, const T & v5, const T & v6 );

template < typename T >
inline T ABS( const T & a );

template < typename T1, typename T2 >
inline T2 SIGN( const T1 & a, const T2 & b );

template < typename T >
inline T SQR( const T & a );

template < typename T >
inline T SQR( const T & a, const T & b );

template < typename T >
inline T SQR( const T & a, const T & b, const T & c );

template < typename T >
inline T SQR( const T & a, const T & b, const T & c, const T & d );

template < typename T >
inline T POWER3( const T & a );

template < typename T >
inline T POWER3( const T & a, const T & b );

template < typename T >
inline T POWER3( const T & a, const T & b, const T & c );

template < typename T >
inline T POWER4( const T & a );

template < typename T >
inline T DIST( const T & a, const T & b );

template < typename T >
inline T DIST( const T & a, const T & b, const T & c );

template < typename T >
inline T DIST( const T & a, const T & b, const T & c, const T & d );

template < typename T >
inline void SWAP( T & a, T & b );

template < typename T >
inline T COUNT( const T & ist, const T & ied );

template < typename T >
inline T COUNT( const T & ist, const T & ied, const T & jst, const T & jed );

template < typename T >
inline T COUNT( const T & ist, const T & ied, const T & jst, const T & jed, const T & kst, const T & ked );

template < typename T >
inline T COUNT( const T & ist, const T & ied, const T & jst, const T & jed, const T & kst, const T & ked, const T & lst, const T & led );


template < typename T >
bool NotANumber( const T & value )
{
    return value != value;
}

template < typename T >
inline T MIN( const T & value1, const T & value2 )
{
    return ( value1 < value2 ) ? value1 : value2;
}

inline double MIN( const double & value1, const float & value2 )
{
    return static_cast< double > ( value1 < value2 ) ? value1 : value2;
}

inline double MIN( const float & value1, const double & value2 )
{
    return static_cast< double > ( value1 < value2 ) ? value1 : value2;
}

template < typename T >
inline T MIN( const T & value1, const T & value2, const T & value3 )
{
    return MIN( MIN( value1, value2 ), value3 );
}

template < typename T >
inline T MAX( const T & a, const T & b )
{
    return ( a > b ) ? a : b;
}

inline double MAX( const double & a, const float & b )
{
    return ( a > b ) ? a : b;
}

inline double MAX( const float & a, const double & b )
{
    return ( a > b ) ? a : b;
}

template < typename T >
inline T MAX( const T & a, const T & b, const T & c )
{
    return MAX( MAX( a, b ), c );
}

template < typename T >
inline T MAX( const T & v1, const T & v2, const T & v3, const T & v4 )
{
    return MAX( MAX( v1, v2, v3 ), v4 );
}

template < typename T >
inline T MAX( const T & v1, const T & v2, const T & v3, const T & v4, const T & v5 )
{
    return MAX( MAX( v1, v2, v3, v4 ), v5 );
}

template < typename T >
inline T MAX( const T & v1, const T & v2, const T & v3, const T & v4, const T & v5, const T & v6 )
{
    return MAX( MAX( v1, v2, v3, v4, v5 ), v6 );
}

template < typename T >
inline T ABS( const T & a )
{
    return ( a < 0 ) ? -a : a;
}

template < typename T1, typename T2 >
inline T2 SIGN( const T1 & a, const T2 & b )
{
    return a * ABS( b ) / ( b + SMALL );
}

template < typename T >
inline T SQR( const T & a )
{
    return a * a;
}

template < typename T >
inline T SQR( const T & a, const T & b )
{
    return a * a + b * b;
}

template < typename T >
inline T SQR( const T & a, const T & b, const T & c )
{
    return a * a + b * b + c * c;
}

template < typename T >
inline T SQR( const T & a, const T & b, const T & c, const T & d )
{
    return a * a + b * b + c * c + d * d;
}

template < typename T >
inline T POWER3( const T & a )
{
    return a * a * a;
}

template < typename T >
inline T POWER3( const T & a, const T & b )
{
    return a * a * a + b * b * b;
}

template < typename T >
inline T POWER3( const T & a, const T & b, const T & c )
{
    return a * a * a + b * b * b + c * c * c;
}

template < typename T >
inline T POWER4( const T & a )
{
    return a * a * a * a;
}

template < typename T >
inline T DIST( const T & a, const T & b )
{
    return sqrt( SQR( a, b ) );
}

template < typename T >
inline T DIST( const T & a, const T & b, const T & c )
{
    return sqrt( SQR( a, b, c ) );
}

template < typename T >
inline T DIST( const T & a, const T & b, const T & c, const T & d )
{
    return sqrt( SQR( a, b, c, d ) );
}

template < typename T >
inline void SWAP( T & a, T & b )
{
    T c = a;
    a = b;
    b = c;
}

template < typename T >
inline T COUNT( const T & ist, const T & ied )
{
    return MAX( ( ied - ist + 1 ), 1 );
}

template < typename T >
inline T COUNT( const T & ist, const T & ied, const T & jst, const T & jed )
{
    return COUNT( ist, ied ) * COUNT( jst, jed );
}

template < typename T >
inline T COUNT( const T & ist, const T & ied, const T & jst, const T & jed, const T & kst, const T & ked )
{
    return COUNT( ist, ied ) * COUNT( jst, jed ) * COUNT( kst, ked );
}

template < typename T >
inline T COUNT( const T & ist, const T & ied, const T & jst, const T & jed, const T & kst, const T & ked, const T & lst, const T & led )
{
    return COUNT( ist, ied ) * COUNT( jst, jed ) * COUNT( kst, ked ) * COUNT( lst, led );
}

EndNameSpace