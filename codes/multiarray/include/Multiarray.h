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
#include "Range.h"
#include "Memop.h"

BeginNameSpace( ONEFLOW )

template < typename T, int N >
class Multiarray;

template < typename T, int N >
class Multiarray
{
public:
    typedef Multiarray< T, N > TArray;
public:
    ArrayPointer< T > arrayPointer;
public:
    ~Multiarray()
    {
    }

    Multiarray( const Range & r0 )
    {
        arrayPointer.Allocate( r0 );
    }

    Multiarray( const Range & r0, const Range & r1 )
    {
        arrayPointer.Allocate( r0, r1 );
    }

    Multiarray( const Range & r0, const Range & r1, const Range & r2 )
    {
        arrayPointer.Allocate( r0, r1, r2 );
    }

    Multiarray( const Range & r0, const Range & r1, const Range & r2, const Range & r3 )
    {
        arrayPointer.Allocate( r0, r1, r2, r3 );
    }

    Multiarray( T * dataPointer, const Range & r0 )
    {
        arrayPointer.Allocate( dataPointer, r0 );
    }

    Multiarray( T * dataPointer, const Range & r0, const Range & r1 )
    {
        arrayPointer.Allocate( dataPointer, r0, r1 );
    }

    Multiarray( T * dataPointer, const Range & r0, const Range & r1, const Range & r2 )
    {
        arrayPointer.Allocate( dataPointer, r0, r1, r2 );
    }

    Multiarray( T * dataPointer, const Range & r0, const Range & r1, const Range & r2, const Range & r3 )
    {
        arrayPointer.Allocate( dataPointer, r0, r1, r2, r3 );
    }
public:
    const T & operator()( int i0 ) const
    { 
        return arrayPointer.datap1[ i0 ];
    }

    T & operator()( int i0 )
    {
        return arrayPointer.datap1[ i0 ];
    }
    
    const T & operator()( int i0, int i1 ) const
    { 
        return arrayPointer.datap2[ i1 ][ i0 ];
    }

    T & operator()( int i0, int i1 )
    {
        return arrayPointer.datap2[ i1 ][ i0 ];
    }

    const T & operator()( int i0, int i1, int i2 ) const
    {
        return arrayPointer.datap3[ i2 ][ i1 ][ i0 ];
    }

    T & operator()( int i0, int i1, int i2 ) 
    {
        return arrayPointer.datap3[ i2 ][ i1 ][ i0 ];
    }

    const T & operator()( int i0, int i1, int i2, int i3 ) const
    {
        return arrayPointer.datap4[ i3 ][ i2 ][ i1 ][ i0 ];
    }

    T & operator()( int i0, int i1, int i2, int i3 ) 
    {
        return arrayPointer.datap4[ i3 ][ i2 ][ i1 ][ i0 ];
    }
};


EndNameSpace