/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
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
#include <fstream>
#include <string>

BeginNameSpace( ONEFLOW )

class DataPage
{
public:
    using CharMemory = std::vector<char>;
public:
    DataPage();
    ~DataPage();
    // A DataPage can hold up to `maxUnitSize` bytes (default ~1GB).
    // Implicit copy would silently deep-copy that buffer on any accidental
    // pass-by-value, assignment, or vector<DataPage> reallocation -- a
    // severe, invisible performance trap. Disable copy explicitly to force
    // callers to be deliberate (e.g. wrap in unique_ptr, as DataBook does).
    DataPage( const DataPage & )            = delete;
    DataPage & operator=( const DataPage & ) = delete;

    // Move is cheap: it's just a pointer/size swap inside std::vector<char>,
    // O(1) regardless of buffer size. Must be explicitly defaulted here,
    // because declaring the deleted copy ctor above already counts as a
    // "user-declared copy constructor", which suppresses implicit move
    // generation just like a user-declared destructor would.
    DataPage( DataPage && )            = default;
    DataPage & operator=( DataPage && ) = default;
public:
    HXSize_t size() const;
    void Read ( void * data, HXSize_t dataSize );
    void Read( void * data, HXSize_t position, HXSize_t dataSize ) const;
    void Write( const void * data, HXSize_t dataSize );
    void Write( const void * data, HXSize_t position, HXSize_t dataSize );
    void ReadFile ( std::fstream & file );
    void WriteFile( std::fstream & file );
    void ToString( std::string & str );

    char * data();
    const char * data() const;
    char * CurrentPtr();

    void MoveToBegin() { MoveToPosition( 0 ); };
    void MoveToEnd  () { currPos = size(); };
    void ReSize( HXSize_t newSize );
    void Send( int pId, int tag );
    void Recv( int pId, int tag );
    void Bcast( int rootid );
protected:
    void MoveToPosition( HXSize_t position );
    void Advance( HXOffset_t offset );
protected:
    HXSize_t currPos;
    CharMemory dataMemory;     // value member, automatically managed
public:
    char * PtrAt( int offset );
};

EndNameSpace
