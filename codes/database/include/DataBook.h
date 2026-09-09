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
#include "DataPage.h"
#include "HXType.h"
#include <vector>
#include <string>
#include <fstream>
#include <memory>

BeginNameSpace( ONEFLOW )

class DataBook;

using DATA_COMPRESS = void( * )( DataBook *& dataBook );
using DATA_DECOMPRESS = void( * )( DataBook *  dataBook );

class DataBook
{
public:
    // unitSize controls how many bytes each internal DataPage holds before
    // data spills into the next page. Defaults to ~1GB for production use;
    // tests can pass a small value to exercise cross-page logic directly.
    explicit DataBook( HXOffset_t unitSize = 1024000000 );
    // Do NOT rely on implicit deletion via the vector<unique_ptr<DataPage>>
    // member: std::vector<T>'s copy constructor is unconditionally declared
    // regardless of whether T is copyable, so std::is_copy_constructible
    // (and naive "try to copy it" code) will NOT reliably reflect the true
    // deletion status -- the actual failure only surfaces deep inside
    // vector's copy-ctor body when unique_ptr's deleted copy ctor is
    // odr-used, which type traits do not detect. Explicit deletion here
    // makes the intent unambiguous and correctly reported by type traits.
    DataBook( const DataBook & )             = delete;
    DataBook & operator=( const DataBook & ) = delete;

    // Declaring the copy ops above suppresses implicit move-op generation,
    // so restore them explicitly. This is cheap: unique_ptr elements are
    // moved via pointer transfer, O(1) regardless of how much data each
    // DataPage holds.
    DataBook( DataBook && )            = default;
    DataBook & operator=( DataBook && ) = default;

public:
    std::vector< std::unique_ptr<DataPage> > pages;   // renamed from dataBook
    HXOffset_t currPos;
    HXOffset_t maxUnitSize;
public:
    HXSize_t GetPageCount();   // moved from protected to public,
    // so tests can assert on page count directly
    DataPage * GetCurrentPage();
    DataPage * GetPage( HXSize_t iPage );
protected:
    void SetPageCount( HXSize_t newPageCount );
    HXOffset_t  GetRemainingSizeOfCurrentPage();
    void Advance( HXOffset_t offset );
public:
    void Read ( void * data, HXOffset_t dataSize );
    void Write( void * data, HXOffset_t dataSize );
    void ReadFile ( std::fstream & file );
    void WriteFile( std::fstream & file );

    void ReadString ( std::string & cs );
    void WriteString( std::string & cs );

    void Write( std::ostringstream * oss );

    HXOffset_t GetSize();
    void ReSize( HXOffset_t nLength );

    void Send( int pid, int tag );
    void Recv( int pid, int tag );

    void Bcast( int rootid );

    void SendRecv( int sendpid, int recvpid, int tag );

    void ToString( std::string & str );
    void Append( void * data, HXOffset_t dataSize );
    void AppendString( std::string & cs );

    void SecureRelativeSpace( HXOffset_t dataSize );
    void SecureAbsoluteSpace( HXOffset_t needSize );
    void MoveToBegin();
    void MoveToEnd();
};

void ToDataBook( DataBook * dataBook, std::ostringstream & oss );

EndNameSpace
