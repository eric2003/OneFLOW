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
#include <string_view>
#include <fstream>
#include <memory>
#include <sstream>


BeginNameSpace( ONEFLOW )

class DataBook;

using DATA_COMPRESS = void( * )( DataBook *& dataBook );
using DATA_DECOMPRESS = void( * )( DataBook *  dataBook );

class DataBook
{
public:
    explicit DataBook( HXOffset_t unitSize = 1024000000 );

    // Explicitly delete copy semantics to match type traits expectations.
    DataBook( const DataBook & )             = delete;
    DataBook & operator=( const DataBook & ) = delete;

    // Default move semantics for O(1) performance.
    DataBook( DataBook && )                  = default;
    DataBook & operator=( DataBook && )       = default;
    ~DataBook()                              = default;
public:
    // Read-only state queries (marked const)
    HXSize_t GetPageCount() const;
    HXOffset_t size() const;

    DataPage * GetCurrentPage();
    const DataPage * GetCurrentPage() const;

    DataPage * GetPage( HXSize_t iPage );
    const DataPage * GetPage( HXSize_t iPage ) const;
public:
    // Data Read/Write interfaces with const-safety
    void Read ( void * data, HXOffset_t dataSize );
    void Write( const void * data, HXOffset_t dataSize );

    void ReadString ( std::string & str );
    void WriteString( std::string_view str );
    void AppendString( std::string_view str );

    void Write( const std::ostringstream * oss );
    void Append( const void * data, HXOffset_t dataSize );

    // Memory management and seek operations
    void Resize( HXOffset_t nLength );
    void Reserve( HXOffset_t needSize );
    void MoveToBegin();
    void MoveToEnd();

    // I/O Operations
    void ReadFile ( std::fstream & file );
    void WriteFile( std::fstream & file ) const;
    void ToString ( std::string & str ) const;

    // MPI / Parallel Communication interfaces
    void Send( int pid, int tag ) const;
    void Recv( int pid, int tag );
    void Bcast( int rootid );
    void SendRecv( int sendpid, int recvpid, int tag );

protected:
    void SetPageCount( HXSize_t newPageCount );
    HXOffset_t  GetRemainingSizeOfCurrentPage() const;
    void Advance( HXOffset_t offset );

private:
    std::vector< std::unique_ptr<DataPage> > pages;
    HXOffset_t currPos{0};
    HXOffset_t maxUnitSize{1024000000};

    void SecureRelativeSpace( HXOffset_t dataSize );

};

void ToDataBook( DataBook * dataBook, std::ostringstream & oss );

EndNameSpace
