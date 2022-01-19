/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2022 He Xin and the OneFLOW contributors.
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

BeginNameSpace( ONEFLOW )

class DataBook;

typedef void ( * DATA_COMPRESS )( DataBook *& dataBook );
typedef void ( * DATA_DECOMPRESS )( DataBook *  dataBook );

class DataBook
{
public:
    DataBook();
    ~DataBook();
public:
    std::vector< DataPage * > * dataBook;
    HXSize_t currPageId;
    HXLongLong_t currPos;
    HXLongLong_t maxUnitSize;
public:
    DataPage * GetCurrentPage();
    DataPage * GetPage( HXSize_t iPage );
    void Destroy( DataPage * dataPage );
    void Erase( HXSize_t startPage, HXSize_t endPage );
protected:
    HXSize_t  GetNPage();
    void ResizeNPage( HXSize_t newNPage );
    HXLongLong_t  GetRemainingSizeOfCurrentPage();
    void MoveForwardPosition( HXLongLong_t dataSize );
public:
    void Read ( void * data, HXLongLong_t dataSize );
    void Write( void * data, HXLongLong_t dataSize );
    void ReadFile ( std::fstream & file );
    void WriteFile( std::fstream & file );

    void ReadString ( std::string & cs );
    void WriteString( std::string & cs );

    void Write( std::ostringstream * oss );

    HXLongLong_t GetSize();
    void ReSize( HXLongLong_t nLength );

    void Send( int pid, int tag );
    void Recv( int pid, int tag );

    void Bcast( int rootid );

    void SendRecv( int sendpid, int recvpid, int tag );

    void ToString( std::string & str );
    void Append( void * data, HXLongLong_t dataSize );
    void AppendString( std::string & cs );

    void SecureRelativeSpace( HXLongLong_t dataSize );
    void SecureAbsoluteSpace( HXLongLong_t needSize );
    void MoveToBegin();
    void MoveToEnd();
};

void ToDataBook( DataBook * dataBook, std::ostringstream & oss );

EndNameSpace
