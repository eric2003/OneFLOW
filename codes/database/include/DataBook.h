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
#include "DataPage.h"
#include "HXType.h"
#include <vector>
#include <string>
#include <fstream>
using namespace std;

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
    vector< DataPage * > * dataBook;
    UInt currPageId;
    LLong currPos;
    LLong maxUnitSize;
public:
    DataPage * GetCurrentPage();
    DataPage * GetPage( UInt iPage );
    void Destroy( DataPage * dataPage );
    void Erase( UInt startPage, UInt endPage );
protected:
    UInt  GetNPage();
    void ResizeNPage( UInt newNPage );
    LLong  GetRemainingSizeOfCurrentPage();
    void MoveForwardPosition( LLong dataSize );
public:
    void Read ( void * data, LLong dataSize );
    void Write( void * data, LLong dataSize );
    void ReadFile ( fstream & file );
    void WriteFile( fstream & file );

    void ReadString ( string & cs );
    void WriteString( string & cs );

    void Write( ostringstream * oss );

    LLong GetSize();
    void ReSize( LLong nLength );

    void Send( int pid, int tag );
    void Recv( int pid, int tag );
    void Bcast( int rootid );

    void SendRecv( int sendpid, int recvpid, int tag );

    void ToString( string & str );
    void Append( void * data, LLong dataSize );
    void AppendString( string & cs );

    void SecureRelativeSpace( LLong dataSize );
    void SecureAbsoluteSpace( LLong needSize );
    void MoveToBegin();
    void MoveToEnd();
};

void ToDataBook( DataBook * dataBook, ostringstream & oss );

EndNameSpace
