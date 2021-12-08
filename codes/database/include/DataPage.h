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
#include "HXType.h"
#include <vector>
#include <fstream>
#include <string>

BeginNameSpace( ONEFLOW )

class DataPage
{
public:
    DataPage();
    ~DataPage();
public:
    typedef std::vector< char > CharMemory;
public:
    UInt GetSize();
    void Read ( void * data, UInt dataSize );
    void Read ( void * data, UInt dataSize, UInt position );
    void Write( void * data, UInt dataSize );
    void Write( void * data, UInt dataSize, UInt position );
    void ReadFile ( std::fstream & file );
    void WriteFile( std::fstream & file );
    void ToString( std::string & str );

    char * GetBeginDataPointer();
    char * GetCurrentDataPointer();

    void MoveToBegin() { MoveToPosition( 0 ); };
    void MoveToEnd  () { currPos = GetSize(); };
    void ReSize( UInt newSize );
    void Send( int pId, int tag );
    void Recv( int pId, int tag );
    void Bcast( int rootid );
protected:
    void MoveToPosition( UInt position );
    void MoveForwardPosition( UInt dataSize );
protected:
    UInt currPos;
    CharMemory * dataMemory;
public:
    char * GetDataPointer( int begin ) { return &( ( * dataMemory )[ begin ] ); }
};

EndNameSpace
