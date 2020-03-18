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
#include "HXDefine.h"
#include <fstream>
using namespace std;

BeginNameSpace( ONEFLOW )

class DataBook;

class VirtualFile
{
public:
    VirtualFile( fstream * file );
    VirtualFile( DataBook * databook );
public:
    int type;
    fstream * file;
    DataBook * databook;
public:
    streamsize sectionBegin, sectionEnd, dataLength;
public:
    void Read ( void * data, streamsize size );
    void Write( void * data, streamsize size );
protected:
    void MarkSectionBegin();
    void MarkSectionEnd  ();
    streamsize GetSectionBegin() { return sectionBegin; };
    streamsize GetSectionEnd() { return sectionEnd; };
    void MoveToPosition( streamsize sectionPosition );
    void MoveToSectionBegin();
    void MoveToSectionEnd();
    void WriteDataLength();
public:
    void BeginWriteWork();
    void EndWriteWork();
    void BeginReadWork();
    void EndReadWork();

    void ReadDataLength();
    void ReservePlaceholder();
    streamsize GetCurrentPosition();
};

EndNameSpace
