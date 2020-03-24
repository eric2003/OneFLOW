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
#include "HXDefine.h"
#include <string>
#include <vector>
using namespace std;

BeginNameSpace( ONEFLOW )

class CgnsBase;

string GetCgnsFileTypeName( int file_type );

class CgnsFile
{
public:
    CgnsFile();
    CgnsFile( const string & fileName, int openMode );
    ~CgnsFile();
public:
    string fileName;
    int fileId;
    int openMode;
    int openStatus;
    int nBases;
    vector< CgnsBase * > baseList;
public:
    int currBaseId;
public:
    void OpenCgnsFile( const string & fileName, int cgnsOpenMode );
    void CloseCgnsFile();
    CgnsBase * WriteBase( const string & baseName );
    CgnsBase * WriteBase( const string & baseName, int celldim, int physdim );
private:
    CgnsBase * AddBase( int fileId, const string & baseName, int celldim, int physdim, int baseId );
    void FreeBaseList();
public:
    void GoPath( const string & path );
    void ReadNumberOfBases();
    void ReadBases();
    void ReadArray();
    void ReadReferenceState();
    void ReadBaseDescriptor();
    void WriteBaseDescriptor();
    void ReadConvergence();
    void ReadFlowEqn();
};


EndNameSpace