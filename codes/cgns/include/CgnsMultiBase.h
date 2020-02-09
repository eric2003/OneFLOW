/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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

BeginNameSpace( ONEFLOW )

#ifdef ENABLE_CGNS

class CgnsBase;
class CgnsZone;
class GridMediator;

class CgnsMultiBase
{
public:
    CgnsMultiBase ();
    ~CgnsMultiBase();
public:
    int fileId, nBases;
    int nTZones;

    HXVector< CgnsBase * > baseVector;
public:
    int GetSystemZoneType();
    void ReadCgnsGrid();
    void DumpCgnsGrid( GridMediator * gridMediator );
    void ReadCgnsGrid( const string & fileName );
    void OpenCgnsFile( const string & fileName, int cgnsOpenMode );
    void CloseCgnsFile();
    void ReadCgnsMultiBase();
    void DumpCgnsMultiBase( GridMediator * gridMediator );
    void ReadNumCgnsBase();

    void ReadCgnsMultiBase( CgnsMultiBase * strCgnsMultiBase );
    void ReadNumCgnsBase( CgnsMultiBase * strCgnsMultiBase );
public:
    void CreateDefaultCgnsZones( int nZones );
    void InitDefaultCgnsBase();
    void InitCgnsBase();
    void ComputeNumberOfTotalZones();
    void ConvertStrCgns2UnsCgnsGrid( CgnsMultiBase * strCgnsMultiBase );
public:
    CgnsBase * GetCgnsBase( int iBase );
    CgnsZone * GetCgnsZone( int globalZoneId );
    CgnsZone * FindGlobalCgnsZone( int globalZoneId );
};

#endif

EndNameSpace