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

class CgnsMultiBase
{
public:
    CgnsMultiBase ();
    ~CgnsMultiBase();
public:
    int fileId, nBases;
    int nTZones;

    int volBcType;

    HXVector< CgnsBase * > baseVector;
    IntField zid1, zid2;
public:
    int GetSystemZoneType();
    void ReadCgnsGrid();
    void ReadCgnsGrid( const string & fileName );
    void OpenCgnsFile( const string & fileName, int cgnsOpenMode );
    void CloseCgnsFile();
    void ReadCgnsMultiBase();
    void ReadNumCgnsBase();

    void ReadCgnsMultiBase( CgnsMultiBase * strCgnsMultiBase );
    void ReadNumCgnsBase( CgnsMultiBase * strCgnsMultiBase );
    void ReadGeneralizedCgnsZoneScale();
    void ReadGeneralizedCgnsZoneScale( CgnsMultiBase * strCgnsMultiBase );
    void ReadAllCgnsZonesInEachCgnsBase();
    void ReadAllCgnsZonesInEachCgnsBase( CgnsMultiBase * strCgnsMultiBase );
public:
    void Create( int nZones );
    void InitDefaultCgnsBase();
    void AllocateCgnsBase();
    void InitCgnsBase();
    void SetGeneralizedCgnsZoneScale( int nZones );
    void ComputeNumberOfTotalZones();
    void AllocateCgnsZonesInEachCgnsBase();
    void InitAllCgnsZonesInEachCgnsBase();
    void ConvertStrCgns2UnsCgnsGrid( CgnsMultiBase * strCgnsMultiBase );
public:
    CgnsZone * GetZone( int iZone );
    int FindBaseId( int iZone );
};

#endif

EndNameSpace