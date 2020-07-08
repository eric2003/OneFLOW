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
#include "HXCgns.h"
#include <string>
#include <map>
using namespace std;

BeginNameSpace( ONEFLOW )

#ifdef ENABLE_CGNS

class CgnsFile;
class CgnsZone;
class CgnsFamilyBc;

class CgnsBase
{
public:
    CgnsBase( CgnsFile * cgnsFile );
    CgnsBase();
    ~CgnsBase();
public:
    CgnsFile * cgnsFile;
    int fileId, baseId;
    int nZones;
    int celldim, phydim;
    string baseName;
    HXVector< CgnsZone * > cgnsZones;
    bool freeFlag;
public:
    void FreeZoneList();
    CgnsZone * GetCgnsZoneByName( const string & zoneName );
    CgnsZone * GetCgnsZone( int iZone );
    void ConstructZoneNameMap();
    map< string, int > zoneNameMap;
    CgnsFamilyBc * familyBc;
public:
    int GetNZone();
    void SetDefaultCgnsBaseBasicInfo();
    void AddCgnsZone( CgnsZone * cgnsZone );
    void AllocateAllCgnsZones();
    void ReadCgnsBaseBasicInfo();
    void DumpCgnsBaseBasicInfo();
    void ReadNumberOfCgnsZones();
    void ReadAllCgnsZones();
public:
    void SetFamilyBc( BCType_t & bcType, const string & bcRegionName );
    BCType_t GetFamilyBcType( const string & bcFamilyName );
    void ReadFamilySpecifiedBc();
public:
    void GoToBase();
    void GoToNode( const string & nodeName, int ith );
public:
    CgnsZone * WriteZoneInfo( const string & zoneName, ZoneType_t zoneType, cgsize_t * isize );
    CgnsZone * WriteZone( const string & zoneName );
    void SetTestISize( cgsize_t * isize );
    void ReadArray();
    void ReadReferenceState();
    void ReadBaseDescriptor();
    void ReadConvergence();
    void ReadCgnsZones();
    void ReadFlowEqn();
};

#endif

EndNameSpace