/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2024 He Xin and the OneFLOW contributors.
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
#include "HXCgns.h"
#include <string>
#include <map>

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
    int baseId;
    int nZones;
    int celldim, phydim;
    std::string baseName;
    HXVector< CgnsZone * > cgnsZones;
    bool freeFlag;
public:
    void FreeZoneList();
    CgnsZone * GetCgnsZoneByName( const std::string & zoneName );
    CgnsZone * GetCgnsZone( int iZone );
    void ConstructZoneNameMap();
    std::map< std::string, int > zoneNameMap;
    CgnsFamilyBc * familyBc;
public:
    int GetNZones();
    void SetDefaultCgnsBaseBasicInfo();
    void AddCgnsZone( CgnsZone * cgnsZone );
    void AllocateAllCgnsZones();
    void ReadCgnsBaseBasicInfo();
    void DumpCgnsBaseBasicInfo();
    void ReadNumberOfCgnsZones();
    CgnsZone * CreateCgnsZone();
    void CreateCgnsZones( int nZones );
    void ReadAllCgnsZones();
    void DumpAllCgnsZones();
    void ProcessCgnsZones();
    void ConvertToInnerDataStandard();
public:
    void SetFamilyBc( BCType_t & bcType, const std::string & bcRegionName );
    BCType_t GetFamilyBcType( const std::string & bcFamilyName );
    void ReadFamilySpecifiedBc();
public:
    void GoToBase();
    void GoToNode( const std::string & nodeName, int ith );
public:
    CgnsZone * WriteZoneInfo( const std::string & zoneName, ZoneType_t zoneType, cgsize_t * isize );
    CgnsZone * WriteZone( const std::string & zoneName );
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
