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
#include "HXCgns.h"
#include <string>
#include <map>
using namespace std;

BeginNameSpace( ONEFLOW )

#ifdef ENABLE_CGNS

class CgnsZone;
class CgnsFamilyBc;

class CgnsBase
{
public:
    CgnsBase();
    ~CgnsBase();
public:
    int fileId, baseId;
    int nZones;
    int celldim, phydim;
    string baseName;
    HXVector< CgnsZone * > cgnsZones;
public:
    CgnsZone * GetCgnsZone( int zoneId );
    CgnsZone * GetCgnsZone( const string & zoneName );
    void ConstructZoneNameMap();
    map< string, int > zoneNameMap;
    CgnsFamilyBc * familyBc;
public:
    void SetDefaultCgnsBaseBasicInformation();
    void AllocateAllCgnsZonesInCurrentCgnsBase();
    void InitAllCgnsZonesInCurrentCgnsBase();
    void ReadCgnsBaseBasicInfo();
    void ReadCgnsBaseBasicInfo( CgnsBase * cgnsBaseIn );
    void ReadNumberOfCgnsZones();
    void ReadNumberOfCgnsZones( CgnsBase * cgnsBaseIn );
    void ReadAllCgnsZones();
    void ReadAllCgnsZones( CgnsBase * cgnsBaseIn );
public:
    void SetFamilyBc( BCType_t & bcType, const string & bcRegionName );
    void ReadFamilySpecifiedBc();
};

#endif

EndNameSpace