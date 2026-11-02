/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
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
#include "NamespaceMacros.h"
#include "HXDefine.h"
#include <string>
#include <map>
#include <vector>
#include <set>

BeginNameSpace( ONEFLOW )

class BC
{
public:
    BC() = default;
    ~BC() = default;
public:
    static int PERIODIC        ;
    static int WAKE            ;
    static int INTERFACE       ;
    static int NO_BOUNDARY     ;
    static int EXTRAPOLATION   ;
    static int SOLID_SURFACE   ;
    static int SYMMETRY        ;
    static int FARFIELD        ;
    static int INFLOW          ;
    static int OUTFLOW         ;
    static int OUTFLOW_CONFINED;
    static int POLE            ;
    static int POLE1           ;
    static int POLE2           ;
    static int POLE3           ;
    static int GENERIC_1       ;
    static int GENERIC_2       ;
    static int GENERIC_3       ;
    static int INTERIOR        ;
    static int OVERSET         ;
public:
    static bool IsInterfaceBc( int bcType );
    static bool IsSlipfaceBc( int bcType );
    static bool IsPoleBc( int bcType );
    static bool IsNotNormalBc( int bcType );
    static bool IsWallBc( int bcType );
};

class BcTypeMap
{
public:
    BcTypeMap() = default;
    ~BcTypeMap() = default;
private:
    int numberOfMaxBoundaryConditions;
    std::map< int, int > cgns2OneFlow;
    std::map< int, int > oneFlow2Cgns;
public:
    void Init();
    int OneFlow2Cgns( int oneflow_bctype ) const;
    int Cgns2OneFlow( int cgns_bctype ) const;

    // Convert each entry from CGNS BC type ids to OneFLOW BC ids in place.
    void MapFromCgnsInPlace( IntField & bcTypes ) const;
};

class CommonNameMap
{
public:
    CommonNameMap() = default;
    ~CommonNameMap() = default;

    void AddName(const std::string& name);
    void Clear() noexcept { stringMap.clear(); }
    int FindNameId(const std::string& name) const;
    const std::map<std::string, int>& GetNameMap() const { return stringMap; }

private:
    std::map<std::string, int> stringMap;  // String -> ID mapping
};

void DumpRegion( const std::string & fileName, CommonNameMap & nameMap );

// Instance-owned region-name table. Prefer this in new code and tests so
// callers are not coupled to process-wide state.
class RegionNameRegistry
{
public:
    RegionNameRegistry() = default;
    ~RegionNameRegistry() = default;

    void Reset();
    void AddRegion( const std::string & regionName );
    int FindRegionId( const std::string & regionName ) const;
    void DumpRegion( const std::string & fileName );
    const CommonNameMap & GetNameMap() const { return nameMap; }

private:
    CommonNameMap nameMap;
};

// Process-wide default registry. Static methods forward to Default().
// Prefer injecting RegionNameRegistry& into APIs that need isolation.
class RegionNameMap
{
public:
    RegionNameMap() = default;
    ~RegionNameMap() = default;

    static RegionNameRegistry & Default();

    static void Reset();
    static void AddRegion( const std::string & regionName );
    static int  FindRegionId( const std::string & regionName );
    static void DumpRegion();
    static void DumpRegion( const std::string & fileName );
};

// Process-wide default for volume-region names. Same registry type as BC
// regions; kept as a separate Default() instance so the two namespaces
// cannot collide. Prefer injecting RegionNameRegistry& when isolation matters.
class VolumeNameMap
{
public:
    VolumeNameMap() = default;
    ~VolumeNameMap() = default;

    static RegionNameRegistry & Default();

    static void Reset();
    static void AddRegion( const std::string & regionName );
    static int  FindRegionId( const std::string & regionName );
    static void DumpRegion();
    static void DumpRegion( const std::string & fileName );
};


// Register each name in order; returns parallel name ids (same length as names).
// Mirrors the CgnsZbcBoco region-registration loop without CGNS types.
IntField RegisterRegionNames( RegionNameRegistry & registry, const std::vector< std::string > & names );

EndNameSpace
