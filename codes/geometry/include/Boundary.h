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
#include "HXSort.h"
#include <string>
#include <map>
#include <set>
using namespace std;
BeginNameSpace( ONEFLOW )

class BC
{
public:
    BC();
    ~BC();
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
    BcTypeMap();
    ~BcTypeMap();
private:
    int numberOfMaxBoundaryConditions;
    map< int, int > cgns2OneFlow;
    map< int, int > oneFlow2Cgns;
public:
    void Init();
    int OneFlow2Cgns( int oneflow_bctype );
    int Cgns2OneFlow( int cgns_bctype );
};

class CommonNameMap
{
public:
    CommonNameMap();
    ~CommonNameMap();
protected:
    std::set< HXSort< std::string > > stringMap;
public:
    void AddName( const std::string & name );
    int  FindNameId( const std::string & name );
    set< HXSort< std::string > > & GetNameMap() { return stringMap; }
};

void DumpRegion( const string & fileName, CommonNameMap & nameMap );

class RegionNameMap
{
public:
    RegionNameMap();
    ~RegionNameMap();
public:
    static CommonNameMap nameMap;
public:
    static void AddRegion( const std::string & regionName );
    static int  FindRegionId( const std::string & regionName );
    static void DumpRegion();
};

class VolumeNameMap
{
public:
    VolumeNameMap();
    ~VolumeNameMap();
public:
    static CommonNameMap nameMap;
public:
    static void AddRegion( const std::string & regionName );
    static int  FindRegionId( const std::string & regionName );
    static void DumpRegion();
};

EndNameSpace