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

#include "Boundary.h"
#include "HXCgns.h"
#include "Dimension.h"
#include "FileUtil.h"
#include "Prj.h"
#include "BcRecord.h"

#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

int BC::PERIODIC         =   99;
int BC::WAKE             = - 2 ;
int BC::INTERFACE        = - 1 ;
int BC::NO_BOUNDARY      =   0 ;
int BC::EXTRAPOLATION    =   1 ;
int BC::SOLID_SURFACE    =   2 ;
int BC::SYMMETRY         =   3 ;
int BC::FARFIELD         =   4 ;
int BC::INFLOW           =   5 ;
int BC::OUTFLOW          =   6 ;
int BC::OUTFLOW_CONFINED =   61;
int BC::POLE             =   7 ;
int BC::POLE1            =   71;
int BC::POLE2            =   72;
int BC::POLE3            =   73;
int BC::GENERIC_1        =   8 ;
int BC::GENERIC_2        =   9 ;
int BC::GENERIC_3        =   10;
int BC::INTERIOR         =   100;
int BC::OVERSET          =   1000;

BC::BC()
{
    ;
}

BC::~BC()
{
    ;
}

bool BC::IsInterfaceBc( int bcType )
{
    if ( bcType == BC::INTERFACE || bcType < 0 )
    {
        return true;
    }
    return false;
}

bool BC::IsSlipfaceBc( int bcType )
{
    if ( bcType == BC::GENERIC_2 )
    {
        return true;
    }
    return false;
}

bool BC::IsPoleBc( int bcType )
{
    int bc10 = ( bcType / 10 );
    if ( ( bcType == BC::POLE ) || ( bc10 == BC::POLE ) )
    {
        return true;
    }
    return false;
}

bool BC::IsNotNormalBc( int bcType )
{
    return BC::IsInterfaceBc( bcType ) || BC::IsPoleBc( bcType );
}

bool BC::IsWallBc( int bcType )
{
    if ( bcType == BC::SOLID_SURFACE )
    {
        return true;
    }
    return false;
}

BcTypeMap::BcTypeMap()
{
    ;
}

BcTypeMap::~BcTypeMap()
{
    ;
}

void BcTypeMap::Init()
{
    typedef pair< int, int > IntPair;

    cgns2OneFlow.insert( IntPair( CGNS_ENUMV( BCTypeUserDefined       ), BC::GENERIC_2       ) );
    cgns2OneFlow.insert( IntPair( CGNS_ENUMV( BCSymmetryPlane         ), BC::SYMMETRY        ) );
    cgns2OneFlow.insert( IntPair( CGNS_ENUMV( BCInflow                ), BC::INFLOW          ) );
    cgns2OneFlow.insert( IntPair( CGNS_ENUMV( BCInflowSubsonic        ), BC::INFLOW          ) );
    cgns2OneFlow.insert( IntPair( CGNS_ENUMV( BCInflowSupersonic      ), BC::INFLOW          ) );
    cgns2OneFlow.insert( IntPair( CGNS_ENUMV( BCTunnelInflow          ), BC::INFLOW          ) );
    cgns2OneFlow.insert( IntPair( CGNS_ENUMV( BCOutflow               ), BC::OUTFLOW         ) );
    cgns2OneFlow.insert( IntPair( CGNS_ENUMV( BCOutflowSubsonic       ), BC::OUTFLOW         ) );
    cgns2OneFlow.insert( IntPair( CGNS_ENUMV( BCOutflowSupersonic     ), BC::OUTFLOW         ) );
    cgns2OneFlow.insert( IntPair( CGNS_ENUMV( BCTunnelOutflow         ), BC::OUTFLOW         ) );
    cgns2OneFlow.insert( IntPair( CGNS_ENUMV( BCWall                  ), BC::SOLID_SURFACE   ) );
    cgns2OneFlow.insert( IntPair( CGNS_ENUMV( BCWallInviscid          ), BC::SOLID_SURFACE   ) );
    cgns2OneFlow.insert( IntPair( CGNS_ENUMV( BCWallViscous           ), BC::SOLID_SURFACE   ) );
    cgns2OneFlow.insert( IntPair( CGNS_ENUMV( BCWallViscousHeatFlux   ), BC::SOLID_SURFACE   ) );
    cgns2OneFlow.insert( IntPair( CGNS_ENUMV( BCWallViscousIsothermal ), BC::SOLID_SURFACE   ) );
    cgns2OneFlow.insert( IntPair( CGNS_ENUMV( BCDegenerateLine        ), BC::POLE            ) );
    cgns2OneFlow.insert( IntPair( CGNS_ENUMV( BCTypeNull              ), BC::INTERFACE       ) );
    cgns2OneFlow.insert( IntPair( CGNS_ENUMV( BCFarfield              ), BC::FARFIELD        ) );
    cgns2OneFlow.insert( IntPair( CGNS_ENUMV( BCGeneral               ), BC::OVERSET         ) );

    oneFlow2Cgns.insert( IntPair( BC::INTERFACE    , CGNS_ENUMV( BCTypeNull        ) ) );
    oneFlow2Cgns.insert( IntPair( BC::EXTRAPOLATION, CGNS_ENUMV( BCExtrapolate     ) ) );
    oneFlow2Cgns.insert( IntPair( BC::SYMMETRY     , CGNS_ENUMV( BCSymmetryPlane   ) ) );
    oneFlow2Cgns.insert( IntPair( BC::INFLOW       , CGNS_ENUMV( BCInflow          ) ) );
    oneFlow2Cgns.insert( IntPair( BC::OUTFLOW      , CGNS_ENUMV( BCOutflow         ) ) );
    oneFlow2Cgns.insert( IntPair( BC::SOLID_SURFACE, CGNS_ENUMV( BCWall            ) ) );
    oneFlow2Cgns.insert( IntPair( BC::POLE         , CGNS_ENUMV( BCDegenerateLine  ) ) );
    oneFlow2Cgns.insert( IntPair( BC::POLE1        , CGNS_ENUMV( BCDegenerateLine  ) ) );
    oneFlow2Cgns.insert( IntPair( BC::POLE2        , CGNS_ENUMV( BCDegenerateLine  ) ) );
    oneFlow2Cgns.insert( IntPair( BC::POLE3        , CGNS_ENUMV( BCDegenerateLine  ) ) );
    oneFlow2Cgns.insert( IntPair( BC::NO_BOUNDARY  , CGNS_ENUMV( BCTypeNull        ) ) );
    oneFlow2Cgns.insert( IntPair( BC::FARFIELD     , CGNS_ENUMV( BCFarfield        ) ) );
    oneFlow2Cgns.insert( IntPair( BC::OVERSET      , CGNS_ENUMV( BCGeneral         ) ) );
    oneFlow2Cgns.insert( IntPair( BC::GENERIC_1    , CGNS_ENUMV( BCGeneral         ) ) );
    oneFlow2Cgns.insert( IntPair( BC::GENERIC_2    , CGNS_ENUMV( BCTypeUserDefined ) ) );
    oneFlow2Cgns.insert( IntPair( BC::GENERIC_3    , CGNS_ENUMV( BCTypeUserDefined ) ) );
}

int BcTypeMap::OneFlow2Cgns( int oneflow_bctype )
{
    map< int, int >::iterator iter;
    iter = oneFlow2Cgns.find( oneflow_bctype );
    int cgns_bctype = BCTypeUserDefined;
    if ( iter != oneFlow2Cgns.end() )
    {
        cgns_bctype = iter->second;
    }
    return cgns_bctype;
}

int BcTypeMap::Cgns2OneFlow( int cgns_bctype )
{
    map< int, int >::iterator iter;
    iter = cgns2OneFlow.find( cgns_bctype );
    int oneflow_bctype = BC::GENERIC_2;
    if ( iter != cgns2OneFlow.end() )
    {
        oneflow_bctype = iter->second;
    }
    return oneflow_bctype;
}

CommonNameMap::CommonNameMap()
{
}

CommonNameMap::~CommonNameMap()
{
}

void CommonNameMap::AddName( const std::string & name )
{
    HXSort< std::string > data;
    data.value = name;
    data.index = 0;

    set< HXSort< std::string > >::iterator iter = this->stringMap.find( data );

    if ( iter == this->stringMap.end() )
    {
        data.index = this->stringMap.size();
        stringMap.insert( data );
    }
}

int CommonNameMap::FindNameId( const std::string & name )
{
    HXSort< std::string > data;
    data.value = name;
    data.index = 0;

    set< HXSort< std::string > >::iterator iter = this->stringMap.find( data );

    if ( iter == this->stringMap.end() )
    {
        return - 1;
    }
    else
    {
        return iter->index;
    }
}

void DumpRegion( const string & fileName, CommonNameMap & nameMap )
{
    fstream file;
    ONEFLOW::OpenPrjFile( file, fileName, ios_base::out );

    set< HXSort< std::string > > & stringMap = nameMap.GetNameMap();

    file << stringMap.size() << endl;

    for ( std::set< HXSort< std::string > >::iterator iter = stringMap.begin(); iter != stringMap.end(); ++ iter )
    {
        file << iter->index << " " << iter->value << endl;
    }
    CloseFile( file );
}

CommonNameMap RegionNameMap::nameMap;

RegionNameMap::RegionNameMap()
{
}

RegionNameMap::~RegionNameMap()
{
}

void RegionNameMap::AddRegion( const std::string & regionName )
{
    RegionNameMap::nameMap.AddName( regionName );
}

int RegionNameMap::FindRegionId( const std::string & regionName )
{
    return RegionNameMap::nameMap.FindNameId( regionName );
}

void RegionNameMap::DumpRegion()
{
    string fileName = "grid/bcRegionMap.txt";
    ONEFLOW::DumpRegion( fileName, RegionNameMap::nameMap );
}

CommonNameMap VolumeNameMap::nameMap;

VolumeNameMap::VolumeNameMap()
{
}

VolumeNameMap::~VolumeNameMap()
{
}

void VolumeNameMap::AddRegion( const std::string & regionName )
{
    VolumeNameMap::nameMap.AddName( regionName );
}

int VolumeNameMap::FindRegionId( const std::string & regionName )
{
    return VolumeNameMap::nameMap.FindNameId( regionName );
}

void VolumeNameMap::DumpRegion()
{
    string fileName = "grid/volumeRegionMap.txt";
    ONEFLOW::DumpRegion( fileName, VolumeNameMap::nameMap );
}


EndNameSpace