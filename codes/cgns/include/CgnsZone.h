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
#include "HXDefine.h"
#include "HXCgns.h"

BeginNameSpace( ONEFLOW )

#ifdef ENABLE_CGNS

class CgnsZone;
class CgnsBase;
class CgnsCoor;
class NodeMesh;
class CgnsZsection;
class CgnsZbc;
class GridElem;
class PointFactory;
class ElemFeature;
class FaceSolver;
class CgnsBase;

class CgnsZone
{
public:
    CgnsZone( CgnsBase * cgnsBase );
    ~CgnsZone();
public:
    CgnsBase * cgnsBase;
    CgnsCoor * cgnsCoor;
    CgnsZsection * cgnsZsection;
    CgnsZbc * cgnsZbc;
public:
    std::string zoneName;
    int zId;

    ZoneType_t cgnsZoneType;

    int volBcType;

    IntField l2g;

    CgInt isize[ 9 ];
public:
    void InitISize();
    void CopyISize( CgInt * isize );
    void SetVolBcType( int volBcType );
    int  GetVolBcType();
public:
    void Create();
    void SetPeriodicBc();
    void ConstructCgnsGridPoints( PointFactory * point_factory );
    void SetElementTypeAndNode( ElemFeature  * elem_feature );
    void InitLgMapping();
    void ConvertToInnerDataStandard();
public:
    void ScanBcFace( FaceSolver * face_solver );
    void GetElementNodeId( CgInt eId, CgIntField & eNodeId );
    void ReadCgnsGrid();
    void DumpCgnsGrid();
    void ReadCgnsZoneAttribute();
    void DumpCgnsZoneAttribute();
    void ReadCgnsZoneType();
    void DumpCgnsZoneType();
    void ReadCgnsZoneNameAndGeneralizedDimension();
    void DumpCgnsZoneNameAndGeneralizedDimension();
    void WriteZoneInfo( const std::string & zoneName, ZoneType_t zoneType, cgsize_t * isize );
    void SetDimension();
    void ReadElementConnectivities();
    void DumpElementConnectivities();
    void ReadNumberOfCgnsSections();
    void CreateCgnsSections();
    void ReadCgnsSections();
    void DumpCgnsSections();
    void ReadCgnsGridCoordinates();
    void DumpCgnsGridCoordinates();
    void ReadCgnsGridBoundary();
    void DumpCgnsGridBoundary();
    void ProcessPeriodicBc();
    void ReadCgnsZoneBasicInfo();
    void ReadCgnsGridCoordinates( CgnsZone * cgnsZoneIn );
public:
    void SetElemPosition();
public:
    CgInt GetNI() const;
    CgInt GetNJ() const;
    CgInt GetNK() const;
public:
    bool ExistSection( const std::string & sectionName );
    void GoToZone();
    void GoToNode( const std::string & nodeName, int ith );
    void GoToNode( const std::string & nodeNamei, int ith, const std::string & nodeNamej, int jth );
    void GoToNode( const std::string & nodeNamei, int ith, const std::string & nodeNamej, int jth, const std::string & nodeNamek, int kth );
public:
    void ReadFlowEqn();
};

#endif

EndNameSpace
