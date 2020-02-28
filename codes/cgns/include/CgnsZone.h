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
#include <vector>
#include <string>
#include <fstream>
#include "HXCgns.h"
#include "HXArray.h"
#include "GridDef.h"
using namespace std;

BeginNameSpace( ONEFLOW )

#ifdef ENABLE_CGNS

class Grid;
class StrGrid;
class CgnsZone;
class CgnsBase;
class NodeMesh;
class CgnsMultiSection;
class CgnsBcRegionProxy;
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
    NodeMesh * nodeMesh;
    CgnsMultiSection * multiSection;
    CgnsBcRegionProxy * bcRegionProxy;

    CgInt nNode, nCell;
    int nCoor;

    ZoneType_t cgnsZoneType;
    int volBcType;

    int zId;
    CgInt irmin[ 3 ], irmax[ 3 ], cellSize[ 3 ];
    CgInt isize[ 9 ];

    string zoneName;

    IntField l2g;

    Real minLen, maxLen;
public:
    void FreeMesh();
    void SetVolBcType( int volBcType );
    int GetVolBcType();
public:
    void Create();
    void SetPeriodicBc();
    void InitElement( GridElem * ge );
    void ConstructCgnsGridPoints( PointFactory * point_factory );
    void SetElementTypeAndNode( ElemFeature  * elem_feature );
    void InitLgMapping();
    void ConvertToInnerDataStandard();
public:
    void ScanBcFace( FaceSolver * face_solver );
    void GetElementNodeId( CgInt eId, CgIntField & eNodeId );
    void ReadCgnsGrid();
    void ReadCgnsGrid( CgnsZone * cgnsZoneIn );
    void ReadCgnsZoneAttribute();
    void DumpCgnsZoneAttribute( Grid * grid );
    void ReadCgnsZoneAttribute( CgnsZone * cgnsZoneIn );
    void ReadCgnsZoneType();
    void DumpCgnsZoneType( Grid * grid );
    void ReadCgnsZoneType( CgnsZone * cgnsZoneIn );
    void ReadCgnsZoneNameAndGeneralizedDimension();
    void DumpCgnsZoneNameAndGeneralizedDimension( Grid * gridIn );
    void ReadCgnsZoneNameAndGeneralizedDimension( CgnsZone * cgnsZoneIn );
    void SetDimension();
    void SetDimension( CgnsZone * cgnsZoneIn );
    void ReadElementConnectivities();
    void ReadElementConnectivities( CgnsZone * cgnsZoneIn );
    void ReadNumberOfCgnsSections();
    void ReadNumberOfCgnsSections( CgnsZone * cgnsZoneIn );
    void CreateCgnsSections();
    void ReadCgnsSections();
    void ReadCgnsGridCoordinates();
    void DumpCgnsGridCoordinates( Grid * grid );
    void ReadCgnsGridCoordinates( CgnsZone * cgnsZoneIn );
    void ReadCgnsGridBoundary();
    void DumpCgnsGridBoundary( Grid * grid );
    void ProcessPeriodicBc();
    void DumpCgnsZone( Grid * grid );
    void FillISize( Grid * gridIn );
    void FillISize( int ni, int nj, int nk, int dimension );
    void PrepareCgnsZone( Grid * grid );
public:
    void AllocateUnsElemConn( CgnsZone * cgnsZoneIn );
    void GenerateUnsVolElemConn( CgnsZone * cgnsZoneIn );
    void GenerateUnsBcElemConn ( CgnsZone * cgnsZoneIn );
    void GenerateUnsBcCondConn ( CgnsZone * cgnsZoneIn );
    void SetElemPosition();
    void CreateCgnsBcRegion( CgnsZone * cgnsZoneIn );
    void InitL2g();
    CgInt GetNI() const;
    CgInt GetNJ() const;
    CgInt GetNK() const;

    void GetStrZonePara( int & s1, int & e1, int & s2, int & e2, int & etype1, int & etype2 );
public:
    bool ExistSection( const string & sectionName );
};

void EncodeIJK( int & index, int i, int j, int k, int ni, int nj, int nk );
void DecodeIJK( int index, int & i, int & j, int & k, int ni, int nj, int nk );
void GetRange( int ni, int nj, int nk, int startShift, int endShift, Range & I, Range & J, Range & K );
void GetIJKRegion( Range & I, Range & J, Range & K, int & ist, int & ied, int & jst, int & jed, int & kst, int & ked );

class PointSearch;
class BcRegion;
void PrepareCgnsZone( Grids & grids, CgnsZone * cgnsZone );
void MergeToSingleZone( Grids & grids, HXVector< Int3D * > & unsIdList, NodeMesh * nodeMesh, int & nNode, int & nCell );
void FillSection( Grids & grids, HXVector< Int3D * > & unsIdList, CgnsZone * cgnsZone );
void ComputeUnsId( StrGrid * grid, PointSearch * pointSearch, Int3D * unsId );
void SetUnsBcConn( BcRegion * bcRegion, CgIntField& conn, int & pos, Int3D & unsId );

#endif

EndNameSpace