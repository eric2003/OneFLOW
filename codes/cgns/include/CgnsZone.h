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
#include <vector>
#include <string>
#include <fstream>
#include "HXCgns.h"
using namespace std;

BeginNameSpace( ONEFLOW )

#ifdef ENABLE_CGNS

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
class CgnsData;
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

	cgsize_t nNode, nCell;
    int nCoor;

	ZoneType_t cgnsZoneType;

	int zId;
	cgsize_t irmin[ 3 ], irmax[ 3 ], cellSize[ 3 ];
	cgsize_t isize[ 3 ][ 3 ];

	string coorName;
	string zoneName;

    IntField l2g;

    Real minLen, maxLen;
public:
	void FreeMesh();
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
    void GetElementNodeId( int eId, IntField & eNodeId );
    void ReadCgnsGrid();
    void ReadCgnsGrid( CgnsZone * cgnsZoneIn );
    void ReadCgnsZoneAttribute();
    void ReadCgnsZoneAttribute( CgnsZone * cgnsZoneIn );
    void ReadCgnsZoneType();
    void ReadCgnsZoneType( CgnsZone * cgnsZoneIn );
    void ReadCgnsZoneNameAndGeneralizedDimension();
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
    void ReadCgnsGridCoordinates( CgnsZone * cgnsZoneIn );
    void ReadCgnsGridBoundary();
	void ProcessPeriodicBc();
public:
    void AllocateUnsElemConn( CgnsZone * cgnsZoneIn );
    void GenerateUnsVolElemConn( CgnsZone * cgnsZoneIn );
    void GenerateUnsBcElemConn ( CgnsZone * cgnsZoneIn );
    void GenerateUnsBcCondConn ( CgnsZone * cgnsZoneIn );
	void SetElemPosition();
    void CreateCgnsBcRegion( CgnsZone * cgnsZoneIn );
    void InitL2g();
	cgsize_t GetNI() const;
	cgsize_t GetNJ() const;
	cgsize_t GetNK() const;
    void FillCgnsData( CgnsData * cgnsData );
};

void EncodeIJK( int & index, int i, int j, int k, int ni, int nj, int nk );
void DecodeIJK( int index, int & i, int & j, int & k, int ni, int nj, int nk );
void GetRange( int ni, int nj, int nk, int startShift, int endShift, Range & I, Range & J, Range & K );
void GetIJKRegion( Range & I, Range & J, Range & K, int & ist, int & ied, int & jst, int & jed, int & kst, int & ked );

#endif

EndNameSpace