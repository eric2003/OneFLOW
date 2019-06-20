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
#include <map>
#include "HXCgns.h"
using namespace std;

BeginNameSpace( ONEFLOW )

#ifdef ENABLE_CGNS

class CgnsZone;
class CgnsBase;
class NodeMesh;
class FaceSolver;

class CgnsBcInterface;

class CgnsBcRegion
{
public:
	CgnsBcRegion( CgnsZone * cgnsZone );
	~CgnsBcRegion();
public:
    int id;
    int nameId;
	string name;
    double bc_double_id;

	BCType_t bcType;
	PointSetType_t pointSetType;
	GridLocation_t gridLocation;
    GridConnectivityType_t gridConnType;  //Overset, Abutting, Abutting1to1
	DataType_t normalDataType;
	cgsize_t normalListSize;
	int normalIndex[ 3 ];
	int nDataSets;

	cgsize_t nElements;

	vector<cgsize_t> connList;

    CgnsBcInterface * bcInterface;

    CgnsZone * cgnsZone;
public:
    void ConvertToInnerDataStandard();
    int  ComputeBase();
    void ShiftBcRegion();
    void ScanBcFace( FaceSolver * face_solver );
public:
    void ProcessVertexBc( IntSet & bcVertex );
    void ProcessFaceBc( IntSet & bcVertex );
public:
    void ReadCgnsOrdinaryBcRegion();
    void ReadCgnsOrdinaryBcRegionInfo();
    void ReadCgnsOrdinaryBcRegionGridLocation();
    void CreateCgnsBcConn();
    void ReadCgnsBcConn();
    void PrintCgnsBcConn();
    void ExtractIJKRegionFromBcConn( IntField & ijkMin, IntField & ijkMax, vector<cgsize_t>& bcConn );
    void ExtractIJKRegionFromBcConn( IntField & ijkMin, IntField & ijkMax );
public:
    void ProcessCgns1to1BcRegion( int bcId );
    void ReadCgns1to1BoundaryRegion( int iCgns1to1BoundaryRegion );
    void CopyStrBcRegion( CgnsBcRegion * strBcRegion, cgsize_t& startId );
    void ReadCgnsBcConn( CgnsBcRegion * strBcRegion, cgsize_t & startId );
	void ReconstructStrRegion( IntField & ijkMin, IntField & ijkMax );
    cgsize_t GetActualNumberOfBoundaryElements();
};

void SetBcConn( CgnsZone * cgnsZone, IntField & ijkMin, IntField & ijkMax, vector<cgsize_t>& conn, int & pos, int & nElem );

#endif

EndNameSpace