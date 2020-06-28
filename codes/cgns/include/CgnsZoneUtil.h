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
#include "HXArray.h"
#include "GridDef.h"
#include <vector>
using namespace std;

BeginNameSpace( ONEFLOW )

#ifdef ENABLE_CGNS

class Grid;
class StrGrid;
class CgnsZone;
class NodeMesh;

void EncodeIJK( int & index, int i, int j, int k, int ni, int nj, int nk );
void DecodeIJK( int index, int & i, int & j, int & k, int ni, int nj, int nk );
void GetRange( int ni, int nj, int nk, int startShift, int endShift, Range & I, Range & J, Range & K );
void GetIJKRegion( Range & I, Range & J, Range & K, int & ist, int & ied, int & jst, int & jed, int & kst, int & ked );

class PointSearch;
class BcRegion;
void PrepareCgnsZoneSub( Grids & grids, CgnsZone * cgnsZone );
void MergeToSingleZone( Grids & grids, HXVector< Int3D * > & unsIdList, NodeMesh * nodeMesh, int & nNode, int & nCell );
void FillSection( Grids & grids, HXVector< Int3D * > & unsIdList, CgnsZone * cgnsZone );
void CalcUnsId( StrGrid * grid, PointSearch * pointSearch, Int3D * unsId );
void SetUnsBcConn( BcRegion * bcRegion, CgIntField& conn, int & pos, Int3D & unsId );

void GenerateUnsBcElemConn( CgnsZone * myZone, CgnsZone * cgnsZoneIn );
void GenerateUnsBcCondConn( CgnsZone * myZone, CgnsZone * cgnsZoneIn );
void GenerateUnsVolElemConn( CgnsZone * myZone, CgnsZone * cgnsZoneIn );
void AllocateUnsElemConn( CgnsZone * myZone, CgnsZone * cgnsZoneIn );
void ReadElementConnectivities( CgnsZone * myZone, CgnsZone * cgnsZoneIn );
void GetStrZonePara( CgnsZone * myZone, int & s1, int & e1, int & s2, int & e2, int & etype1, int & etype2 );
void ReadCgnsZoneType( CgnsZone * myZone, CgnsZone * cgnsZoneIn );
void ReadCgnsZoneNameAndGeneralizedDimension( CgnsZone * myZone, CgnsZone * cgnsZoneIn );
void SetDimension( CgnsZone * myZone, CgnsZone * cgnsZoneIn );
void ReadCgnsGridCoordinates( CgnsZone * myZone, CgnsZone * cgnsZoneIn );
void ReadCgnsZoneAttribute( CgnsZone * myZone, CgnsZone * cgnsZoneIn );
void ReadCgnsGrid( CgnsZone * myZone, CgnsZone * cgnsZoneIn );
void DumpCgnsZoneType( CgnsZone * myZone, Grid * grid );
void FillISize( CgInt *isize, int ni, int nj, int nk, int dimension );
void FillISize( CgnsZone * myZone, Grid * gridIn );
void DumpCgnsZoneNameAndGeneralizedDimension( CgnsZone * myZone, Grid * gridIn );
void DumpCgnsZoneAttribute( CgnsZone * myZone, Grid * grid );
void DumpCgnsGridBoundary( CgnsZone * myZone, Grid * grid );
void DumpCgnsGridCoordinates( CgnsZone * myZone, Grid * grid );
void DumpCgnsZone( CgnsZone * myZone, Grid * grid );
void PrepareCgnsZone( CgnsZone * myZone, Grid * grid );

#endif

EndNameSpace