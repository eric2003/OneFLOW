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
#include "HXSort.h"
#include "GridDef.h"
#include <set>
using namespace std;

BeginNameSpace( ONEFLOW )

class Grid;
class FaceSearch;
class PointSearch;
class NodeMesh;

class IFaceLink
{
public:
    IFaceLink( Grids & grids );
    ~IFaceLink();
public:
    //faceListis used primarily as a search list for global face
    set < HXSort< IntField > > inFaceList;
    //In general, each interface is made up of two different blocks of surface.
    //This requires each surface to have a block number and the serial number of the surface in this block
    LinkField gI2Zid;
    LinkField g2l;
    LinkField l2g;

    LinkField gI2ZidNew;
    LinkField g2lNew;
    LinkField l2gNew;

    LinkField nChild;

    FaceSearch * face_search;

    PointSearch * point_search;

    Grids grids;
public:
    void Init( Grid * grid );
    Grid * GetGrid( int zoneIndex ) { return grids[ zoneIndex ]; }
public:
    void CreateLink( IntField & faceNode, int zid, int lCount );
    void MatchInterfaceTopology( Grid * grid );
    void MatchPeoridicInterface( Grid * grid );
    void ReconstructInterFace();
protected:
    void AddFace( const IntField & facePointIndexes );
public:
    void UpdateLgMapping();
    void InitNewLgMapping();
};

void GetFaceCoorList( IntField & faceNode, RealField & xList, RealField & yList, RealField & zList, NodeMesh * nodeMesh );
void GetCoorIdList( IFaceLink * iFaceLink, RealField & xList, RealField & yList, RealField & zList, int nPoint, IntField & pointId );

EndNameSpace