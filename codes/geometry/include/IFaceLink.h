/*---------------------------------------------------------------------------*\\
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

\\*---------------------------------------------------------------------------*/

#pragma once
#include "HXDefine.h"
#include "HXLookup.h"
#include "GridHandles.h"
#include <memory>
#include <set>

BeginNameSpace( ONEFLOW )

class Grid;
class FaceSearch;
class PointLocator;
class NodeMesh;

class IFaceLink
{
public:
    explicit IFaceLink( Grids & gridsIn );
    ~IFaceLink();
public:
    HXLookup<int> faceLookup;
    LinkField gI2Zid;
    LinkField g2l;
    LinkField l2g;

    LinkField gI2ZidNew;
    LinkField g2lNew;
    LinkField l2gNew;

    LinkField nChild;

    std::unique_ptr< FaceSearch > face_search;
    std::unique_ptr< PointLocator > point_search;

private:
    Grids & grids;

public:
    [[nodiscard]] Grid & GetGrid( int zoneIndex );
    void Init( Grid & grid );
public:
    void CreateLink( IntField & faceNode, int zid, int lCount );
    void MatchInterfaceTopology( Grid & grid );
    void MatchPeriodicInterface( Grid & grid );
    void ReconstructInterFace();
protected:
    void AddFace( const IntField & facePointIndexes );
public:
    void UpdateLgMapping();
    void InitNewLgMapping();
};

void GetFaceCoorList( const IntField & faceNode, RealField & xList, RealField & yList, RealField & zList, const NodeMesh & nodeMesh );
void GetCoorIdList( IFaceLink & iFaceLink, RealField & xList, RealField & yList, RealField & zList, int nPoint, IntField & pointId );

EndNameSpace
