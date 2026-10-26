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

    // Interface-face relation, stored in two directions:
    // l2g[zone][localFace] -> global interface-face ID.
    // gI2Zid means global Interface-face to zone IDs:
    // gI2Zid[globalFace][slot] -> participating zone ID.
    // g2l[globalFace][slot] -> that zone's local interface-face ID.
    // The slot index is shared by gI2Zid and g2l; their row lengths must match.
    LinkField gI2Zid;
    LinkField g2l;
    LinkField l2g;

    // Candidate tables used while interface faces are split and remapped.
    // They replace the live tables only after CommitInterfaceMappingBuffers
    // validates the paired inverse-map rows.
    LinkField gI2ZidNew;
    LinkField g2lNew;
    LinkField l2gNew;

    // Legacy per-zone child-count bookkeeping for interface-face rebuilding.
    LinkField nChild;

    std::unique_ptr< FaceSearch > face_search;
    std::unique_ptr< PointLocator > point_search;

private:
    Grids & grids;
    // PointLocator ADT build needs valid nodeMesh bounding boxes. Defer until a
    // code path actually searches/adds points so topology-only clients (and
    // unit tests) can construct IFaceLink without a full mesh.
    bool pointLocatorReady = false;
    void ValidateGridIndex( const Grid & grid, const char * operation ) const;

public:
    // Builds the PointLocator ADT on first use (requires valid nodeMesh boxes).
    void EnsurePointLocatorReady();
    [[nodiscard]] bool IsPointLocatorReady() const { return pointLocatorReady; }
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

// Validate one global interface face: id in range, and zoneIds[id] aligns
// with localFaceIds[id] (same length). operation is included in error text.
void ValidateGlobalFaceMapping(
    int globalFaceId,
    const LinkField & zoneIds,
    const LinkField & localFaceIds,
    const char * operation );

// Validate table sizes match, then each entry via ValidateGlobalFaceMapping.
void ValidateGlobalFaceMappings(
    const LinkField & zoneIds,
    const LinkField & localFaceIds,
    const char * operation );

// Register local interface face (zid, localFaceId) on global face gIid.
// isNew: append a new parallel row to gI2Zid/g2l; otherwise push onto the
// existing row. Always writes l2g[zid][localFaceId] = gIid.
// Throws on inconsistent table sizes or out-of-range indices.
void RegisterLocalInterfaceFace(
    int gIid,
    bool isNew,
    int zid,
    int localFaceId,
    LinkField & l2g,
    LinkField & gI2Zid,
    LinkField & g2l );

// Among zone/local pairs on globalFaceId, find the first entry that is not
// (selfZoneId, selfLocalFaceId). Returns false if none.
bool FindRemoteInterfacePartner(
    int globalFaceId,
    int selfZoneId,
    int selfLocalFaceId,
    const LinkField & zoneIdsPerGlobalFace,
    const LinkField & localFaceIdsPerGlobalFace,
    int & outRemoteZoneId,
    int & outRemoteLocalFaceId );

// Count local interface faces whose global face has fewer than two zone
// references (unpaired). Faces with more than two zones are not counted here
// (legacy MatchInterfaceTopology behavior).
int CountInterfaceFacesMissingPartner(
    const IntField & localToGlobal,
    const LinkField & zoneIdsPerGlobalFace,
    const LinkField & localFaceIdsPerGlobalFace );

// Copy current inverse maps into *New buffers and clear/resize l2gNew to
// match the number of zones in l2g (empty rows ready for rebuild).
void InitNewInterfaceMappingBuffers(
    const LinkField & gI2Zid,
    const LinkField & g2l,
    const LinkField & l2g,
    LinkField & gI2ZidNew,
    LinkField & g2lNew,
    LinkField & l2gNew );

// Validate *New inverse maps then swap them into the live tables (including l2g).
void CommitInterfaceMappingBuffers(
    LinkField & gI2Zid,
    LinkField & g2l,
    LinkField & l2g,
    LinkField & gI2ZidNew,
    LinkField & g2lNew,
    LinkField & l2gNew,
    const char * operation );

EndNameSpace

