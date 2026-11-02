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
struct InterfaceFaceSplitTable;

class IFaceLink
{
public:
    explicit IFaceLink( Grids & gridsIn );
    ~IFaceLink();
private:
    // Interface-face relation, stored in two directions:
    // l2g[zone][localFace] -> global interface-face ID.
    // gI2Zid means global Interface-face to zone IDs:
    // gI2Zid[globalFace][slot] -> participating zone ID.
    // g2l[globalFace][slot] -> that zone's local interface-face ID.
    // The slot index is shared by gI2Zid and g2l; their row lengths must match.
    // Access only through Ensure*/Register*/Get*/Mutable* APIs.
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

    Grids & grids;
    // Interface face geometry search + split tables. Prefer GetFaceSearch /
    // MutableInterfaceSplit over direct field access.
    std::unique_ptr< FaceSearch > face_search;
    // Geometry identity for interface faces (CreateLink only). Not part of the
    // dual-map tables; kept private so registration cannot bypass CreateLink.
    HXLookup<int> faceLookup;
    // Built lazily via EnsurePointLocatorReady (needs nodeMesh bounds).
    std::unique_ptr< PointLocator > point_search;
    // PointLocator ADT build needs valid nodeMesh bounding boxes. Defer until a
    // code path actually searches/adds points so topology-only clients (and
    // unit tests) can construct IFaceLink without a full mesh.
    bool pointLocatorReady = false;
    void ValidateGridIndex( const Grid & grid, const char * operation ) const;

public:
    // Builds the PointLocator ADT on first use (requires valid nodeMesh boxes).
    void EnsurePointLocatorReady();
    [[nodiscard]] bool IsPointLocatorReady() const { return pointLocatorReady; }
    // Non-null locator object (may not be Initialize'd yet). Throws if missing.
    [[nodiscard]] PointLocator & GetPointLocator();
    [[nodiscard]] const PointLocator & GetPointLocator() const;
    // Take ownership of a locator (production ADT or test double).
    // markReady=true skips EnsurePointLocatorReady's grid-bound ADT build so
    // topology-only / synthetic-point tests need no nodeMesh boxes.
    // Null locator is rejected.
    void ReplacePointLocator( std::unique_ptr< PointLocator > locator, bool markReady );
    [[nodiscard]] FaceSearch & GetFaceSearch();
    [[nodiscard]] const FaceSearch & GetFaceSearch() const;
    [[nodiscard]] InterfaceFaceSplitTable & MutableInterfaceSplit();
    [[nodiscard]] const InterfaceFaceSplitTable & GetInterfaceSplit() const;
    [[nodiscard]] Grid & GetGrid( int zoneIndex );
    [[nodiscard]] const Grid & GetGrid( int zoneIndex ) const;
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

    // Read-only views of the live dual maps. Prefer these over touching fields
    // when only observing topology (tests and match helpers).
    [[nodiscard]] const LinkField & GetLocalToGlobalMap() const { return l2g; }
    [[nodiscard]] const LinkField & GetGlobalToZoneIds() const { return gI2Zid; }
    [[nodiscard]] const LinkField & GetGlobalToLocalIds() const { return g2l; }

    // Ensure l2g / l2gNew have at least nZones rows (zone-parallel tables).
    // Does not touch global-face tables (gI2Zid / g2l).
    void EnsureInterfaceZoneCount( std::size_t nZones );

    // Ensure l2g[zoneId] can hold localFaceCount slots (resize up only).
    // Zone row must already exist (EnsureInterfaceZoneCount). Does not register
    // global faces; pair with RegisterLocalInterfaceFace for seeding tests.
    void EnsureLocalInterfaceFaceCount( int zoneId, std::size_t localFaceCount );

    // Register (zid, localFaceId) on global face gIid using the live dual maps.
    // Thin wrapper around the free-function helper; keeps fields private.
    void RegisterLocalInterfaceFace( int gIid, bool isNew, int zid, int localFaceId );

    // Split-rebuild bookkeeping: nChild must have one row per live zone (l2g).
    // Call after InitNewLgMapping so l2gNew already has empty per-zone rows.
    void EnsureSplitBookkeepingMatchesZones();

    // Bounds-checked mutable access for the active zone during split rebuild.
    // Prefer these over raw l2gNew[zoneId] / nChild[zoneId] at call sites.
    [[nodiscard]] IntField & MutableL2gNewForZone( int zoneId );
    [[nodiscard]] IntField & MutableNChildForZone( int zoneId );
    [[nodiscard]] const IntField & GetL2gForZone( int zoneId ) const;

    // Zone count for split bookkeeping table (may lag l2g until Ensure*).
    [[nodiscard]] std::size_t GetSplitBookkeepingZoneCount() const { return nChild.size(); }

    // Staging inverse maps (*New). Prefer these over raw gI2ZidNew / g2lNew.
    [[nodiscard]] LinkField & MutableGlobalToZoneIdsNew() { return gI2ZidNew; }
    [[nodiscard]] LinkField & MutableGlobalToLocalIdsNew() { return g2lNew; }
    [[nodiscard]] const LinkField & GetGlobalToZoneIdsNew() const { return gI2ZidNew; }
    [[nodiscard]] const LinkField & GetGlobalToLocalIdsNew() const { return g2lNew; }
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
// isNew: append a new parallel row to gI2Zid/g2l; otherwise append to the
// existing row, which supports at most two zone participants. The existing-row
// path rejects duplicate participants and a third zone before mutating any map.
// Always writes l2g[zid][localFaceId] = gIid on success.
// Named *OnMaps to avoid colliding with IFaceLink::RegisterLocalInterfaceFace.
void RegisterLocalInterfaceFaceOnMaps(
    int gIid,
    bool isNew,
    int zid,
    int localFaceId,
    LinkField & l2g,
    LinkField & gI2Zid,
    LinkField & g2l );

// Select the first participant on a periodic partner face that is not the
// current local face. A periodic lookup can return a face whose mapping row
// contains the current face; never create a self-link in that case.
bool FindPeriodicInterfacePartner(
    int selfZoneId,
    int selfLocalFaceId,
    const IntField & zoneIds,
    const IntField & localFaceIds,
    int & outRemoteZoneId,
    int & outRemoteLocalFaceId );

// Find the first non-self participant on globalFaceId, but only when the
// exact (selfZoneId, selfLocalFaceId) pair is registered. Returns false when
// self is absent or no other participant exists; output arguments stay unchanged.
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

