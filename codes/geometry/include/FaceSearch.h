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
#include <set>
#include <vector>


BeginNameSpace( ONEFLOW )

// Unique face registry: identity via HXLookup (sorted nodes); first-seen
// node ring kept privately for orientation-sensitive consumers.
class FaceSearchBasic
{
public:
    FaceSearchBasic() = default;
    ~FaceSearchBasic() = default;

    // Add a face. Returns the unique face ID (always >= 0).
    // Same node set in any order yields the same id.
    int AddFace( const IntField & faceNode );

    // Only find. Returns INVALID_INDEX if not found (order-invariant).
    int FindFace( const IntField & faceNode ) const;

    // First-seen node ring for faceId (bounds-checked). Not overwritten on
    // later AddFace of the same identity.
    const IntField & GetFaceNodeRing( int faceId ) const;

    std::size_t Size() const { return lookup_.Size(); }
    void Clear();

private:
    HXLookup<int> lookup_;   // responsible for uniqueness and ID assignment
    // Original node order of each face (index == faceId).
    HXVector< IntField > faceNodeRings;
};

class IFaceLink;

// A registered child face stays paired with its original candidate position.
struct RegisteredChildFace
{
    int globalFaceId;
    int candidateIndex;
};

// Resolve candidate node rings to registered global face IDs.
// Missing candidates are skipped; each result retains its original candidate index.
std::vector< RegisteredChildFace > ResolveRegisteredChildFaces(
    const FaceSearchBasic & faceSearch,
    const LinkField & candidateFaces );

// Tables produced when interface faces are split (e.g. quad->tri).
// Separate from FaceSearchBasic (unique face registry).
struct InterfaceFaceSplitTable
{
    LinkField childFaceIdsOfParent;
    LinkField childFaceNodeIds;
    LinkField childFaceNodeFlags;
    std::vector< LinkField > childFaceNodeIdsByParent;
    std::vector< LinkField > childFaceNodeFlagsByParent;

    void Clear()
    {
        childFaceIdsOfParent.clear();
        childFaceNodeIds.clear();
        childFaceNodeFlags.clear();
        childFaceNodeIdsByParent.clear();
        childFaceNodeFlagsByParent.clear();
    }

    // The first three tables use global face IDs; parent-specific rows use the
    // global parent face ID and align by child position within that parent's list.
    void ResetForFaceCount( std::size_t faceCount )
    {
        childFaceIdsOfParent.assign( faceCount, IntField{} );
        childFaceNodeIds.assign( faceCount, IntField{} );
        childFaceNodeFlags.assign( faceCount, IntField{} );
        childFaceNodeIdsByParent.assign( faceCount, LinkField{} );
        childFaceNodeFlagsByParent.assign( faceCount, LinkField{} );
    }

    void ClearParentMapping( int parentFaceId );
    // Reject conflicting metadata for a shared child instead of silently overwriting it.
    void SetChildMapping(
        int parentFaceId,
        int childFaceId,
        const IntField & localNodeIds,
        const IntField & localNodeFlags );
};

class FaceSearch : public FaceSearchBasic
{
public:
    FaceSearch() = default;
    ~FaceSearch() = default;

    InterfaceFaceSplitTable interfaceSplit;

    // Clear both the unique-face registry and its derived split mappings.
    void Clear();
    void CalcNewFaceId( IFaceLink & iFaceLink );
    void SplitQuad2Tri( int faceId );
    // Split a 2-node edge face using mid-point insertion (needs IFaceLink geometry).
    void SplitLine( IFaceLink & iFaceLink, int globalFaceId );

private:
    // Try to split edge (nodeId size 2) at midpoint if the midpoint already
    // exists in the point locator. globalFaceId indexes the dual-map global face tables (GetGlobalToZoneIds).
    bool TrySplitLineAtMidpoint(
        IFaceLink & iFaceLink,
        int globalFaceId,
        const IntField & nodeId,
        LinkField & localLineId,
        LinkField & localLineFlag,
        std::vector< RegisteredChildFace > & registeredChildren );
};

// Build directed global child-edge candidates. Outputs are replaced and order
// follows the parent edge direction.
// @doc-algorithm line-interface-split
// @doc-stage 2
// @doc-pure
// @doc-formula m = (1/2) * (a + b) component-wise in R^3
// Arithmetic midpoint of a segment. No mesh, locator, or tolerance.
void ComputeSegmentMidpoint(
    Real x0, Real y0, Real z0,
    Real x1, Real y1, Real z1,
    Real & mx, Real & my, Real & mz );

void BuildGlobalLineSplitCandidates(
    int firstNode,
    int secondNode,
    int midpointGlobalNode,
    LinkField & globalLineIds );

// @doc-algorithm line-interface-split
// @doc-stage 1
// @doc-pure
// Line-split policy: the global interface face must reference exactly one
// zone. Returns that zone ID or throws.
int RequireSoleZoneIdForGlobalFace(
    const LinkField & zoneIdsPerGlobalFace,
    int globalFaceId );

// @doc-algorithm line-interface-split
// @doc-stage 7
// @doc-pure
// Build local child-edge indices and endpoint flags after a midpoint is found.
// The local midpoint index is zone-specific; global node identity is handled
// separately by BuildGlobalLineSplitCandidates.
void BuildLocalLineSplitMapping(
    int midpointLocalNode,
    LinkField & localLineId,
    LinkField & localLineFlag );

// Pure helpers for quad diagonal split (testable without IFaceLink).
// BuildQuadDiagonalLocalTris clears outputs then fills four candidate tris:
// diagonal 0-2 (tris 0-1-2 and 2-3-0) and diagonal 1-3 (tris 3-0-1 and 1-2-3).
void BuildQuadDiagonalLocalTris( LinkField & localTriId, LinkField & localTriFlag );

// Map local corner indices in localTriId through nodeId into global node ids.
void MapLocalTrisToGlobalNodes(
    const IntField & nodeId,
    const LinkField & localTriId,
    LinkField & triId );

EndNameSpace
