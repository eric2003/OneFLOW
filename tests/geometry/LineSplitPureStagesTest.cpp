// tests/geometry/LineSplitPureStagesTest.cpp
// Composes the pure stages of a line interface split without CGNS, ADT, or
// zone nodeMesh. Production still needs PointLocator + nodeMesh after these
// steps; this test locks the order and table contracts of the pure prefix.

#include <gtest/gtest.h>

#include "FaceSearch.h"
#include "HXDefine.h"

using namespace ONEFLOW;

TEST( LineSplitPureStagesTest, SoleZoneMidpointCandidatesAndResolve )
{
    // Stage 1: sole-zone policy on the dual-map inverse row.
    LinkField zoneIdsPerGlobalFace;
    zoneIdsPerGlobalFace.push_back( IntField{ 3 } );
    EXPECT_EQ( RequireSoleZoneIdForGlobalFace( zoneIdsPerGlobalFace, 0 ), 3 );

    // Stage 2: arithmetic midpoint (no tolerance).
    Real mx = 0, my = 0, mz = 0;
    ComputeSegmentMidpoint( 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, mx, my, mz );
    EXPECT_DOUBLE_EQ( mx, 1.0 );
    EXPECT_DOUBLE_EQ( my, 0.0 );
    EXPECT_DOUBLE_EQ( mz, 0.0 );

    // Production would FindPoint(mx,my,mz) on the locator here (hard edge).
    // Assume global mid node id 10 already exists.
    const int p0 = 4;
    const int p1 = 5;
    const int mid = 10;

    // Stage 4: directed child edge candidates.
    LinkField candidates;
    BuildGlobalLineSplitCandidates( p0, p1, mid, candidates );
    ASSERT_EQ( candidates.size(), static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( candidates[ 0 ], ( IntField{ 4, 10 } ) );
    EXPECT_EQ( candidates[ 1 ], ( IntField{ 10, 5 } ) );

    // Stage 5: only registered faces survive resolution.
    FaceSearchBasic registry;
    const int childA = registry.AddFace( IntField{ 4, 10 } );
    const int childB = registry.AddFace( IntField{ 10, 5 } );
    const std::vector< RegisteredChildFace > resolved =
        ResolveRegisteredChildFaces( registry, candidates );
    ASSERT_EQ( resolved.size(), static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( resolved[ 0 ].globalFaceId, childA );
    EXPECT_EQ( resolved[ 0 ].candidateIndex, 0 );
    EXPECT_EQ( resolved[ 1 ].globalFaceId, childB );
    EXPECT_EQ( resolved[ 1 ].candidateIndex, 1 );

    // Stage 7 (local mapping only; zone mesh insert is the other hard edge).
    LinkField localIds, localFlags;
    BuildLocalLineSplitMapping( /*midpointLocalNode=*/7, localIds, localFlags );
    EXPECT_EQ( localIds[ 0 ], ( IntField{ 0, 7 } ) );
    EXPECT_EQ( localIds[ 1 ], ( IntField{ 7, 1 } ) );
}

TEST( LineSplitPureStagesTest, MissingRegisteredChildAbortsPrefix )
{
    LinkField candidates;
    BuildGlobalLineSplitCandidates( 1, 2, 3, candidates );
    FaceSearchBasic registry;
    registry.AddFace( IntField{ 1, 3 } ); // only first child registered
    const std::vector< RegisteredChildFace > resolved =
        ResolveRegisteredChildFaces( registry, candidates );
    ASSERT_EQ( resolved.size(), static_cast< std::size_t >( 1 ) );
    EXPECT_EQ( resolved[ 0 ].candidateIndex, 0 );
    // Production TrySplitLineAtMidpoint requires non-empty resolve; a single
    // surviving child is allowed by Resolve* but SplitLine still only maps
    // what was resolved. Empty resolve is the hard abort.
    FaceSearchBasic emptyRegistry;
    EXPECT_TRUE( ResolveRegisteredChildFaces( emptyRegistry, candidates ).empty() );
}
