// tests/geometry/FaceSearchBasicTest.cpp
// FaceSearchBasic pairs undirected identity (HXLookup) with directed storage
// (faceNodeRings keeps the node order of the first AddFace).

#include <gtest/gtest.h>
#include <stdexcept>

#include "FaceSearch.h"
#include "Constant.h"
#include "HXDefine.h"

using namespace ONEFLOW;

TEST( FaceSearchBasicTest, AddFaceAssignsSequentialIds )
{
    FaceSearchBasic search;
    EXPECT_EQ( search.AddFace( IntField{ 1, 2, 3 } ), 0 );
    EXPECT_EQ( search.AddFace( IntField{ 4, 5, 6 } ), 1 );
    EXPECT_EQ( search.Size(), static_cast< std::size_t >( 2 ) );
    ASSERT_EQ( search.Size(), static_cast< std::size_t >( 2 ) );
}

TEST( FaceSearchBasicTest, DuplicateNodesSameIdKeepsFirstOrder )
{
    // Intent: identity ignores order; orientation memory keeps first ring.
    FaceSearchBasic search;
    const int id0 = search.AddFace( IntField{ 1, 2, 3 } );
    const int id1 = search.AddFace( IntField{ 3, 1, 2 } );

    EXPECT_EQ( id0, id1 );
    EXPECT_EQ( search.Size(), static_cast< std::size_t >( 1 ) );
    ASSERT_EQ( search.Size(), static_cast< std::size_t >( 1 ) );
    EXPECT_EQ( search.GetFaceNodeRing( 0 ), ( IntField{ 1, 2, 3 } ) );
}

TEST( FaceSearchBasicTest, FindFaceIsOrderInvariant )
{
    FaceSearchBasic search;
    search.AddFace( IntField{ 2, 0, 1 } );
    EXPECT_EQ( search.FindFace( IntField{ 0, 1, 2 } ), 0 );
    EXPECT_EQ( search.FindFace( IntField{ 9, 8, 7 } ), INVALID_INDEX );
}

TEST( FaceSearchBasicTest, ClearResetsLookupAndFaceArray )
{
    FaceSearchBasic search;
    search.AddFace( IntField{ 1, 2 } );
    search.Clear();
    EXPECT_EQ( search.Size(), static_cast< std::size_t >( 0 ) );
    EXPECT_EQ( search.Size(), static_cast< std::size_t >( 0 ) );
    EXPECT_EQ( search.AddFace( IntField{ 1, 2 } ), 0 );
}

TEST( FaceSearchBasicTest, ClearResetsDerivedSplitMappings )
{
    // Intent: clearing a FaceSearch invalidates mappings derived from its face IDs.
    FaceSearch search;
    ASSERT_EQ( search.AddFace( IntField{ 1, 2, 3, 4 } ), 0 );

    search.interfaceSplit.childFaceIdsOfParent = LinkField{ IntField{ 4, 5 } };
    search.interfaceSplit.childFaceNodeIds = LinkField{ IntField{ 0, 1, 2 } };
    search.interfaceSplit.childFaceNodeFlags = LinkField{ IntField{ 1, -1, 1 } };

    search.Clear();

    EXPECT_EQ( search.Size(), static_cast< std::size_t >( 0 ) );
    EXPECT_TRUE( search.interfaceSplit.childFaceIdsOfParent.empty() );
    EXPECT_TRUE( search.interfaceSplit.childFaceNodeIds.empty() );
    EXPECT_TRUE( search.interfaceSplit.childFaceNodeFlags.empty() );
    EXPECT_TRUE( search.interfaceSplit.childFaceNodeIdsByParent.empty() );
    EXPECT_TRUE( search.interfaceSplit.childFaceNodeFlagsByParent.empty() );
    EXPECT_EQ( search.AddFace( IntField{ 1, 2 } ), 0 );
}

TEST( FaceSearchBasicTest, NonQuadSplitClearsStaleParentMapping )
{
    FaceSearch search;
    ASSERT_EQ( search.AddFace( IntField{ 1, 2, 3 } ), 0 );
    search.interfaceSplit.ResetForFaceCount( search.Size() );
    search.interfaceSplit.childFaceIdsOfParent[ 0 ] = IntField{ 7, 8 };

    search.SplitQuad2Tri( 0 );

    EXPECT_TRUE( search.interfaceSplit.childFaceIdsOfParent[ 0 ].empty() );
}

TEST( FaceSearchBasicTest, RepeatedQuadSplitReplacesParentChildMapping )
{
    // Intent: rebuilding one parent's split mapping must not append duplicate child IDs.
    FaceSearch search;
    ASSERT_EQ( search.AddFace( IntField{ 10, 20, 30, 40 } ), 0 );
    EXPECT_EQ( search.AddFace( IntField{ 10, 20, 30 } ), 1 );
    EXPECT_EQ( search.AddFace( IntField{ 30, 40, 10 } ), 2 );
    EXPECT_EQ( search.AddFace( IntField{ 40, 10, 20 } ), 3 );
    EXPECT_EQ( search.AddFace( IntField{ 20, 30, 40 } ), 4 );
    search.interfaceSplit.ResetForFaceCount( search.Size() );

    search.SplitQuad2Tri( 0 );
    search.SplitQuad2Tri( 0 );

    ASSERT_EQ( search.interfaceSplit.childFaceIdsOfParent[ 0 ].size(),
               static_cast< std::size_t >( 4 ) );
    EXPECT_EQ( search.interfaceSplit.childFaceIdsOfParent[ 0 ],
               ( IntField{ 1, 2, 3, 4 } ) );
}

TEST( FaceSearchBasicTest, QuadSplitUsesOnlyRegisteredTrianglesForOneDiagonal )
{
    FaceSearch search;
    ASSERT_EQ( search.AddFace( IntField{ 10, 20, 30, 40 } ), 0 );

    // Register only the two triangles using diagonal 10-30.
    EXPECT_EQ( search.AddFace( IntField{ 10, 20, 30 } ), 1 );
    EXPECT_EQ( search.AddFace( IntField{ 30, 40, 10 } ), 2 );

    search.interfaceSplit.ResetForFaceCount( search.Size() );
    search.SplitQuad2Tri( 0 );

    ASSERT_EQ( search.interfaceSplit.childFaceIdsOfParent[ 0 ].size(),
               static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( search.interfaceSplit.childFaceIdsOfParent[ 0 ],
               ( IntField{ 1, 2 } ) );
}

TEST( FaceSearchBasicTest, QuadSplitUsesOriginalCandidateIndexForPartialMatch )
{
    // Only candidate 1 is registered; its local mapping must not shift to candidate 0.
    FaceSearch search;
    ASSERT_EQ( search.AddFace( IntField{ 10, 20, 30, 40 } ), 0 );
    ASSERT_EQ( search.AddFace( IntField{ 30, 40, 10 } ), 1 );

    search.interfaceSplit.ResetForFaceCount( search.Size() );
    search.SplitQuad2Tri( 0 );

    ASSERT_EQ( search.interfaceSplit.childFaceIdsOfParent[ 0 ], ( IntField{ 1 } ) );
    ASSERT_EQ( search.interfaceSplit.childFaceNodeIdsByParent[ 0 ].size(),
               static_cast< std::size_t >( 1 ) );
    EXPECT_EQ( search.interfaceSplit.childFaceNodeIdsByParent[ 0 ][ 0 ],
               ( IntField{ 2, 3, 0 } ) );
    EXPECT_EQ( search.interfaceSplit.childFaceNodeFlagsByParent[ 0 ][ 0 ],
               ( IntField{ 1, 1, 1 } ) );
}

TEST( FaceSearchBasicTest, QuadSplitPreservesLastCandidateIndexForPartialMatch )
{
    // Three missing candidates must not compact the final candidate's local map.
    FaceSearch search;
    ASSERT_EQ( search.AddFace( IntField{ 10, 20, 30, 40 } ), 0 );
    ASSERT_EQ( search.AddFace( IntField{ 20, 30, 40 } ), 1 );

    search.interfaceSplit.ResetForFaceCount( search.Size() );
    search.SplitQuad2Tri( 0 );

    ASSERT_EQ( search.interfaceSplit.childFaceIdsOfParent[ 0 ], ( IntField{ 1 } ) );
    ASSERT_EQ( search.interfaceSplit.childFaceNodeIdsByParent[ 0 ].size(),
               static_cast< std::size_t >( 1 ) );
    EXPECT_EQ( search.interfaceSplit.childFaceNodeIdsByParent[ 0 ][ 0 ],
               ( IntField{ 1, 2, 3 } ) );
    EXPECT_EQ( search.interfaceSplit.childFaceNodeFlagsByParent[ 0 ][ 0 ],
               ( IntField{ 1, 1, 1 } ) );
}

TEST( FaceSearchBasicTest, QuadSplitMatchesRegisteredTrianglesRegardlessOfNodeOrder )
{
    FaceSearch search;
    ASSERT_EQ( search.AddFace( IntField{ 10, 20, 30, 40 } ), 0 );

    // Face identity is node-set based, but the first registered ring is preserved.
    const IntField firstTriangleRing{ 30, 20, 10 };
    const IntField secondTriangleRing{ 10, 40, 30 };
    EXPECT_EQ( search.AddFace( firstTriangleRing ), 1 );
    EXPECT_EQ( search.AddFace( secondTriangleRing ), 2 );

    search.interfaceSplit.ResetForFaceCount( search.Size() );
    search.SplitQuad2Tri( 0 );

    ASSERT_EQ( search.interfaceSplit.childFaceIdsOfParent[ 0 ].size(),
               static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( search.interfaceSplit.childFaceIdsOfParent[ 0 ],
               ( IntField{ 1, 2 } ) );
    EXPECT_EQ( search.GetFaceNodeRing( 1 ), firstTriangleRing );
    EXPECT_EQ( search.GetFaceNodeRing( 2 ), secondTriangleRing );
}

TEST( FaceSearchBasicTest, InvalidQuadCandidateLeavesMappedOutputUnchanged )
{
    // Intent: validate every candidate before publishing any mapped triangles.
    const IntField quadNodes{ 10, 20, 30, 40 };
    const LinkField localTriangles{
        IntField{ 0, 1, 2 },
        IntField{ 2, 3, 4 }
    };
    LinkField mappedTriangles{ IntField{ 90, 91, 92 } };

    EXPECT_THROW(
        MapLocalTrisToGlobalNodes( quadNodes, localTriangles, mappedTriangles ),
        std::out_of_range );

    EXPECT_EQ( mappedTriangles, ( LinkField{ IntField{ 90, 91, 92 } } ) );
}

TEST( FaceSearchBasicTest, RegisteredChildResolutionReturnsOnlyKnownChildren )
{
    FaceSearchBasic search;
    EXPECT_EQ( search.AddFace( IntField{ 10, 20 } ), 0 );
    EXPECT_EQ( search.AddFace( IntField{ 20, 30 } ), 1 );

    const LinkField candidates{
        IntField{ 10, 20 },
        IntField{ 20, 30 },
        IntField{ 30, 40 }
    };

    const std::vector< RegisteredChildFace > children =
        ResolveRegisteredChildFaces( search, candidates );

    ASSERT_EQ( children.size(), static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( children[ 0 ].globalFaceId, 0 );
    EXPECT_EQ( children[ 0 ].candidateIndex, 0 );
    EXPECT_EQ( children[ 1 ].globalFaceId, 1 );
    EXPECT_EQ( children[ 1 ].candidateIndex, 1 );
}

TEST( FaceSearchBasicTest, RegisteredChildResolutionPreservesCandidateIndexForPartialMatch )
{
    // Intent: a missing first candidate must not shift the local mapping for a registered second child.
    FaceSearchBasic search;
    EXPECT_EQ( search.AddFace( IntField{ 20, 30 } ), 0 );

    const LinkField candidates{
        IntField{ 10, 20 },
        IntField{ 20, 30 }
    };

    const std::vector< RegisteredChildFace > children =
        ResolveRegisteredChildFaces( search, candidates );

    ASSERT_EQ( children.size(), static_cast< std::size_t >( 1 ) );
    EXPECT_EQ( children[ 0 ].globalFaceId, 0 );
    EXPECT_EQ( children[ 0 ].candidateIndex, 1 );
}

TEST( FaceSearchBasicTest, RegisteredChildResolutionReturnsEmptyWhenNoChildExists )
{
    FaceSearchBasic search;
    search.AddFace( IntField{ 10, 20 } );

    const LinkField candidates{
        IntField{ 20, 30 },
        IntField{ 30, 40 }
    };

    const std::vector< RegisteredChildFace > children =
        ResolveRegisteredChildFaces( search, candidates );

    EXPECT_TRUE( children.empty() );
}

TEST( FaceSearchBasicTest, RegisteredChildResolutionRejectsDuplicateIdentity )
{
    FaceSearchBasic search;
    search.AddFace( IntField{ 10, 20 } );

    const LinkField candidates{
        IntField{ 10, 20 },
        IntField{ 20, 10 }
    };

    EXPECT_THROW(
        ResolveRegisteredChildFaces( search, candidates ),
        std::logic_error );
}

TEST( FaceSearchBasicTest, GlobalLineCandidatesPreserveEdgeDirection )
{
    // Global child identity follows the directed parent ring A -> B.
    LinkField candidates{ IntField{ 99, 98 } };

    BuildGlobalLineSplitCandidates( 12, 18, 25, candidates );

    EXPECT_EQ( candidates, ( LinkField{ IntField{ 12, 25 }, IntField{ 25, 18 } } ) );
}

TEST( FaceSearchBasicTest, GlobalLineCandidatesPreserveReversedEdgeDirection )
{
    LinkField candidates;

    BuildGlobalLineSplitCandidates( 18, 12, 25, candidates );

    EXPECT_EQ( candidates, ( LinkField{ IntField{ 18, 25 }, IntField{ 25, 12 } } ) );
}

TEST( FaceSearchBasicTest, GlobalLineCandidatesReplacePreviousOutput )
{
    LinkField candidates{ IntField{ 95, 94 } };

    BuildGlobalLineSplitCandidates( 3, 4, 8, candidates );

    EXPECT_EQ( candidates, ( LinkField{ IntField{ 3, 8 }, IntField{ 8, 4 } } ) );
}

TEST( FaceSearchBasicTest, LocalLineSplitMappingPreservesEndpointFlags )
{
    LinkField localLineId{ IntField{ 9, 9 } };
    LinkField localLineFlag{ IntField{ 9, 9 } };

    BuildLocalLineSplitMapping( 7, localLineId, localLineFlag );

    EXPECT_EQ( localLineId, ( LinkField{ IntField{ 0, 7 }, IntField{ 7, 1 } } ) );
    EXPECT_EQ( localLineFlag, ( LinkField{ IntField{ 1, -1 }, IntField{ -1, 1 } } ) );
}

TEST( FaceSearchBasicTest, LocalLineSplitMappingReplacesPreviousOutput )
{
    LinkField localLineId{ IntField{ 99, 98 } };
    LinkField localLineFlag{ IntField{ 97, 96 } };

    BuildLocalLineSplitMapping( 2, localLineId, localLineFlag );

    EXPECT_EQ( localLineId, ( LinkField{ IntField{ 0, 2 }, IntField{ 2, 1 } } ) );
    EXPECT_EQ( localLineFlag, ( LinkField{ IntField{ 1, -1 }, IntField{ -1, 1 } } ) );
}

TEST( FaceSearchBasicTest, GetFaceNodeRingRejectsOutOfRange )
{
    FaceSearchBasic search;
    EXPECT_THROW( search.GetFaceNodeRing( 0 ), std::out_of_range );
    search.AddFace( IntField{ 1, 2, 3 } );
    EXPECT_NO_THROW( search.GetFaceNodeRing( 0 ) );
    EXPECT_THROW( search.GetFaceNodeRing( 1 ), std::out_of_range );
}

TEST( FaceSearchBasicTest, SharedChildKeepsParentSpecificMappingRows )
{
    // Intent: global child identity is shared, but local rows belong to each parent.
    InterfaceFaceSplitTable split;
    split.ResetForFaceCount( 3 );
    split.childFaceIdsOfParent[ 0 ] = IntField{ 2 };
    split.childFaceIdsOfParent[ 1 ] = IntField{ 2 };

    split.SetChildMapping( 0, 2, IntField{ 0, 1, 2 }, IntField{ 1, 1, 1 } );
    split.SetChildMapping( 1, 2, IntField{ 2, 0, 1 }, IntField{ -1, 1, 1 } );

    ASSERT_EQ( split.childFaceNodeIdsByParent.size(), static_cast< std::size_t >( 3 ) );
    EXPECT_EQ( split.childFaceNodeIdsByParent[ 0 ][ 0 ], ( IntField{ 0, 1, 2 } ) );
    EXPECT_EQ( split.childFaceNodeIdsByParent[ 1 ][ 0 ], ( IntField{ 2, 0, 1 } ) );
    EXPECT_EQ( split.childFaceNodeFlagsByParent[ 0 ][ 0 ], ( IntField{ 1, 1, 1 } ) );
    EXPECT_EQ( split.childFaceNodeFlagsByParent[ 1 ][ 0 ], ( IntField{ -1, 1, 1 } ) );

    // Legacy child-global metadata is compatibility-only and must not be
    // overwritten by a later parent's different local-corner mapping or flags.
    EXPECT_EQ( split.childFaceNodeIds[ 2 ], ( IntField{ 0, 1, 2 } ) );
    EXPECT_EQ( split.childFaceNodeFlags[ 2 ], ( IntField{ 1, 1, 1 } ) );
}

TEST( FaceSearchBasicTest, ResetForFaceCountClearsAllSplitTables )
{
    // Intent: every table uses the global face ID domain and reset removes stale rows.
    InterfaceFaceSplitTable split;
    split.ResetForFaceCount( 3 );
    split.childFaceIdsOfParent[ 0 ] = IntField{ 1, 2 };
    split.childFaceNodeIds[ 1 ] = IntField{ 0, 1 };
    split.childFaceNodeFlags[ 1 ] = IntField{ 1, -1 };
    split.childFaceNodeIdsByParent[ 0 ] = LinkField{ IntField{ 0, 1 } };
    split.childFaceNodeFlagsByParent[ 0 ] = LinkField{ IntField{ 1, -1 } };

    split.ResetForFaceCount( 2 );

    ASSERT_EQ( split.childFaceIdsOfParent.size(), static_cast< std::size_t >( 2 ) );
    ASSERT_EQ( split.childFaceNodeIds.size(), static_cast< std::size_t >( 2 ) );
    ASSERT_EQ( split.childFaceNodeFlags.size(), static_cast< std::size_t >( 2 ) );
    EXPECT_TRUE( split.childFaceIdsOfParent[ 0 ].empty() );
    EXPECT_TRUE( split.childFaceNodeIds[ 1 ].empty() );
    EXPECT_TRUE( split.childFaceNodeFlags[ 1 ].empty() );
    ASSERT_EQ( split.childFaceNodeIdsByParent.size(), static_cast< std::size_t >( 2 ) );
    ASSERT_EQ( split.childFaceNodeFlagsByParent.size(), static_cast< std::size_t >( 2 ) );
    EXPECT_TRUE( split.childFaceNodeIdsByParent[ 0 ].empty() );
    EXPECT_TRUE( split.childFaceNodeFlagsByParent[ 0 ].empty() );
}

TEST( FaceSearchBasicTest, ClearParentMappingOnlyClearsRequestedParent )
{
    InterfaceFaceSplitTable split;
    // Face IDs 0 and 1 are parents; face ID 2 is their shared child.
    split.ResetForFaceCount( 3 );
    split.childFaceIdsOfParent[ 0 ] = IntField{ 2 };
    split.childFaceIdsOfParent[ 1 ] = IntField{ 2 };
    split.childFaceNodeIds[ 2 ] = IntField{ 0, 1, 2 };
    split.childFaceNodeFlags[ 2 ] = IntField{ 1, 1, 1 };
    split.childFaceNodeIdsByParent[ 0 ] = LinkField{ IntField{ 1, 2, 3 } };
    split.childFaceNodeFlagsByParent[ 0 ] = LinkField{ IntField{ 1, 1, 1 } };
    split.childFaceNodeIdsByParent[ 1 ] = LinkField{ IntField{ 0, 1, 2 } };
    split.childFaceNodeFlagsByParent[ 1 ] = LinkField{ IntField{ 1, 1, 1 } };

    split.ClearParentMapping( 0 );

    EXPECT_TRUE( split.childFaceIdsOfParent[ 0 ].empty() );
    EXPECT_EQ( split.childFaceIdsOfParent[ 1 ], ( IntField{ 2 } ) );
    EXPECT_TRUE( split.childFaceNodeIdsByParent[ 0 ].empty() );
    EXPECT_TRUE( split.childFaceNodeFlagsByParent[ 0 ].empty() );
    EXPECT_EQ( split.childFaceNodeIdsByParent[ 1 ], ( LinkField{ IntField{ 0, 1, 2 } } ) );
    EXPECT_EQ( split.childFaceNodeFlagsByParent[ 1 ], ( LinkField{ IntField{ 1, 1, 1 } } ) );

    // The global child metadata is shared and must survive clearing one parent.
    EXPECT_EQ( split.childFaceNodeIds[ 2 ], ( IntField{ 0, 1, 2 } ) );
    EXPECT_EQ( split.childFaceNodeFlags[ 2 ], ( IntField{ 1, 1, 1 } ) );
}

TEST( FaceSearchBasicTest, ClearParentMappingRejectsInvalidParentId )
{
    InterfaceFaceSplitTable split;
    split.ResetForFaceCount( 1 );

    EXPECT_THROW( split.ClearParentMapping( -1 ), std::out_of_range );
    EXPECT_THROW( split.ClearParentMapping( 1 ), std::out_of_range );
}

TEST( FaceSearchBasicTest, ClearParentMappingRejectsMisalignedParentTables )
{
    InterfaceFaceSplitTable split;
    split.ResetForFaceCount( 2 );
    split.childFaceIdsOfParent[ 0 ] = IntField{ 1 };
    split.childFaceNodeFlagsByParent.clear();

    EXPECT_THROW( split.ClearParentMapping( 0 ), std::logic_error );
    EXPECT_EQ( split.childFaceIdsOfParent[ 0 ], ( IntField{ 1 } ) );
}

TEST( FaceSearchBasicTest, QuadSplitRejectsDuplicateChildBeforeAppendingIt )
{
    // Intent: a degenerate quad must not register the same global child twice.
    FaceSearch search;
    ASSERT_EQ( search.AddFace( IntField{ 10, 20, 30, 30 } ), 0 );
    ASSERT_EQ( search.AddFace( IntField{ 10, 20, 30 } ), 1 );
    search.interfaceSplit.ResetForFaceCount( search.Size() );

    EXPECT_THROW( search.SplitQuad2Tri( 0 ), std::logic_error );

    // Duplicate candidates are rejected before any child mapping is published.
    EXPECT_TRUE( search.interfaceSplit.childFaceIdsOfParent[ 0 ].empty() );
    EXPECT_TRUE( search.interfaceSplit.childFaceNodeIdsByParent[ 0 ].empty() );
    EXPECT_TRUE( search.interfaceSplit.childFaceNodeFlagsByParent[ 0 ].empty() );
}

TEST( FaceSearchBasicTest, SetChildMappingRejectsDuplicateChildRegistration )
{
    // Intent: each parent child position must identify exactly one global child face.
    InterfaceFaceSplitTable split;
    split.ResetForFaceCount( 3 );
    split.childFaceIdsOfParent[ 0 ] = IntField{ 2, 2 };

    EXPECT_THROW(
        split.SetChildMapping( 0, 2, IntField{ 0, 1, 2 }, IntField{ 1, 1, 1 } ),
        std::logic_error );

    EXPECT_TRUE( split.childFaceNodeIdsByParent[ 0 ].empty() );
    EXPECT_TRUE( split.childFaceNodeFlagsByParent[ 0 ].empty() );
}

TEST( FaceSearchBasicTest, SetChildMappingRejectsConflictingMappingForSameParent )
{
    // Intent: one parent-child relation cannot change its local corner mapping silently.
    InterfaceFaceSplitTable split;
    split.ResetForFaceCount( 3 );
    split.childFaceIdsOfParent[ 0 ] = IntField{ 2 };

    split.SetChildMapping( 0, 2, IntField{ 0, 1, 2 }, IntField{ 1, 1, 1 } );
    EXPECT_NO_THROW(
        split.SetChildMapping( 0, 2, IntField{ 0, 1, 2 }, IntField{ 1, 1, 1 } ) );

    EXPECT_THROW(
        split.SetChildMapping( 0, 2, IntField{ 2, 1, 0 }, IntField{ 1, 1, 1 } ),
        std::logic_error );

    ASSERT_EQ( split.childFaceNodeIdsByParent[ 0 ].size(), static_cast< std::size_t >( 1 ) );
    EXPECT_EQ( split.childFaceNodeIdsByParent[ 0 ][ 0 ], ( IntField{ 0, 1, 2 } ) );
    EXPECT_EQ( split.childFaceNodeFlagsByParent[ 0 ][ 0 ], ( IntField{ 1, 1, 1 } ) );
}

TEST( FaceSearchBasicTest, SetChildMappingRejectsMisalignedGlobalRowsWithoutMutation )
{
    InterfaceFaceSplitTable split;
    split.ResetForFaceCount( 3 );
    split.childFaceIdsOfParent[ 0 ] = IntField{ 2 };
    split.childFaceNodeIds[ 2 ] = IntField{ 9, 8 };
    split.childFaceNodeFlags[ 2 ] = IntField{ 1 };

    EXPECT_THROW(
        split.SetChildMapping( 0, 2, IntField{ 0, 1 }, IntField{ 1, 1 } ),
        std::logic_error );

    EXPECT_TRUE( split.childFaceNodeIdsByParent[ 0 ].empty() );
    EXPECT_TRUE( split.childFaceNodeFlagsByParent[ 0 ].empty() );
    EXPECT_EQ( split.childFaceNodeIds[ 2 ], ( IntField{ 9, 8 } ) );
    EXPECT_EQ( split.childFaceNodeFlags[ 2 ], ( IntField{ 1 } ) );
}


TEST( FaceSearchBasicTest, SetChildMappingRejectsMissingParentSpecificRow )
{
    InterfaceFaceSplitTable split;
    split.ResetForFaceCount( 2 );
    split.childFaceIdsOfParent[ 0 ] = IntField{ 1 };
    split.childFaceNodeFlagsByParent.clear();

    EXPECT_THROW(
        split.SetChildMapping( 0, 1, IntField{ 0, 1 }, IntField{ 1, 1 } ),
        std::logic_error );
}

TEST( FaceSearchBasicTest, SetChildMappingRejectsMisalignedParentRowsWithoutMutation )
{
    InterfaceFaceSplitTable split;
    split.ResetForFaceCount( 3 );
    split.childFaceIdsOfParent[ 0 ] = IntField{ 1, 2 };
    split.childFaceNodeIdsByParent[ 0 ] = LinkField{ IntField{ 0, 1 } };

    EXPECT_THROW(
        split.SetChildMapping( 0, 2, IntField{ 1, 2 }, IntField{ 1, 1 } ),
        std::logic_error );

    EXPECT_EQ(
        split.childFaceNodeIdsByParent[ 0 ],
        ( LinkField{ IntField{ 0, 1 } } ) );
    EXPECT_TRUE( split.childFaceNodeFlagsByParent[ 0 ].empty() );
}

