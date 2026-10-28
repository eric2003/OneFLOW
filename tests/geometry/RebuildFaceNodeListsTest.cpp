// tests/geometry/RebuildFaceNodeListsTest.cpp
// Locks face-node list rebuild when interface boundary faces split.

#include <gtest/gtest.h>
#include <stdexcept>

#include "FaceTopo.h"
#include "FaceSearch.h"
#include "Boundary.h"
#include "HXDefine.h"

using namespace ONEFLOW;

TEST( RebuildFaceNodeListsTest, NonInterfaceAndInteriorCopied )
{
    LinkField faces;
    faces.push_back( IntField{ 1, 2 } ); // boundary wall
    faces.push_back( IntField{ 3, 4 } ); // interior
    const IntField bcType{ BC::SOLID_SURFACE };
    LinkField out;
    RebuildFaceNodeListsAfterInterfaceSplit(
        faces, bcType, /*nBFaces=*/1, /*nFaces=*/2,
        LinkField{}, LinkField{}, LinkField{},
        false, 0, out );

    ASSERT_EQ( out.size(), static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( out[ 0 ], ( IntField{ 1, 2 } ) );
    EXPECT_EQ( out[ 1 ], ( IntField{ 3, 4 } ) );
}

TEST( RebuildFaceNodeListsTest, UnsplitInterfaceKeepsParentNodes )
{
    LinkField faces;
    faces.push_back( IntField{ 10, 11, 12 } );
    const IntField bcType{ BC::INTERFACE };
    LinkField children( 1 ); // empty children
    LinkField out;
    RebuildFaceNodeListsAfterInterfaceSplit(
        faces, bcType, 1, 1, children, LinkField{}, LinkField{},
        false, 0, out );

    ASSERT_EQ( out.size(), static_cast< std::size_t >( 1 ) );
    EXPECT_EQ( out[ 0 ], ( IntField{ 10, 11, 12 } ) );
}

TEST( RebuildFaceNodeListsTest, SplitInterfaceResolvesChildNodeLists )
{
    // Parent face nodes 100,101,102,103. Child 0 uses local corners 0,1,2 (flag 1).
    LinkField faces;
    faces.push_back( IntField{ 100, 101, 102, 103 } );
    const IntField bcType{ BC::INTERFACE };

    LinkField childGlobals( 1 );
    childGlobals[ 0 ] = IntField{ 0 };

    LinkField childNodeIds( 1 );
    LinkField childNodeFlags( 1 );
    childNodeIds[ 0 ] = IntField{ 0, 1, 2 };
    childNodeFlags[ 0 ] = IntField{ 1, 1, 1 };

    LinkField out;
    RebuildFaceNodeListsAfterInterfaceSplit(
        faces, bcType, 1, 1, childGlobals, childNodeIds, childNodeFlags,
        false, 0, out );

    ASSERT_EQ( out.size(), static_cast< std::size_t >( 1 ) );
    EXPECT_EQ( out[ 0 ], ( IntField{ 100, 101, 102 } ) );
}

TEST( RebuildFaceNodeListsTest, SharedChildUsesEachParentLocalMapping )
{
    // Both parents refer to global child 7, but use different local corner indices.
    LinkField faces{
        IntField{ 100, 101, 102, 103 },
        IntField{ 101, 102, 100, 104 }
    };
    const IntField bcType{ BC::INTERFACE, BC::INTERFACE };
    const LinkField childGlobals{ IntField{ 7 }, IntField{ 7 } };
    const std::vector< LinkField > childNodeIds{
        LinkField{ IntField{ 0, 1, 2 } },
        LinkField{ IntField{ 2, 0, 1 } }
    };
    const std::vector< LinkField > childNodeFlags{
        LinkField{ IntField{ 1, 1, 1 } },
        LinkField{ IntField{ 1, 1, 1 } }
    };

    LinkField out;
    RebuildFaceNodeListsAfterInterfaceSplit(
        faces, bcType, 2, 2, childGlobals, childNodeIds, childNodeFlags,
        false, 0, out );

    ASSERT_EQ( out.size(), static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( out[ 0 ], ( IntField{ 100, 101, 102 } ) );
    EXPECT_EQ( out[ 1 ], ( IntField{ 100, 101, 102 } ) );
}

TEST( RebuildFaceNodeListsTest, RejectsInterfaceChildListShortfall )
{
    LinkField faces;
    faces.push_back( IntField{ 1, 2 } );
    const IntField bcType{ BC::INTERFACE };
    LinkField out;
    EXPECT_THROW(
        RebuildFaceNodeListsAfterInterfaceSplit(
            faces, bcType, 1, 1, LinkField{}, LinkField{}, LinkField{},
            false, 0, out ),
        std::out_of_range );
}

TEST( RebuildFaceNodeListsTest, MappingCollectorDoesNotPublishPartialOutputs )
{
    // The second parent has malformed parent-specific rows. Collection must
    // fail without exposing mappings collected from the first parent.
    LinkField zoneToLocalToGlobal{ IntField{ 0, 1 } };
    InterfaceFaceSplitTable split;
    split.childFaceIdsOfParent = LinkField{
        IntField{ 0 },
        IntField{ 1 }
    };
    split.childFaceNodeIdsByParent.resize( 2 );
    split.childFaceNodeFlagsByParent.resize( 2 );
    split.childFaceNodeIdsByParent[ 0 ] = LinkField{ IntField{ 0, 1 } };
    split.childFaceNodeFlagsByParent[ 0 ] = LinkField{ IntField{ 1, 1 } };
    split.childFaceNodeIdsByParent[ 1 ] = LinkField{
        IntField{ 0, 1 },
        IntField{ 1, 2 }
    };
    split.childFaceNodeFlagsByParent[ 1 ] = LinkField{
        IntField{ 1, 1 },
        IntField{ 1, 1 }
    };

    std::vector< LinkField > nodeIds{ LinkField{ IntField{ 90 } } };
    std::vector< LinkField > nodeFlags{ LinkField{ IntField{ 91 } } };
    const auto originalNodeIds = nodeIds;
    const auto originalNodeFlags = nodeFlags;

    EXPECT_THROW(
        CollectBoundaryInterfaceChildNodeMappings(
            IntField{ BC::INTERFACE, BC::INTERFACE },
            2,
            zoneToLocalToGlobal,
            split,
            0,
            nodeIds,
            nodeFlags ),
        std::logic_error );

    EXPECT_EQ( nodeIds, originalNodeIds );
    EXPECT_EQ( nodeFlags, originalNodeFlags );
}

TEST( RebuildFaceNodeListsTest, SplitProducerCollectorAndRebuildPreserveSharedChildMapping )
{
    // Two parent rings share one global child, but need different local corner rows.
    FaceSearch search;
    const int firstParent = search.AddFace( IntField{ 100, 101, 102, 103 } );
    const int secondParent = search.AddFace( IntField{ 101, 102, 103, 104 } );
    const int sharedChild = search.AddFace( IntField{ 101, 102, 103 } );
    search.interfaceSplit.ResetForFaceCount( search.Size() );

    search.SplitQuad2Tri( firstParent );
    search.SplitQuad2Tri( secondParent );

    ASSERT_EQ(
        search.interfaceSplit.childFaceIdsOfParent[ firstParent ],
        ( IntField{ sharedChild } ) );
    ASSERT_EQ(
        search.interfaceSplit.childFaceIdsOfParent[ secondParent ],
        ( IntField{ sharedChild } ) );
    ASSERT_EQ(
        search.interfaceSplit.childFaceNodeIdsByParent[ firstParent ][ 0 ],
        ( IntField{ 1, 2, 3 } ) );
    ASSERT_EQ(
        search.interfaceSplit.childFaceNodeIdsByParent[ secondParent ][ 0 ],
        ( IntField{ 0, 1, 2 } ) );

    LinkField zoneToLocalToGlobal{
        IntField{ firstParent },
        IntField{ secondParent }
    };
    const IntField interfaceBc{ BC::INTERFACE };

    LinkField firstChildIds;
    CollectBoundaryInterfaceChildFaceIds(
        interfaceBc, 1, zoneToLocalToGlobal, search.interfaceSplit, 0, firstChildIds );
    std::vector< LinkField > firstNodeIds, firstNodeFlags;
    CollectBoundaryInterfaceChildNodeMappings(
        interfaceBc, 1, zoneToLocalToGlobal, search.interfaceSplit, 0,
        firstNodeIds, firstNodeFlags );

    LinkField secondChildIds;
    CollectBoundaryInterfaceChildFaceIds(
        interfaceBc, 1, zoneToLocalToGlobal, search.interfaceSplit, 1, secondChildIds );
    std::vector< LinkField > secondNodeIds, secondNodeFlags;
    CollectBoundaryInterfaceChildNodeMappings(
        interfaceBc, 1, zoneToLocalToGlobal, search.interfaceSplit, 1,
        secondNodeIds, secondNodeFlags );

    LinkField firstOut, secondOut;
    RebuildFaceNodeListsAfterInterfaceSplit(
        LinkField{ IntField{ 100, 101, 102, 103 } }, interfaceBc, 1, 1,
        firstChildIds, firstNodeIds, firstNodeFlags, false, 0, firstOut );
    RebuildFaceNodeListsAfterInterfaceSplit(
        LinkField{ IntField{ 101, 102, 103, 104 } }, interfaceBc, 1, 1,
        secondChildIds, secondNodeIds, secondNodeFlags, false, 0, secondOut );

    ASSERT_EQ( firstOut.size(), static_cast< std::size_t >( 1 ) );
    ASSERT_EQ( secondOut.size(), static_cast< std::size_t >( 1 ) );
    EXPECT_EQ( firstOut[ 0 ], ( IntField{ 101, 102, 103 } ) );
    EXPECT_EQ( secondOut[ 0 ], ( IntField{ 101, 102, 103 } ) );
}


TEST( RebuildFaceNodeListsTest, ParentSpecificMappingTakesPrecedenceOverLegacyRows )
{
    // The legacy row is child-global; it cannot represent both parents safely.
    // Collection must use the row owned by this parent when it is available.
    LinkField zoneToLocalToGlobal{ IntField{ 0 } };
    InterfaceFaceSplitTable split;
    split.childFaceIdsOfParent = LinkField{ IntField{ 1 } };
    split.childFaceNodeIds = LinkField{ IntField{ 99, 98, 97 } };
    split.childFaceNodeFlags = LinkField{ IntField{ 0, 0, 0 } };
    split.childFaceNodeIdsByParent = std::vector< LinkField >{
        LinkField{ IntField{ 2, 0, 1 } },
        LinkField{}
    };
    split.childFaceNodeFlagsByParent = std::vector< LinkField >{
        LinkField{ IntField{ 1, 1, 1 } },
        LinkField{}
    };

    std::vector< LinkField > nodeIds;
    std::vector< LinkField > nodeFlags;
    CollectBoundaryInterfaceChildNodeMappings(
        IntField{ BC::INTERFACE },
        1,
        zoneToLocalToGlobal,
        split,
        0,
        nodeIds,
        nodeFlags );

    ASSERT_EQ( nodeIds.size(), static_cast< std::size_t >( 1 ) );
    ASSERT_EQ( nodeFlags.size(), static_cast< std::size_t >( 1 ) );
    EXPECT_EQ( nodeIds[ 0 ], ( LinkField{ IntField{ 2, 0, 1 } } ) );
    EXPECT_EQ( nodeFlags[ 0 ], ( LinkField{ IntField{ 1, 1, 1 } } ) );
}

TEST( RebuildFaceNodeListsTest, FailedRebuildDoesNotPublishPartialFaceLists )
{
    LinkField faces{
        IntField{ 10, 11, 12, 13 },
        IntField{ 20, 21, 22, 23 }
    };
    const IntField bcType{ BC::INTERFACE, BC::INTERFACE };
    const LinkField childGlobals{ IntField{ 7 }, IntField{ 8 } };
    const std::vector< LinkField > childNodeIds{
        LinkField{ IntField{ 0, 1, 2 } },
        LinkField{ IntField{ 0, 1, 4 } }
    };
    const std::vector< LinkField > childNodeFlags{
        LinkField{ IntField{ 1, 1, 1 } },
        LinkField{ IntField{ 1, 1, 1 } }
    };
    LinkField out{ IntField{ 90, 91, 92 } };
    const LinkField originalOut = out;

    EXPECT_THROW(
        RebuildFaceNodeListsAfterInterfaceSplit(
            faces, bcType, 2, 2, childGlobals, childNodeIds, childNodeFlags,
            false, 0, out ),
        std::out_of_range );

    EXPECT_EQ( out, originalOut );
}

TEST( RebuildFaceNodeListsTest, LegacyMappingAdapterUsesCanonicalRebuildAtomically )
{
    LinkField faces{
        IntField{ 10, 11, 12, 13 },
        IntField{ 20, 21, 22, 23 }
    };
    const IntField bcType{ BC::INTERFACE, BC::INTERFACE };
    const LinkField childGlobals{ IntField{ 0 }, IntField{ 1 } };
    const LinkField globalNodeIds{
        IntField{ 0, 1, 2 },
        IntField{ 0, 1, 4 }
    };
    const LinkField globalNodeFlags{
        IntField{ 1, 1, 1 },
        IntField{ 1, 1, 1 }
    };
    LinkField out{ IntField{ 90, 91, 92 } };
    const LinkField originalOut = out;

    EXPECT_THROW(
        RebuildFaceNodeListsAfterInterfaceSplit(
            faces, bcType, 2, 2, childGlobals, globalNodeIds, globalNodeFlags,
            false, 0, out ),
        std::out_of_range );

    EXPECT_EQ( out, originalOut );
}
