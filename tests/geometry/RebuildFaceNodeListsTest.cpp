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

TEST( RebuildFaceNodeListsTest, SplitProducerCollectorAndRebuildPreserveSharedChildMapping )
{
    // Two parent rings share one global child, but need different local corner rows.
    FaceSearch search;
    const int firstParent = search.AddFace( IntField{ 100, 101, 102, 103 } );
    const int secondParent = search.AddFace( IntField{ 101, 102, 103, 100 } );
    const int sharedChild = search.AddFace( IntField{ 100, 101, 102 } );
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
        ( IntField{ 0, 1, 2 } ) );
    ASSERT_EQ(
        search.interfaceSplit.childFaceNodeIdsByParent[ secondParent ][ 0 ],
        ( IntField{ 3, 0, 1 } ) );

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
        LinkField{ IntField{ 101, 102, 103, 100 } }, interfaceBc, 1, 1,
        secondChildIds, secondNodeIds, secondNodeFlags, false, 0, secondOut );

    ASSERT_EQ( firstOut.size(), static_cast< std::size_t >( 1 ) );
    ASSERT_EQ( secondOut.size(), static_cast< std::size_t >( 1 ) );
    EXPECT_EQ( firstOut[ 0 ], ( IntField{ 100, 101, 102 } ) );
    EXPECT_EQ( secondOut[ 0 ], ( IntField{ 100, 101, 102 } ) );
}

