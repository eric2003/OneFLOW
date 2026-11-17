// tests/geometry/CollectBoundaryInterfaceChildCountsTest.cpp
// Locks boundary-order interface child counting used by SetNewFace2Cell /
// ResetNumberOfBoundaryCondition.

#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

#include "FaceTopo.h"
#include "FaceSearch.h"
#include "Boundary.h"
#include "HXDefine.h"

using namespace ONEFLOW;

namespace {

// zone 0: two local interface parents -> global 0,1
// parent 0 has two children; parent 1 has none
void Seed( LinkField & l2g, InterfaceFaceSplitTable & split )
{
    l2g.clear();
    l2g.push_back( IntField{ 0, 1 } );
    split.Clear();
    split.childFaceIdsOfParent.resize( 2 );
    split.childFaceIdsOfParent[ 0 ] = IntField{ 0, 1 };
    split.childFaceIdsOfParent[ 1 ] = IntField{};
    split.childFaceNodeIds.resize( 2 );
    split.childFaceNodeFlags.resize( 2 );
    split.childFaceNodeIds[ 0 ] = IntField{ 1, 2, 3 };
    split.childFaceNodeFlags[ 0 ] = IntField{ 0, 0, 0 };
    split.childFaceNodeIds[ 1 ] = IntField{ 2, 3, 4 };
    split.childFaceNodeFlags[ 1 ] = IntField{ 0, 0, 0 };
}

} // namespace

TEST( CollectBoundaryInterfaceChildNodeMappingsTest, RejectsMissingParentRowInsteadOfUsingSharedChildGlobalRow )
{
    LinkField l2g{ IntField{ 0, 1 } };
    InterfaceFaceSplitTable split;
    split.ResetForFaceCount( 3 );
    split.childFaceIdsOfParent[ 0 ] = IntField{ 2 };
    split.childFaceIdsOfParent[ 1 ] = IntField{ 2 };
    split.SetChildMapping(
        0, 2, IntField{ 0, 1, 2 }, IntField{ 1, 1, 1 } );

    std::vector< LinkField > nodeIds{ LinkField{ IntField{ 90 } } };
    std::vector< LinkField > nodeFlags{ LinkField{ IntField{ -90 } } };

    EXPECT_THROW(
        CollectBoundaryInterfaceChildNodeMappings(
            IntField{ BC::INTERFACE, BC::INTERFACE },
            2,
            l2g,
            split,
            0,
            nodeIds,
            nodeFlags ),
        std::logic_error );

    EXPECT_EQ( nodeIds, ( std::vector< LinkField >{ LinkField{ IntField{ 90 } } } ) );
    EXPECT_EQ( nodeFlags, ( std::vector< LinkField >{ LinkField{ IntField{ -90 } } } ) );
}


TEST( CollectBoundaryInterfaceChildCountsTest, SkipsNonInterfaceBoundaryFaces )
{
    LinkField l2g;
    InterfaceFaceSplitTable split;
    Seed( l2g, split );

    const IntField bcType{ BC::SOLID_SURFACE, BC::INTERFACE, BC::INFLOW, BC::INTERFACE };
    IntField counts;
    CollectBoundaryInterfaceChildCounts(
        bcType, /*nBFaces=*/4, l2g, split, /*zoneId=*/0, counts );

    EXPECT_EQ( counts, ( IntField{ 2, 0 } ) );
}

TEST( CollectBoundaryInterfaceChildCountsTest, EmptyWhenNoInterface )
{
    LinkField l2g;
    InterfaceFaceSplitTable split;
    Seed( l2g, split );
    const IntField bcType{ BC::SOLID_SURFACE, BC::FARFIELD };
    IntField counts;
    CollectBoundaryInterfaceChildCounts(
        bcType, 2, l2g, split, 0, counts );
    EXPECT_TRUE( counts.empty() );
}

TEST( CollectBoundaryInterfaceChildCountsTest, RejectsInvalidBoundaryCount )
{
    LinkField l2g;
    InterfaceFaceSplitTable split;
    Seed( l2g, split );
    const IntField bcType{ BC::INTERFACE };
    IntField counts;
    EXPECT_THROW(
        CollectBoundaryInterfaceChildCounts( bcType, 5, l2g, split, 0, counts ),
        std::logic_error );
}


TEST( CollectBoundaryInterfaceChildFaceIdsTest, ReturnsChildIdListsInOrder )
{
    LinkField l2g;
    InterfaceFaceSplitTable split;
    Seed( l2g, split );

    const IntField bcType{ BC::SOLID_SURFACE, BC::INTERFACE, BC::INTERFACE };
    LinkField lists;
    CollectBoundaryInterfaceChildFaceIds(
        bcType, 3, l2g, split, 0, lists );

    ASSERT_EQ( lists.size(), static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( lists[ 0 ], ( IntField{ 0, 1 } ) );
    EXPECT_TRUE( lists[ 1 ].empty() );
}

TEST( CollectBoundaryInterfaceChildCountsTest, CountsMatchFaceIdListSizes )
{
    LinkField l2g;
    InterfaceFaceSplitTable split;
    Seed( l2g, split );

    const IntField bcType{ BC::INTERFACE, BC::INTERFACE };
    LinkField lists;
    IntField counts;
    CollectBoundaryInterfaceChildFaceIds( bcType, 2, l2g, split, 0, lists );
    CollectBoundaryInterfaceChildCounts( bcType, 2, l2g, split, 0, counts );

    ASSERT_EQ( counts.size(), lists.size() );
    EXPECT_EQ( counts[ 0 ], static_cast< int >( lists[ 0 ].size() ) );
    EXPECT_EQ( counts[ 1 ], static_cast< int >( lists[ 1 ].size() ) );
}
