// tests/geometry/MultiZoneInterfaceSplitE2ETest.cpp
// Table-level multi-zone interface split (no CGNS, no IFaceLink required).
// Calls RunSymmetricInterfaceSplitMappingPipeline — the same entry a future
// workbench/GUI should call with editable tables.

#include <gtest/gtest.h>

#include "FaceTopo.h"
#include "IFaceLink.h"
#include "HXDefine.h"
#include "Boundary.h"

using namespace ONEFLOW;

namespace {

void SeedTwoZoneTwoFaceMaps( LinkField & l2g, LinkField & gI2Zid, LinkField & g2l )
{
    l2g.assign( 2, IntField() );
    l2g[ 0 ].assign( 2, -1 );
    l2g[ 1 ].assign( 2, -1 );
    gI2Zid.clear();
    g2l.clear();
    RegisterLocalInterfaceFaceOnMaps( 0, true,  0, 0, l2g, gI2Zid, g2l );
    RegisterLocalInterfaceFaceOnMaps( 0, false, 1, 0, l2g, gI2Zid, g2l );
    RegisterLocalInterfaceFaceOnMaps( 1, true,  0, 1, l2g, gI2Zid, g2l );
    RegisterLocalInterfaceFaceOnMaps( 1, false, 1, 1, l2g, gI2Zid, g2l );
}

} // namespace

TEST( MultiZoneInterfaceSplitE2ETest, SplitOneSharedParentViaTablePipeline )
{
    LinkField l2g, gI2Zid, g2l;
    SeedTwoZoneTwoFaceMaps( l2g, gI2Zid, g2l );

    EXPECT_EQ( l2g[ 0 ], ( IntField{ 0, 1 } ) );
    EXPECT_EQ( l2g[ 1 ], ( IntField{ 0, 1 } ) );
    EXPECT_EQ( gI2Zid[ 0 ], ( IntField{ 0, 1 } ) );
    EXPECT_EQ( g2l[ 0 ], ( IntField{ 0, 0 } ) );

    const IntField parents{ 0, 1 };
    LinkField children( 2 );
    children[ 0 ] = IntField{ 2, 3 };
    children[ 1 ] = IntField{};

    const IntField bcType{ BC::INTERFACE, BC::INTERFACE };
    const IntField bcNameId{ 100, 200 };
    const IntField childCounts{ 2, 0 };

    InterfaceSplitMappingPipelineResult out;
    RunSymmetricInterfaceSplitMappingPipeline(
        l2g, gI2Zid, g2l,
        parents, children,
        bcType, bcNameId, childCounts,
        out );

    EXPECT_EQ( out.l2g[ 0 ], ( IntField{ 2, 3, 1 } ) );
    EXPECT_EQ( out.l2g[ 1 ], ( IntField{ 2, 3, 1 } ) );
    EXPECT_EQ( out.nChild[ 0 ], ( IntField{ 2, 0, 0, 0 } ) );
    EXPECT_EQ( out.gI2Zid[ 2 ], ( IntField{ 0, 1 } ) );
    EXPECT_EQ( out.g2l[ 2 ], ( IntField{ 0, 0 } ) );
    EXPECT_EQ( out.gI2Zid[ 3 ], ( IntField{ 0, 1 } ) );
    EXPECT_EQ( out.g2l[ 3 ], ( IntField{ 1, 1 } ) );
    EXPECT_EQ( out.gI2Zid[ 1 ], ( IntField{ 0, 1 } ) );
    EXPECT_EQ( out.g2l[ 1 ], ( IntField{ 2, 2 } ) );
    EXPECT_EQ( out.outBcType, ( IntField{ BC::INTERFACE, BC::INTERFACE, BC::INTERFACE } ) );
    EXPECT_EQ( out.outBcNameId, ( IntField{ 100, 100, 200 } ) );

    // Legacy InitNew copy: retired parent global 0 remains as an inverse row.
    ASSERT_GE( out.gI2Zid.size(), static_cast< std::size_t >( 1 ) );
    EXPECT_EQ( out.gI2Zid[ 0 ], ( IntField{ 0, 1 } ) );

    IntField unreferenced;
    CollectUnreferencedGlobalInterfaceFaces(
        out.l2g, out.gI2Zid.size(), unreferenced );
    EXPECT_EQ( unreferenced, ( IntField{ 0 } ) );
}

TEST( MultiZoneInterfaceSplitE2ETest, FailedPipelineDoesNotPublishPartialResult )
{
    LinkField l2g, gI2Zid, g2l;
    SeedTwoZoneTwoFaceMaps( l2g, gI2Zid, g2l );

    const IntField parents{ 0, 1 };
    LinkField children( 2 );
    children[ 0 ] = IntField{ 2, 3 };
    children[ 1 ] = IntField{};
    const IntField bcType{ BC::INTERFACE, BC::INTERFACE };
    const IntField bcNameId{ 100, 200 };
    // Only one count is supplied for two interface BC entries.
    const IntField childCounts{ 2 };

    InterfaceSplitMappingPipelineResult out;
    out.l2g = LinkField{ IntField{ 70 } };
    out.gI2Zid = LinkField{ IntField{ 80 } };
    out.g2l = LinkField{ IntField{ 90 } };
    out.nChild = LinkField{ IntField{ 100 } };
    out.outBcType = IntField{ 110 };
    out.outBcNameId = IntField{ 120 };

    EXPECT_THROW(
        RunSymmetricInterfaceSplitMappingPipeline(
            l2g, gI2Zid, g2l,
            parents, children,
            bcType, bcNameId, childCounts,
            out ),
        std::out_of_range );

    EXPECT_EQ( out.l2g, ( LinkField{ IntField{ 70 } } ) );
    EXPECT_EQ( out.gI2Zid, ( LinkField{ IntField{ 80 } } ) );
    EXPECT_EQ( out.g2l, ( LinkField{ IntField{ 90 } } ) );
    EXPECT_EQ( out.nChild, ( LinkField{ IntField{ 100 } } ) );
    EXPECT_EQ( out.outBcType, ( IntField{ 110 } ) );
    EXPECT_EQ( out.outBcNameId, ( IntField{ 120 } ) );
}

