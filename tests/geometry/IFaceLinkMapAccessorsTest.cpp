// tests/geometry/IFaceLinkMapAccessorsTest.cpp
// Locks read-only dual-map accessors on IFaceLink without full mesh scan.

#include <gtest/gtest.h>

#include "IFaceLink.h"
#include "GridHandles.h"
#include "HXDefine.h"

using namespace ONEFLOW;

TEST( IFaceLinkMapAccessorsTest, ConstAccessorsReferToLiveTables )
{
    Grids grids;
    IFaceLink link( grids );

    EXPECT_EQ( &link.GetLocalToGlobalMap(), &link.l2g );
    EXPECT_EQ( &link.GetGlobalToZoneIds(), &link.gI2Zid );
    EXPECT_EQ( &link.GetGlobalToLocalIds(), &link.g2l );
}

TEST( IFaceLinkMapAccessorsTest, AccessorsSeeManualTableEdits )
{
    Grids grids;
    IFaceLink link( grids );
    link.l2g.push_back( IntField{ 0 } );
    link.gI2Zid.push_back( IntField{ 0 } );
    link.g2l.push_back( IntField{ 0 } );

    ASSERT_EQ( link.GetLocalToGlobalMap().size(), static_cast< std::size_t >( 1 ) );
    EXPECT_EQ( link.GetLocalToGlobalMap()[ 0 ], ( IntField{ 0 } ) );
    EXPECT_EQ( link.GetGlobalToZoneIds()[ 0 ], ( IntField{ 0 } ) );
    EXPECT_EQ( link.GetGlobalToLocalIds()[ 0 ], ( IntField{ 0 } ) );
}
