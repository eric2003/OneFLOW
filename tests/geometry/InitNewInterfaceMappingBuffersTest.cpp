// tests/geometry/InitNewInterfaceMappingBuffersTest.cpp
// Locks InitNewLgMapping buffer setup without IFaceLink.

#include <gtest/gtest.h>
#include <stdexcept>

#include "IFaceLink.h"
#include "HXDefine.h"

using namespace ONEFLOW;

TEST( InitNewInterfaceMappingBuffersTest, CopiesInverseMapsAndResizesL2gNew )
{
    LinkField gI2Zid;
    LinkField g2l;
    gI2Zid.push_back( IntField{ 0, 1 } );
    g2l.push_back( IntField{ 2, 3 } );

    LinkField l2g;
    l2g.push_back( IntField{ 0 } );
    l2g.push_back( IntField{ 0 } );
    l2g.push_back( IntField{ 0 } ); // three zones

    LinkField gI2ZidNew, g2lNew, l2gNew;
    // Pre-fill to ensure clear/resize
    l2gNew.push_back( IntField{ 99 } );

    InitNewInterfaceMappingBuffers(
        gI2Zid, g2l, l2g, gI2ZidNew, g2lNew, l2gNew );

    ASSERT_EQ( gI2ZidNew.size(), static_cast< std::size_t >( 1 ) );
    EXPECT_EQ( gI2ZidNew[ 0 ], ( IntField{ 0, 1 } ) );
    EXPECT_EQ( g2lNew[ 0 ], ( IntField{ 2, 3 } ) );
    ASSERT_EQ( l2gNew.size(), static_cast< std::size_t >( 3 ) );
    EXPECT_TRUE( l2gNew[ 0 ].empty() );
    EXPECT_TRUE( l2gNew[ 1 ].empty() );
    EXPECT_TRUE( l2gNew[ 2 ].empty() );
}

TEST( InitNewInterfaceMappingBuffersTest, RejectsMismatchedInverseMaps )
{
    LinkField gI2Zid;
    LinkField g2l;
    gI2Zid.push_back( IntField{ 0 } );
    // g2l empty
    LinkField l2g;
    LinkField a, b, c;
    EXPECT_THROW(
        InitNewInterfaceMappingBuffers( gI2Zid, g2l, l2g, a, b, c ),
        std::logic_error );
}
