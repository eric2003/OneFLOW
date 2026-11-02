// tests/geometry/InitNewInterfaceMappingBuffersTest.cpp
// Locks InitNewLgMapping / UpdateLgMapping buffer helpers without IFaceLink.

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
    l2g.push_back( IntField{ 0 } );

    LinkField gI2ZidNew, g2lNew, l2gNew;
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
    LinkField l2g;
    LinkField a, b, c;
    EXPECT_THROW(
        InitNewInterfaceMappingBuffers( gI2Zid, g2l, l2g, a, b, c ),
        std::logic_error );
}

TEST( CommitInterfaceMappingBuffersTest, SwapsValidatedNewTablesIntoLive )
{
    LinkField gI2Zid;
    LinkField g2l;
    LinkField l2g;
    gI2Zid.push_back( IntField{ 0 } );
    g2l.push_back( IntField{ 0 } );
    l2g.push_back( IntField{ 0 } );

    LinkField gI2ZidNew;
    LinkField g2lNew;
    LinkField l2gNew;
    gI2ZidNew.push_back( IntField{ 0, 1 } );
    g2lNew.push_back( IntField{ 0, 1 } );
    l2gNew.push_back( IntField{ 0 } );
    l2gNew.push_back( IntField{ 0 } );

    CommitInterfaceMappingBuffers(
        gI2Zid, g2l, l2g, gI2ZidNew, g2lNew, l2gNew, "test" );

    ASSERT_EQ( gI2Zid.size(), static_cast< std::size_t >( 1 ) );
    EXPECT_EQ( gI2Zid[ 0 ], ( IntField{ 0, 1 } ) );
    EXPECT_EQ( g2l[ 0 ], ( IntField{ 0, 1 } ) );
    ASSERT_EQ( l2g.size(), static_cast< std::size_t >( 2 ) );
}

TEST( CommitInterfaceMappingBuffersTest, RejectsInconsistentNewInverseMaps )
{
    LinkField gI2Zid, g2l, l2g;
    LinkField gI2ZidNew, g2lNew, l2gNew;
    gI2ZidNew.push_back( IntField{ 0, 1 } );
    g2lNew.push_back( IntField{ 0 } ); // length mismatch
    EXPECT_THROW(
        CommitInterfaceMappingBuffers(
            gI2Zid, g2l, l2g, gI2ZidNew, g2lNew, l2gNew, "test" ),
        std::logic_error );
}

TEST( CommitInterfaceMappingBuffersTest, FailedValidationLeavesLiveTablesUnchanged )
{
    LinkField gI2Zid{ IntField{ 4 } };
    LinkField g2l{ IntField{ 7 } };
    LinkField l2g{ IntField{ 2, 3 } };

    LinkField gI2ZidNew{ IntField{ 0, 1 } };
    LinkField g2lNew{ IntField{ 9 } };
    LinkField l2gNew{ IntField{ 5 } };

    EXPECT_THROW(
        CommitInterfaceMappingBuffers(
            gI2Zid, g2l, l2g, gI2ZidNew, g2lNew, l2gNew, "test" ),
        std::logic_error );

    EXPECT_EQ( gI2Zid, ( LinkField{ IntField{ 4 } } ) );
    EXPECT_EQ( g2l, ( LinkField{ IntField{ 7 } } ) );
    EXPECT_EQ( l2g, ( LinkField{ IntField{ 2, 3 } } ) );

    EXPECT_EQ( gI2ZidNew, ( LinkField{ IntField{ 0, 1 } } ) );
    EXPECT_EQ( g2lNew, ( LinkField{ IntField{ 9 } } ) );
    EXPECT_EQ( l2gNew, ( LinkField{ IntField{ 5 } } ) );
}
