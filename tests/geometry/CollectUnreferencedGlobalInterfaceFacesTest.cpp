// tests/geometry/CollectUnreferencedGlobalInterfaceFacesTest.cpp
// Locks detection of inverse-map rows not referenced by any l2g entry.

#include <gtest/gtest.h>
#include <stdexcept>

#include "FaceTopo.h"
#include "HXDefine.h"

using namespace ONEFLOW;

TEST( CollectUnreferencedGlobalInterfaceFacesTest, EmptyWhenAllReferenced )
{
    LinkField l2g;
    l2g.push_back( IntField{ 0, 1 } );
    l2g.push_back( IntField{ 0, 1 } );
    IntField out;
    CollectUnreferencedGlobalInterfaceFaces( l2g, /*globalFaceRowCount=*/2, out );
    EXPECT_TRUE( out.empty() );
}

TEST( CollectUnreferencedGlobalInterfaceFacesTest, FindsOrphanParentRow )
{
    LinkField l2g;
    l2g.push_back( IntField{ 2, 3, 1 } );
    l2g.push_back( IntField{ 2, 3, 1 } );
    IntField out;
    CollectUnreferencedGlobalInterfaceFaces( l2g, /*globalFaceRowCount=*/4, out );
    EXPECT_EQ( out, ( IntField{ 0 } ) );
}

TEST( CollectUnreferencedGlobalInterfaceFacesTest, RejectsNegativeOrOor )
{
    LinkField l2g;
    l2g.push_back( IntField{ -1 } );
    IntField out{ 9 };
    EXPECT_THROW(
        CollectUnreferencedGlobalInterfaceFaces( l2g, 2, out ),
        std::out_of_range );
    // Output must be unchanged when validation fails.
    EXPECT_EQ( out, ( IntField{ 9 } ) );

    l2g[ 0 ] = IntField{ 2 };
    EXPECT_THROW(
        CollectUnreferencedGlobalInterfaceFaces( l2g, 2, out ),
        std::out_of_range );
    EXPECT_EQ( out, ( IntField{ 9 } ) );
}

TEST( CollectUnreferencedGlobalInterfaceFacesTest, ZeroRowCount )
{
    LinkField l2g;
    l2g.push_back( IntField{} );
    IntField out{ 1, 2 };
    CollectUnreferencedGlobalInterfaceFaces( l2g, 0, out );
    EXPECT_TRUE( out.empty() );
}
