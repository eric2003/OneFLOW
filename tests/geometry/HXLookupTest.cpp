// tests/geometry/HXLookupTest.cpp
// FaceSolver dedup depends on HXKey normalizing node lists by sort.
// These tests lock that policy without constructing FaceSolver.

#include <gtest/gtest.h>

#include "HXLookup.h"
#include "Constant.h"

using namespace ONEFLOW;

TEST( HXLookupTest, FindOrAddAssignsSequentialIds )
{
    HXLookup< int > lookup;
    auto [ id0, new0 ] = lookup.FindOrAdd( std::vector< int >{ 1, 2, 3 } );
    auto [ id1, new1 ] = lookup.FindOrAdd( std::vector< int >{ 4, 5, 6 } );

    EXPECT_TRUE( new0 );
    EXPECT_TRUE( new1 );
    EXPECT_EQ( id0, 0 );
    EXPECT_EQ( id1, 1 );
    EXPECT_EQ( lookup.Size(), static_cast< std::size_t >( 2 ) );
}

TEST( HXLookupTest, SameNodesDifferentOrderAreSameFace )
{
    // Intent: undirected face identity; order/orientation must not create a
    // second face entry (HXKey::Normalize sorts node ids).
    HXLookup< int > lookup;
    auto [ id0, new0 ] = lookup.FindOrAdd( std::vector< int >{ 3, 1, 2 } );
    auto [ id1, new1 ] = lookup.FindOrAdd( std::vector< int >{ 1, 2, 3 } );
    auto [ id2, new2 ] = lookup.FindOrAdd( std::vector< int >{ 2, 3, 1 } );

    EXPECT_TRUE( new0 );
    EXPECT_FALSE( new1 );
    EXPECT_FALSE( new2 );
    EXPECT_EQ( id0, id1 );
    EXPECT_EQ( id1, id2 );
    EXPECT_EQ( lookup.Size(), static_cast< std::size_t >( 1 ) );
}

TEST( HXLookupTest, FindReturnsInvalidWhenMissing )
{
    HXLookup< int > lookup;
    lookup.FindOrAdd( std::vector< int >{ 1, 2 } );
    EXPECT_EQ( lookup.Find( std::vector< int >{ 9, 8 } ), INVALID_INDEX );
    EXPECT_EQ( lookup.Find( std::vector< int >{ 2, 1 } ), 0 );
}

TEST( HXLookupTest, ClearRemovesAllEntries )
{
    HXLookup< int > lookup;
    lookup.FindOrAdd( std::vector< int >{ 1, 2, 3 } );
    lookup.Clear();
    EXPECT_EQ( lookup.Size(), static_cast< std::size_t >( 0 ) );
    auto [ id, isNew ] = lookup.FindOrAdd( std::vector< int >{ 1, 2, 3 } );
    EXPECT_TRUE( isNew );
    EXPECT_EQ( id, 0 );
}
