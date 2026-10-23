// tests/geometry/BcTypeMapTest.cpp
// Locks CGNS -> OneFLOW BC type mapping without GridElem or a full mesh.

#include <gtest/gtest.h>

#include "Boundary.h"
#include "HXCgns.h"
#include "HXDefine.h"

using namespace ONEFLOW;

TEST( BcTypeMapTest, MapFromCgnsInPlaceConvertsKnownTypes )
{
    BcTypeMap map;
    map.Init();

    IntField types{
        static_cast< int >( CGNS_ENUMV( BCWall ) ),
        static_cast< int >( CGNS_ENUMV( BCInflow ) ),
        static_cast< int >( CGNS_ENUMV( BCFarfield ) ),
        static_cast< int >( CGNS_ENUMV( BCTypeNull ) )
    };

    map.MapFromCgnsInPlace( types );

    EXPECT_EQ( types[ 0 ], BC::SOLID_SURFACE );
    EXPECT_EQ( types[ 1 ], BC::INFLOW );
    EXPECT_EQ( types[ 2 ], BC::FARFIELD );
    EXPECT_EQ( types[ 3 ], BC::INTERFACE );
}

TEST( BcTypeMapTest, UnknownCgnsTypeMapsToGeneric2 )
{
    BcTypeMap map;
    map.Init();

    // Arbitrary id not present in the table.
    EXPECT_EQ( map.Cgns2OneFlow( 999999 ), BC::GENERIC_2 );
}

TEST( BcTypeMapTest, MapFromCgnsInPlaceEmptyIsNoOp )
{
    BcTypeMap map;
    map.Init();
    IntField empty;
    EXPECT_NO_THROW( map.MapFromCgnsInPlace( empty ) );
    EXPECT_TRUE( empty.empty() );
}

TEST( BcTypeMapTest, RoundTripCommonOneFlowTypes )
{
    BcTypeMap map;
    map.Init();

    EXPECT_EQ( map.Cgns2OneFlow( map.OneFlow2Cgns( BC::INFLOW ) ), BC::INFLOW );
    EXPECT_EQ( map.Cgns2OneFlow( map.OneFlow2Cgns( BC::OUTFLOW ) ), BC::OUTFLOW );
    EXPECT_EQ( map.Cgns2OneFlow( map.OneFlow2Cgns( BC::SOLID_SURFACE ) ), BC::SOLID_SURFACE );
}
