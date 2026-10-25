// tests/geometry/BcClassificationTest.cpp
// Locks OneFLOW BC classification predicates and BcRecord counts that depend
// on them (interface / wall). No mesh required.

#include <gtest/gtest.h>

#include "Boundary.h"
#include "BcRecord.h"

using namespace ONEFLOW;

TEST( BcClassificationTest, InterfaceIsNegativeOrExplicitInterfaceCode )
{
    // Intent: BC::INTERFACE (-1) and any negative code count as interface.
    // Positive physical BCs do not.
    EXPECT_TRUE( BC::IsInterfaceBc( BC::INTERFACE ) );
    EXPECT_TRUE( BC::IsInterfaceBc( -5 ) );
    EXPECT_FALSE( BC::IsInterfaceBc( BC::SOLID_SURFACE ) );
    EXPECT_FALSE( BC::IsInterfaceBc( BC::INFLOW ) );
    EXPECT_FALSE( BC::IsInterfaceBc( BC::NO_BOUNDARY ) );
}

TEST( BcClassificationTest, WallIsSolidSurfaceOnly )
{
    EXPECT_TRUE( BC::IsWallBc( BC::SOLID_SURFACE ) );
    EXPECT_FALSE( BC::IsWallBc( BC::FARFIELD ) );
    EXPECT_FALSE( BC::IsWallBc( BC::INTERFACE ) );
}

TEST( BcClassificationTest, PoleFamilyIncludesPoleAndPoleN )
{
    EXPECT_TRUE( BC::IsPoleBc( BC::POLE ) );
    EXPECT_TRUE( BC::IsPoleBc( BC::POLE1 ) );
    EXPECT_TRUE( BC::IsPoleBc( BC::POLE2 ) );
    EXPECT_TRUE( BC::IsPoleBc( BC::POLE3 ) );
    EXPECT_FALSE( BC::IsPoleBc( BC::SOLID_SURFACE ) );
}

TEST( BcClassificationTest, NotNormalCombinesInterfaceAndPole )
{
    EXPECT_TRUE( BC::IsNotNormalBc( BC::INTERFACE ) );
    EXPECT_TRUE( BC::IsNotNormalBc( BC::POLE1 ) );
    EXPECT_FALSE( BC::IsNotNormalBc( BC::INFLOW ) );
}

TEST( BcClassificationTest, SlipfaceIsGeneric2 )
{
    EXPECT_TRUE( BC::IsSlipfaceBc( BC::GENERIC_2 ) );
    EXPECT_FALSE( BC::IsSlipfaceBc( BC::GENERIC_1 ) );
}

TEST( BcRecordCountTest, CalcNIFaceCountsInterfaceTypes )
{
    BcRecord record;
    record.bcType = IntField{ BC::SOLID_SURFACE, BC::INTERFACE, -3, BC::INFLOW };
    record.bcNameId = IntField{ 0, 1, 2, 3 };

    EXPECT_EQ( record.GetBoundaryFaceCount(), 4 );
    EXPECT_EQ( record.CalcNIFace(), 2 ); // INTERFACE and -3
}

TEST( BcRecordCountTest, CalcNumWallFaceCountsSolidSurface )
{
    BcRecord record;
    record.bcType = IntField{ BC::SOLID_SURFACE, BC::FARFIELD, BC::SOLID_SURFACE };
    record.bcNameId = IntField{ 0, 1, 2 };

    EXPECT_EQ( record.CalcNumWallFace(), 2 );
}
