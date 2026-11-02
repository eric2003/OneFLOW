// tests/geometry/RebuildZoneL2gNewTest.cpp
// Locks l2gNew / nChild rebuild after interface parent faces split.

#include <gtest/gtest.h>
#include <stdexcept>

#include "FaceTopo.h"
#include "HXDefine.h"

using namespace ONEFLOW;

TEST( RebuildZoneL2gNewTest, UnsplitParentsCopyGlobalIds )
{
    const IntField parents{ 10, 20 };
    LinkField children( 2 ); // both empty
    IntField l2gNew;
    IntField nChild;

    RebuildZoneL2gNewAfterInterfaceSplit( parents, children, l2gNew, nChild );

    EXPECT_EQ( l2gNew, ( IntField{ 10, 20 } ) );
    EXPECT_EQ( nChild, ( IntField{ 0, 0 } ) );
}

TEST( RebuildZoneL2gNewTest, SplitParentEmitsChildrenAndChildSlots )
{
    const IntField parents{ 5 };
    LinkField children( 1 );
    children[ 0 ] = IntField{ 50, 51, 52 };
    IntField l2gNew;
    IntField nChild;

    RebuildZoneL2gNewAfterInterfaceSplit( parents, children, l2gNew, nChild );

    EXPECT_EQ( l2gNew, ( IntField{ 50, 51, 52 } ) );
    // parent slot holds count 3; one trailing 0 per child entry
    EXPECT_EQ( nChild, ( IntField{ 3, 0, 0, 0 } ) );
}

TEST( RebuildZoneL2gNewTest, MixedSplitAndUnsplit )
{
    const IntField parents{ 1, 2, 3 };
    LinkField children( 3 );
    children[ 0 ] = IntField{};        // keep 1
    children[ 1 ] = IntField{ 20, 21 }; // split 2
    children[ 2 ] = IntField{};        // keep 3
    IntField l2gNew;
    IntField nChild;

    RebuildZoneL2gNewAfterInterfaceSplit( parents, children, l2gNew, nChild );

    EXPECT_EQ( l2gNew, ( IntField{ 1, 20, 21, 3 } ) );
    EXPECT_EQ( nChild, ( IntField{ 0, 2, 0, 0, 0 } ) );
}

TEST( RebuildZoneL2gNewTest, RejectsNegativeParentIdWithoutChangingOutputs )
{
    const IntField parents{ 1, -2 };
    LinkField children( 2 );
    IntField l2gNew{ 90, 91 };
    IntField nChild{ 7, 8 };

    EXPECT_THROW(
        RebuildZoneL2gNewAfterInterfaceSplit( parents, children, l2gNew, nChild ),
        std::out_of_range );

    EXPECT_EQ( l2gNew, ( IntField{ 90, 91 } ) );
    EXPECT_EQ( nChild, ( IntField{ 7, 8 } ) );
}

TEST( RebuildZoneL2gNewTest, RejectsNegativeChildIdWithoutChangingOutputs )
{
    const IntField parents{ 1, 2 };
    LinkField children( 2 );
    children[ 0 ] = IntField{ 10, 11 };
    children[ 1 ] = IntField{ 20, -21 };
    IntField l2gNew{ 90, 91 };
    IntField nChild{ 7, 8 };

    EXPECT_THROW(
        RebuildZoneL2gNewAfterInterfaceSplit( parents, children, l2gNew, nChild ),
        std::out_of_range );

    EXPECT_EQ( l2gNew, ( IntField{ 90, 91 } ) );
    EXPECT_EQ( nChild, ( IntField{ 7, 8 } ) );
}

TEST( RebuildZoneL2gNewTest, RejectsSizeMismatch )
{
    const IntField parents{ 1, 2 };
    LinkField children( 1 );
    IntField l2gNew, nChild;
    EXPECT_THROW(
        RebuildZoneL2gNewAfterInterfaceSplit( parents, children, l2gNew, nChild ),
        std::logic_error );
}
