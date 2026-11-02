// tests/geometry/FillUnsetStagedBcSlotsTest.cpp
// Locks residual BC fill (ScanInterfaceBc) and empty childFid counting.

#include <gtest/gtest.h>
#include <stdexcept>

#include "FaceSolver.h"
#include "HXDefine.h"

using namespace ONEFLOW;

TEST( FillUnsetStagedBcSlotsTest, FillsOnlyInvalidSlots )
{
    IntField types{ INVALID_INDEX, 10, INVALID_INDEX };
    IntField names{ -9, 2, -9 };

    EXPECT_EQ(
        FillUnsetStagedBcSlots( types, names, /*fillBcType=*/0, /*fillNameId=*/-1 ),
        2 );
    EXPECT_EQ( types, ( IntField{ 0, 10, 0 } ) );
    EXPECT_EQ( names, ( IntField{ -1, 2, -1 } ) );
}

TEST( FillUnsetStagedBcSlotsTest, NoOpWhenAllSet )
{
    IntField types{ 1, 2 };
    IntField names{ 3, 4 };
    EXPECT_EQ( FillUnsetStagedBcSlots( types, names, 0, -1 ), 0 );
    EXPECT_EQ( types, ( IntField{ 1, 2 } ) );
}

TEST( FillUnsetStagedBcSlotsTest, RejectsLengthMismatch )
{
    IntField types{ INVALID_INDEX };
    IntField names{ 1, 2 };
    EXPECT_THROW(
        FillUnsetStagedBcSlots( types, names, 0, -1 ),
        std::logic_error );
}

TEST( CountEmptyChildFaceListsTest, CountsEmptyLists )
{
    LinkField childFid;
    childFid.push_back( IntField{} );
    childFid.push_back( IntField{ 1, 2 } );
    childFid.push_back( IntField{} );
    EXPECT_EQ( CountEmptyChildFaceLists( childFid ), 2 );
}

TEST( CountEmptyChildFaceListsTest, EmptyTopLevelIsZero )
{
    EXPECT_EQ( CountEmptyChildFaceLists( LinkField{} ), 0 );
}
