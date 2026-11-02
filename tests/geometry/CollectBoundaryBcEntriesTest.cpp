// tests/geometry/CollectBoundaryBcEntriesTest.cpp
// Locks boundary BC export order: open right cell faces in face order.

#include <gtest/gtest.h>
#include <stdexcept>

#include "FaceSolver.h"
#include "HXDefine.h"

using namespace ONEFLOW;

TEST( CollectBoundaryBcEntriesTest, CollectsOnlyOpenRightFacesInOrder )
{
    // faces: B, I, B  (right INVALID, 5, INVALID)
    const IntField right{ INVALID_INDEX, 5, INVALID_INDEX };
    const IntField types{ 10, 20, 30 };
    const IntField names{ 1, 2, 3 };
    IntField outT, outN;

    EXPECT_EQ( CollectBoundaryBcEntries( right, types, names, outT, outN ), 2 );
    EXPECT_EQ( outT, ( IntField{ 10, 30 } ) );
    EXPECT_EQ( outN, ( IntField{ 1, 3 } ) );
}

TEST( CollectBoundaryBcEntriesTest, AppendsToExistingOutput )
{
    const IntField right{ INVALID_INDEX };
    const IntField types{ 7 };
    const IntField names{ 8 };
    IntField outT{ 99 };
    IntField outN{ 98 };

    EXPECT_EQ( CollectBoundaryBcEntries( right, types, names, outT, outN ), 1 );
    EXPECT_EQ( outT, ( IntField{ 99, 7 } ) );
    EXPECT_EQ( outN, ( IntField{ 98, 8 } ) );
}

TEST( CollectBoundaryBcEntriesTest, NoBoundaryReturnsZero )
{
    const IntField right{ 1, 2 };
    const IntField types{ 1, 2 };
    const IntField names{ 1, 2 };
    IntField outT, outN;
    EXPECT_EQ( CollectBoundaryBcEntries( right, types, names, outT, outN ), 0 );
    EXPECT_TRUE( outT.empty() );
}

TEST( CollectBoundaryBcEntriesTest, RejectsLengthMismatch )
{
    const IntField right{ INVALID_INDEX, INVALID_INDEX };
    const IntField types{ 1 };
    const IntField names{ 1, 2 };
    IntField outT, outN;
    EXPECT_THROW(
        CollectBoundaryBcEntries( right, types, names, outT, outN ),
        std::logic_error );
}
