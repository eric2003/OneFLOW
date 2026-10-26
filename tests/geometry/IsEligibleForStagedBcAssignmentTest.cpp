// tests/geometry/IsEligibleForStagedBcAssignmentTest.cpp
// Locks which staged faces may receive BOCO BC assignment.

#include <gtest/gtest.h>

#include "FaceSolver.h"
#include "HXDefine.h"

using namespace ONEFLOW;

TEST( IsEligibleForStagedBcAssignmentTest, UnsetBoundaryEligibleInBoundaryOnlyMode )
{
    EXPECT_TRUE(
        IsEligibleForStagedBcAssignment(
            INVALID_INDEX, INVALID_INDEX, /*boundaryFacesOnly=*/true ) );
}

TEST( IsEligibleForStagedBcAssignmentTest, InteriorSkippedInBoundaryOnlyMode )
{
    EXPECT_FALSE(
        IsEligibleForStagedBcAssignment(
            INVALID_INDEX, /*rightCell=*/5, /*boundaryFacesOnly=*/true ) );
}

TEST( IsEligibleForStagedBcAssignmentTest, InteriorEligibleInDetailMode )
{
    EXPECT_TRUE(
        IsEligibleForStagedBcAssignment(
            INVALID_INDEX, /*rightCell=*/5, /*boundaryFacesOnly=*/false ) );
}

TEST( IsEligibleForStagedBcAssignmentTest, AlreadySetNeverEligible )
{
    EXPECT_FALSE(
        IsEligibleForStagedBcAssignment(
            10, INVALID_INDEX, true ) );
    EXPECT_FALSE(
        IsEligibleForStagedBcAssignment(
            10, 5, false ) );
}
