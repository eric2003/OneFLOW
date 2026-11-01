#include <gtest/gtest.h>

#include "TurbViscosityPolicy.h"

using namespace ONEFLOW;

TEST( TurbViscosityPolicy, CapsEddyViscosityAtLaminarRatio )
{
    EXPECT_DOUBLE_EQ(
        TurbCappedEddyViscosity( 100.0, 2.0, 10.0 ), 20.0 );
}

TEST( TurbViscosityPolicy, PreservesEddyViscosityBelowLaminarRatio )
{
    EXPECT_DOUBLE_EQ(
        TurbCappedEddyViscosity( 5.0, 2.0, 10.0 ), 5.0 );
}
