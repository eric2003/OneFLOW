#include <gtest/gtest.h>
#include <cmath>

#include "TurbSaPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( SaProdStrain, ClassicalZeroWhenNoStrain )
{
    EXPECT_DOUBLE_EQ( SaStrainProductionClassical( 0.0, 0.0, 0.0 ), 0.0 );
}

TEST( SaProdStrain, ClassicalSimple )
{
    // work1=0, work2=0.5, work3=0 => stp = 1, S = 1
    EXPECT_DOUBLE_EQ( SaStrainProductionClassical( 0.0, 0.5, 0.0 ), 1.0 );
}

TEST( SaProdStrain, RotationalV4ReducesWhenOmegaExceedsRoot )
{
    // omega=4, stp=1 => root=1 => 4 - 2*max(0,3) = 4-6 = -2
    EXPECT_DOUBLE_EQ( SaStrainProductionRotationalV4( 4.0, 1.0 ), -2.0 );
}

TEST( SaProdStrain, RotationalV5AtEqualOmegaAndStp )
{
    // omega^2 = stp => max term 0 => result omega
    EXPECT_DOUBLE_EQ( SaStrainProductionRotationalV5( 2.0, 4.0 ), 2.0 );
}
