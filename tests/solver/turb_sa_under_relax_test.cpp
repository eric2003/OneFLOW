#include <gtest/gtest.h>

#include "TurbSaPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( SaUnderRelax, FullRelaxIsExplicitUpdate )
{
    EXPECT_DOUBLE_EQ( SaUnderRelaxedUpdate( 1.0, 2.0, 1.0 ), 3.0 );
}

TEST( SaUnderRelax, ZeroRelaxKeepsOld )
{
    EXPECT_DOUBLE_EQ( SaUnderRelaxedUpdate( 1.0, 2.0, 0.0 ), 1.0 );
}

TEST( SaUnderRelax, HalfRelax )
{
    // (1-0.5)*1 + 0.5*(1+2) = 0.5 + 1.5 = 2
    EXPECT_DOUBLE_EQ( SaUnderRelaxedUpdate( 1.0, 2.0, 0.5 ), 2.0 );
}
