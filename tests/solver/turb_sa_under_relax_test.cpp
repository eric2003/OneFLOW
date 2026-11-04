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

// Intermediate: SaUnderRelaxedUpdate on a relax-coefficient grid.
// Ref: qOld=1, deltaQ=2; q = 1 + 2*alpha.
TEST( SaUnderRelax, SaUnderRelaxedUpdateCurveMatchesClosedForm )
{
    const Real qOld = 1.0;
    const Real deltaQ = 2.0;
    const Real alphas[] = { 0.0, 0.25, 0.5, 0.75, 1.0 };
    for ( Real a : alphas )
    {
        EXPECT_NEAR(
            SaUnderRelaxedUpdate( qOld, deltaQ, a ),
            qOld + a * deltaQ,
            1.0e-12 );
    }
}

TEST( SaUnderRelax, SaUnderRelaxedUpdateIsAffineInDelta )
{
    EXPECT_NEAR(
        SaUnderRelaxedUpdate( 0.0, 4.0, 0.5 ),
        2.0 * SaUnderRelaxedUpdate( 0.0, 2.0, 0.5 ),
        1.0e-12 );
}
