#include <gtest/gtest.h>

#include "TimeStepPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( TimeStepPolicy, InviscidFaceSpectralRadiusMatchesFormula )
{
    EXPECT_DOUBLE_EQ(
        InviscidFaceSpectralRadius( -3.0, 1.0, 4.0 ),
        half * 4.0 * ( 3.0 + 1.0 ) );
}

TEST( TimeStepPolicy, ViscousFaceSpectralRadiusModel3IsArea2OverReRho )
{
    EXPECT_DOUBLE_EQ(
        ViscousFaceSpectralRadiusModel3( 2.0, 100.0, 3.0 ),
        9.0 / ( 100.0 * 2.0 ) );
}

TEST( TimeStepPolicy, CellTimeStepFromSpectrumIsCflVolOverLambda )
{
    EXPECT_DOUBLE_EQ( CellTimeStepFromSpectrum( 0.5, 4.0, 2.0 ), 1.0 );
}

TEST( TimeStepPolicy, CombineInvAndVisTimeStepIsHarmonic )
{
    EXPECT_DOUBLE_EQ( CombineInvAndVisTimeStep( 2.0, 2.0 ), 1.0 );
}
