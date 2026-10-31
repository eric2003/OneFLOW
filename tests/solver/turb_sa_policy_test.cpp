#include <gtest/gtest.h>
#include <cmath>

#include "TurbSaPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( TurbSaPolicy, SaChiIsProduct )
{
    EXPECT_NEAR( SaChi( 2.0, 0.5 ), 1.0 + SMALL, 1.0e-12 );
}

TEST( TurbSaPolicy, SaFv1AtLargeChiApproachesOne )
{
    const Real chi = 100.0;
    const Real cv13 = 7.1 * 7.1 * 7.1;
    EXPECT_NEAR( SaFv1( chi, cv13 ), 1.0, 1.0e-3 );
}

TEST( TurbSaPolicy, SaFv2AtZeroChiIsOne )
{
    EXPECT_DOUBLE_EQ( SaFv2( 0.0, 0.0 ), 1.0 );
}

TEST( TurbSaPolicy, SaModifiedStrainSimpleBranch )
{
    EXPECT_DOUBLE_EQ( SaModifiedStrainRate( 2.0, 1.0, 0.7, 0.9 ), 3.0 );
}

TEST( TurbSaPolicy, SaRRatioIsClamped )
{
    EXPECT_DOUBLE_EQ( SaRRatio( 100.0, 1.0 ), 10.0 );
}

TEST( TurbSaPolicy, SaGAndFwAtRZero )
{
    EXPECT_DOUBLE_EQ( SaG( 0.0, 0.3 ), 0.0 );
}

TEST( TurbSaPolicy, SaFt2DisabledIsZero )
{
    EXPECT_DOUBLE_EQ( SaFt2( 1.0, 1.2, 0.5, false ), 0.0 );
}

TEST( TurbSaPolicy, SaEddyViscosityIsRhoNuFv1 )
{
    EXPECT_DOUBLE_EQ( SaEddyViscosity( 1.2, 0.5, 0.8 ), 0.48 );
}

TEST( TurbSaPolicy, SaCappedEddyViscosityRespectsRatio )
{
    EXPECT_DOUBLE_EQ( SaCappedEddyViscosity( 100.0, 2.0, 10.0 ), 20.0 );
    EXPECT_DOUBLE_EQ( SaCappedEddyViscosity( 5.0, 2.0, 10.0 ), 5.0 );
}

// Intermediate verification: fv1(chi) vs closed-form theory (not a full CFD case).
// cv1 = 7.1 is the standard SA constant (NASA TMR / Spalart-Allmaras).
TEST( TurbSaPolicy, Fv1ChiCurveMatchesClosedForm )
{
    const Real cv1 = 7.1;
    const Real cv13 = cv1 * cv1 * cv1;
    const Real chis[] = { 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0 };
    for ( Real chi : chis )
    {
        const Real chi3 = chi * chi * chi;
        const Real expected = chi3 / ( chi3 + cv13 );
        EXPECT_NEAR( SaFv1( chi, cv13 ), expected, 1.0e-12 );
    }
}
