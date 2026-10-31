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

// Source-term components are checked against independently evaluated
// closed-form reference values from the equations documented in sa_source.rst.
TEST( TurbSaPolicy, SaProductionMatchesClosedForm )
{
    EXPECT_NEAR( SaProduction( 0.1355, 0.2, 3.0, 0.4 ), 0.13008, 1.0e-14 );
}

TEST( TurbSaPolicy, SaDestructionMatchesClosedForm )
{
    EXPECT_NEAR( SaDestruction( 3.0, 0.5, 0.1355, 0.2, 0.4 ), 0.58916, 1.0e-14 );
}

TEST( TurbSaPolicy, SaDiffusionMatchesClosedForm )
{
    EXPECT_NEAR( SaDiffusion( 1.3, 0.25, 0.5 ), 0.1625, 1.0e-14 );
}

TEST( TurbSaPolicy, SaLimitedDiffusionSelectsTheSmallerContribution )
{
    EXPECT_DOUBLE_EQ( SaLimitedDiffusion( 0.4, 0.5, 0.3 ), 0.15 );
    EXPECT_DOUBLE_EQ( SaLimitedDiffusion( 0.1, 0.5, 0.3 ), 0.1 );
}

TEST( TurbSaPolicy, SaSourceResidualCombinesProductionDestructionAndDiffusion )
{
    EXPECT_NEAR( SaSourceResidual( 0.13008, 0.58916, 0.1625 ), -0.29658, 1.0e-14 );
}

TEST( TurbSaPolicy, SaNegativePartMatchesPiecewiseDefinition )
{
    EXPECT_DOUBLE_EQ( SaNegativePart( 2.0 ), 0.0 );
    EXPECT_DOUBLE_EQ( SaNegativePart( 0.0 ), 0.0 );
    EXPECT_DOUBLE_EQ( SaNegativePart( -2.0 ), 2.0 );
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
