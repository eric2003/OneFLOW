#include <gtest/gtest.h>
#include "NsVisFluxPolicy.h"
#include "HXMath.h"
using namespace ONEFLOW;

TEST( VisFluxPolicy, StokesStressTensorShearOnly )
{
    Real txx, tyy, tzz, txy, txz, tyz;
    StokesStressTensor( 2.0, 0,1,0, 0,0,0, 0,0,0, txx,tyy,tzz,txy,txz,tyz );
    EXPECT_DOUBLE_EQ( txy, 2.0 );
}

TEST( VisFluxPolicy, ProjectStressMomentumFluxDotNormal )
{
    Real fx, fy, fz;
    ProjectStressMomentumFlux( 1,0,0, 2,0,0, 3,4,0, fx,fy,fz );
    EXPECT_DOUBLE_EQ( fx, 2.0 );
    EXPECT_DOUBLE_EQ( fy, 3.0 );
    EXPECT_DOUBLE_EQ( fz, 4.0 );
}

TEST( VisFluxPolicy, ViscousEnergyFluxAddsHeat )
{
    EXPECT_DOUBLE_EQ( ViscousEnergyFlux( 1,0,0, 2,0,0, 5 ), 7.0 );
}

TEST( VisFluxPolicy, EffectiveThermalConductivityCombinesLaminarAndTurbulent )
{
    EXPECT_DOUBLE_EQ( EffectiveThermalConductivity( 2.0, 4.0, 0.5, 0.25, 1.0 ), 2.0 );
}

TEST( VisFluxPolicy, ConductiveHeatFluxNormalIsProduct )
{
    EXPECT_DOUBLE_EQ( ConductiveHeatFluxNormal( 3.0, -2.0 ), -6.0 );
}

TEST( VisFluxPolicy, FaceNormalDotIsLinear )
{
    EXPECT_DOUBLE_EQ( FaceNormalDot( 0.0, 1.0, 0.0, 1.0, 5.0, 1.0 ), 5.0 );
}

TEST( VisFluxPolicy, ViscousFaceFluxCoefficientIsMinusReInvTimesArea )
{
    EXPECT_DOUBLE_EQ( ViscousFaceFluxCoefficient( 0.5, 4.0 ), -2.0 );
}

TEST( VisFluxPolicy, DepositFaceScalarToCellsAddsBothSides )
{
    Real L = 1.0, R = 2.0;
    DepositFaceScalarToCells( 3.0, L, R );
    EXPECT_DOUBLE_EQ( L, 4.0 );
    EXPECT_DOUBLE_EQ( R, 5.0 );
}

// Intermediate: EffectiveThermalConductivity on a mut grid.
// Ref: mu=1, invPr_l=1/0.72, invPr_t=1/0.9, cp=1004.5.
TEST( VisFluxPolicy, EffectiveThermalConductivityCurveMatchesClosedForm )
{
    const Real mu = 1.0;
    const Real invPrL = 1.0 / 0.72;
    const Real invPrT = 1.0 / 0.9;
    const Real cp = 1004.5;
    const Real muts[] = { 0.0, 0.5, 1.0, 2.0, 5.0, 10.0 };
    for ( Real mut : muts )
    {
        const Real expected = ( mu * invPrL + mut * invPrT ) * cp;
        EXPECT_NEAR(
            EffectiveThermalConductivity( mu, mut, invPrL, invPrT, cp ),
            expected,
            1.0e-10 );
    }
}

TEST( VisFluxPolicy, EffectiveThermalConductivityLinearInMut )
{
    const Real k0 = EffectiveThermalConductivity( 1.0, 0.0, 1.0, 1.0, 2.0 );
    const Real k1 = EffectiveThermalConductivity( 1.0, 1.0, 1.0, 1.0, 2.0 );
    EXPECT_NEAR( k1 - k0, 2.0, 1.0e-12 ); // + mut*invPrT*cp = 1*1*2
}

// Intermediate: StokesStressTensor pure shear on a viscosity grid.
// Gradients: du/dy=1, all others 0 => txy=mu, normals 0.
TEST( VisFluxPolicy, StokesStressTensorPureShearCurveMatchesClosedForm )
{
    const Real mus[] = { 0.1, 0.5, 1.0, 2.0, 5.0, 10.0 };
    for ( Real mu : mus )
    {
        Real txx, tyy, tzz, txy, txz, tyz;
        StokesStressTensor(
            mu,
            0.0, 1.0, 0.0,
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            txx, tyy, tzz, txy, txz, tyz );
        EXPECT_NEAR( txy, mu, 1.0e-12 );
        EXPECT_NEAR( txx, 0.0, 1.0e-12 );
        EXPECT_NEAR( tyy, 0.0, 1.0e-12 );
        EXPECT_NEAR( tzz, 0.0, 1.0e-12 );
        EXPECT_NEAR( txz, 0.0, 1.0e-12 );
        EXPECT_NEAR( tyz, 0.0, 1.0e-12 );
    }
}

TEST( VisFluxPolicy, StokesStressTensorBulkStrainHasFourThirds )
{
    // Pure extension du/dx=1, others 0:
    // divv = 1, divv*2/3 = 2/3
    // txx = mu*(2 - 2/3) = mu*(4/3)
    Real txx, tyy, tzz, txy, txz, tyz;
    StokesStressTensor(
        3.0,
        1.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        txx, tyy, tzz, txy, txz, tyz );
    EXPECT_NEAR( txx, 3.0 * ( 4.0 / 3.0 ), 1.0e-12 );
    EXPECT_NEAR( tyy, 3.0 * ( 0.0 - two3rd ), 1.0e-12 );
    EXPECT_NEAR( tzz, 3.0 * ( 0.0 - two3rd ), 1.0e-12 );
}
