#include <gtest/gtest.h>
#include "NsVisFluxPolicy.h"
#include "ViscousFluxArithmetic.h"
#include "Sutherland.h"
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

// Intermediate: ProjectStressMomentumFlux on an nx grid.
// Ref stress: txx=2, txy=3, txz=4, others 0; n=(nx,0,0).
TEST( VisFluxPolicy, ProjectStressMomentumFluxCurveMatchesClosedForm )
{
    const Real txx = 2.0, tyy = 0.0, tzz = 0.0;
    const Real txy = 3.0, txz = 4.0, tyz = 0.0;
    const Real nxs[] = { -1.0, -0.5, 0.0, 0.25, 0.5, 1.0 };
    for ( Real nx : nxs )
    {
        Real fx, fy, fz;
        ProjectStressMomentumFlux(
            nx, 0.0, 0.0,
            txx, tyy, tzz, txy, txz, tyz,
            fx, fy, fz );
        EXPECT_NEAR( fx, nx * txx, 1.0e-12 );
        EXPECT_NEAR( fy, nx * txy, 1.0e-12 );
        EXPECT_NEAR( fz, nx * txz, 1.0e-12 );
    }
}

TEST( VisFluxPolicy, ProjectStressMomentumFluxUsesFullNormal )
{
    Real fx, fy, fz;
    ProjectStressMomentumFlux(
        0.0, 1.0, 0.0,
        0.0, 5.0, 0.0, 0.0, 0.0, 7.0,
        fx, fy, fz );
    // n=(0,1,0): fy = tyy = 5, fz = tyz = 7, fx = txy = 0
    EXPECT_NEAR( fx, 0.0, 1.0e-12 );
    EXPECT_NEAR( fy, 5.0, 1.0e-12 );
    EXPECT_NEAR( fz, 7.0, 1.0e-12 );
}

// Intermediate: ConductiveHeatFluxNormal on a temperature-gradient grid.
// Ref: k=2.5; q_n = k * dT_dn.
TEST( VisFluxPolicy, ConductiveHeatFluxNormalCurveMatchesClosedForm )
{
    const Real k = 2.5;
    const Real grads[] = { -4.0, -1.0, 0.0, 0.5, 1.0, 2.0, 4.0 };
    for ( Real dT : grads )
    {
        EXPECT_NEAR( ConductiveHeatFluxNormal( k, dT ), k * dT, 1.0e-12 );
    }
}

TEST( VisFluxPolicy, ConductiveHeatFluxNormalLinearInConductivity )
{
    EXPECT_NEAR(
        ConductiveHeatFluxNormal( 4.0, 1.5 ),
        2.0 * ConductiveHeatFluxNormal( 2.0, 1.5 ),
        1.0e-12 );
}

// Intermediate: ViscousEnergyFlux on a heat-flux grid.
// Ref: u=(1,0,0), f_m=(2,0,0) => work = 2; E = 2 + qn.
TEST( VisFluxPolicy, ViscousEnergyFluxCurveMatchesClosedForm )
{
    const Real qns[] = { -4.0, -1.0, 0.0, 1.0, 2.5, 5.0 };
    for ( Real qn : qns )
    {
        EXPECT_NEAR(
            ViscousEnergyFlux( 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, qn ),
            2.0 + qn,
            1.0e-12 );
    }
}

TEST( VisFluxPolicy, ViscousEnergyFluxIncludesAllWorkComponents )
{
    EXPECT_NEAR(
        ViscousEnergyFlux( 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 ),
        39.0,
        1.0e-12 );
}

// Intermediate: ViscousFaceFluxCoefficient on a face-area grid.
// Ref: invReynolds=0.01; coef = -invRe * A.
TEST( VisFluxPolicy, ViscousFaceFluxCoefficientCurveMatchesClosedForm )
{
    const Real invRe = 0.01;
    const Real areas[] = { 0.5, 1.0, 2.0, 4.0, 8.0 };
    for ( Real A : areas )
    {
        EXPECT_NEAR(
            ViscousFaceFluxCoefficient( invRe, A ),
            -invRe * A,
            1.0e-12 );
    }
}

TEST( VisFluxPolicy, ViscousFaceFluxCoefficientLinearInInvReynolds )
{
    EXPECT_NEAR(
        ViscousFaceFluxCoefficient( 0.2, 3.0 ),
        2.0 * ViscousFaceFluxCoefficient( 0.1, 3.0 ),
        1.0e-12 );
}

// The device port evaluates the same physics from a header shared by every
// execution space. These tests pin the shared arithmetic to the host policy
// function it must reproduce; without them the two transcriptions can drift
// apart silently and only a full cross-backend run would notice.
TEST( VisFluxPolicy, SharedArithmeticMatchesHostStokesStress )
{
    // A non-trivial state: every velocity derivative non-zero, so the
    // divergence correction and all three shear terms are exercised.
    const Real viscosity = 3.7e-4;
    const Real dudx = 0.11, dudy = -0.23, dudz = 0.37;
    const Real dvdx = 0.05, dvdy = 0.19, dvdz = -0.41;
    const Real dwdx = -0.07, dwdy = 0.29, dwdz = 0.13;

    Real txx, tyy, tzz, txy, txz, tyz;
    StokesStressTensor( viscosity, dudx, dudy, dudz, dvdx, dvdy, dvdz,
                        dwdx, dwdy, dwdz, txx, tyy, tzz, txy, txz, tyz );

    Real sxx, syy, szz, sxy, sxz, syz;
    ViscousFluxArithmetic::StokesStress(
        viscosity, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz,
        sxx, syy, szz, sxy, sxz, syz );

    EXPECT_DOUBLE_EQ( sxx, txx );
    EXPECT_DOUBLE_EQ( syy, tyy );
    EXPECT_DOUBLE_EQ( szz, tzz );
    EXPECT_DOUBLE_EQ( sxy, txy );
    EXPECT_DOUBLE_EQ( sxz, txz );
    EXPECT_DOUBLE_EQ( syz, tyz );
}

TEST( VisFluxPolicy, SharedArithmeticMatchesHostProjectionAndEnergy )
{
    const Real xNormal = 0.6, yNormal = -0.48, zNormal = 0.64;
    const Real txx = 1.3, tyy = -0.7, tzz = 0.4;
    const Real txy = 0.22, txz = -0.35, tyz = 0.11;

    Real fx, fy, fz;
    ProjectStressMomentumFlux( xNormal, yNormal, zNormal,
                               txx, tyy, tzz, txy, txz, tyz, fx, fy, fz );
    Real sx, sy, sz;
    ViscousFluxArithmetic::ProjectStressMomentum(
        xNormal, yNormal, zNormal, txx, tyy, tzz, txy, txz, tyz, sx, sy, sz );
    EXPECT_DOUBLE_EQ( sx, fx );
    EXPECT_DOUBLE_EQ( sy, fy );
    EXPECT_DOUBLE_EQ( sz, fz );

    const Real u = 2.1, v = -0.4, w = 0.9, heat = 0.027;
    EXPECT_DOUBLE_EQ(
        ViscousFluxArithmetic::ViscousEnergy( u, v, w, fx, fy, fz, heat ),
        ViscousEnergyFlux( u, v, w, fx, fy, fz, heat ) );
}

TEST( VisFluxPolicy, SharedArithmeticMatchesHostCoefficients )
{
    EXPECT_DOUBLE_EQ(
        ViscousFluxArithmetic::ThermalConductivity( 1.7e-4, 2.3e-5, 1.0 / 0.72,
                                                    1.0 / 0.9, 1004.5 ),
        EffectiveThermalConductivity( 1.7e-4, 2.3e-5, 1.0 / 0.72, 1.0 / 0.9,
                                      1004.5 ) );
    EXPECT_DOUBLE_EQ(
        ViscousFluxArithmetic::ConductiveHeatFlux( 0.031, -2.7 ),
        ConductiveHeatFluxNormal( 0.031, -2.7 ) );
    EXPECT_DOUBLE_EQ(
        ViscousFluxArithmetic::NormalDot( 0.6, -0.48, 0.64, 1.0, 2.0, 3.0 ),
        FaceNormalDot( 0.6, -0.48, 0.64, 1.0, 2.0, 3.0 ) );
    EXPECT_DOUBLE_EQ(
        ViscousFluxArithmetic::FaceFluxCoefficient( 2.5e5, 0.013 ),
        ViscousFaceFluxCoefficient( 2.5e5, 0.013 ) );
}

// The Sutherland form is the classic place to be "mathematically right but
// numerically different": sqrt(t^3) and pow(t, 1.5) disagree in the last bits,
// and the host uses the sqrt form.
TEST( VisFluxPolicy, SharedArithmeticMatchesHostSutherland )
{
    const Real ts = 0.5;
    for ( Real t = 0.6; t <= 1.6; t += 0.2 )
    {
        EXPECT_DOUBLE_EQ(
            ViscousFluxArithmetic::SutherlandViscosityNondim( t, ts ),
            SutherlandViscosityNondim( t, ts ) )
            << "at nondimensional temperature " << t;
    }
    // And it must not equal the pow(t,1.5) transcription, or the test above
    // would be passing for the wrong reason.
    const Real t = 1.1;
    const Real viaPow = std::pow( t, 1.5 ) * ( 1.0 + ts ) / ( t + ts );
    EXPECT_NE( ViscousFluxArithmetic::SutherlandViscosityNondim( t, ts ),
               viaPow );
}
