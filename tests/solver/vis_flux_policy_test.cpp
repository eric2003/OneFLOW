#include <gtest/gtest.h>

#include "NsVisFluxPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( VisFluxPolicy, StokesStressTensorShearOnly )
{
    Real txx, tyy, tzz, txy, txz, tyz;
    // pure shear dudy=1, mu=2
    StokesStressTensor(
        2.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        txx, tyy, tzz, txy, txz, tyz );
    EXPECT_DOUBLE_EQ( txx, 0.0 );
    EXPECT_DOUBLE_EQ( txy, 2.0 );
    EXPECT_DOUBLE_EQ( tyz, 0.0 );
}

TEST( VisFluxPolicy, StokesStressTensorTraceFreeForIsotropicStrain )
{
    Real txx, tyy, tzz, txy, txz, tyz;
    // dudx=dvdy=dwdz=1, mu=3 -> div=3, divv2p3=2
    StokesStressTensor(
        3.0,
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
        txx, tyy, tzz, txy, txz, tyz );
    EXPECT_DOUBLE_EQ( txx, 3.0 * ( 2.0 - 2.0 ) );
    EXPECT_DOUBLE_EQ( tyy, 0.0 );
    EXPECT_DOUBLE_EQ( tzz, 0.0 );
}

TEST( VisFluxPolicy, ProjectStressMomentumFluxDotNormal )
{
    Real fx, fy, fz;
    ProjectStressMomentumFlux(
        1.0, 0.0, 0.0,
        2.0, 0.0, 0.0,
        3.0, 4.0, 0.0,
        fx, fy, fz );
    EXPECT_DOUBLE_EQ( fx, 2.0 );
    EXPECT_DOUBLE_EQ( fy, 3.0 );
    EXPECT_DOUBLE_EQ( fz, 4.0 );
}

TEST( VisFluxPolicy, ViscousEnergyFluxAddsHeat )
{
    EXPECT_DOUBLE_EQ(
        ViscousEnergyFlux( 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 5.0 ), 7.0 );
}
