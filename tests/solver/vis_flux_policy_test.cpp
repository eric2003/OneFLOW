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
