/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
\*---------------------------------------------------------------------------*/

// Host-runnable parity test for the shared Spalart-Allmaras arithmetic.
//
// TurbSaArithmetic.h is what the HIP kernels compile.  These tests compile the
// same header on the host next to the production host chain and compare the
// two, so the shared code is pinned in the ordinary CPU CTest suite rather than
// only inside a DCU-only binary.  The device-side contract lives in
// tests/euler/device_sa_contract_test.cpp and needs a card.
//
// The host reference is the real thing: TurbCom::CalcSrcSa driven through the
// `turbcom` singleton exactly as UTurbSrcFlux::PrepareCellValue1Equ loads it,
// SaComputeEddyViscosity, and the TurbVisFluxPolicy / TurbInvFluxPolicy
// helpers.  Nothing here re-derives a formula.

#include "TurbSaArithmetic.h"
#include "TurbSaPolicy.h"
#include "TurbViscosityPolicy.h"
#include "TurbVisFluxPolicy.h"
#include "TurbInvFluxPolicy.h"
#include "NsInvFluxPolicy.h"
#include "UGradPolicy.h"
#include "TurbCom.h"
#include "Com.h"

#include <gtest/gtest.h>

#include <cmath>

using namespace ONEFLOW;

namespace
{

constexpr double kTolerance = 1.0e-14;

// TurbCom::InitConst for "1eq-sa" (TurbCom.cpp:174-190) plus the runtime
// inputs the shipped case files supply.  Kept in one place so the host load
// below and the shared-arithmetic inputs cannot drift.
SaModelConstants MakeConstants()
{
    SaModelConstants c;
    const double karm = 0.41;
    const double karm2 = karm * karm;
    const double okarm2 = 1.0 / karm2;
    const double sigma = 2.0 / 3.0;
    const double osigma = 1.0 / sigma;
    const double cv1 = 7.1;

    c.cb1 = 0.1355;
    c.cb2OverSigma = osigma * 0.622;
    c.cv1Cubed = cv1 * cv1 * cv1;
    c.cw2 = 0.3;
    c.cw3ToSix = 64.0;
    c.ct3 = 1.2;
    c.ct4 = 0.5;
    c.karm2 = karm2;
    c.okarm2 = okarm2;
    c.cw1 = c.cb1 * okarm2 + ( 1.0 + 0.622 ) / sigma;
    c.cw1K = c.cw1 * karm2;
    c.sac2 = 0.7;
    c.sac3 = 0.9;
    c.fwStar = 0.424;
    c.productionRatio = 1.0;
    c.invReynolds = 1.0 / 1.0e6;
    c.invSigma = osigma;
    c.maxViscosityRatio = 1.0e6;
    c.ft2Enabled = false;
    c.productionOption = 2;
    return c;
}

void LoadHostConstants( const SaModelConstants & c )
{
    turbcom.cb1 = c.cb1;
    turbcom.cb2 = 0.622;
    turbcom.cb2s = c.cb2OverSigma;
    turbcom.cv1 = 7.1;
    turbcom.cv13 = c.cv1Cubed;
    turbcom.cw1 = c.cw1;
    turbcom.cw2 = c.cw2;
    turbcom.cw3 = 2.0;
    turbcom.cw36 = c.cw3ToSix;
    turbcom.cw1k = c.cw1K;
    turbcom.ct3 = c.ct3;
    turbcom.ct4 = c.ct4;
    turbcom.karm = 0.41;
    turbcom.karm2 = c.karm2;
    turbcom.okarm2 = c.okarm2;
    turbcom.sac2 = c.sac2;
    turbcom.sac3 = c.sac3;
    turbcom.fwStar = c.fwStar;
    turbcom.rprod = c.productionRatio;
    turbcom.oreynolds = c.invReynolds;
    turbcom.sigma = 2.0 / 3.0;
    turbcom.osigma = c.invSigma;
    turbcom.max_vis_ratio = c.maxViscosityRatio;
    turbcom.ft2_flag = c.ft2Enabled ? 1 : 0;
    turbcom.iprod_sa = c.productionOption;
    turbcom.nEqu = 1;
}

struct Cell
{
    double density = 1.0;
    double workingVariable = 0.0;
    double laminarViscosity = 0.0;
    double velocityGradient[ 9 ] = {};
    double workingVariableGradient[ 3 ] = {};
    double lengthScale = 0.0;
    double cellVolume = 1.0;
};

void HostSource( const Cell & cell, const SaModelConstants & c,
                 double & residual, double & implicitSpectrum )
{
    LoadHostConstants( c );
    turbcom.rho = cell.density;
    turbcom.nuet = cell.workingVariable;
    turbcom.visl = cell.laminarViscosity;
    turbcom.dudx = cell.velocityGradient[ 0 ];
    turbcom.dudy = cell.velocityGradient[ 1 ];
    turbcom.dudz = cell.velocityGradient[ 2 ];
    turbcom.dvdx = cell.velocityGradient[ 3 ];
    turbcom.dvdy = cell.velocityGradient[ 4 ];
    turbcom.dvdz = cell.velocityGradient[ 5 ];
    turbcom.dwdx = cell.velocityGradient[ 6 ];
    turbcom.dwdy = cell.velocityGradient[ 7 ];
    turbcom.dwdz = cell.velocityGradient[ 8 ];
    turbcom.dqdxSa = cell.workingVariableGradient[ 0 ];
    turbcom.dqdySa = cell.workingVariableGradient[ 1 ];
    turbcom.dqdzSa = cell.workingVariableGradient[ 2 ];
    turbcom.dist = cell.lengthScale;
    turbcom.len_scale = cell.lengthScale;
    turbcom.cross_term = 0.0;
    gcom.cvol = cell.cellVolume;

    turbcom.CalcSrcSa();
    residual = turbcom.res_sa;
    implicitSpectrum = turbcom.spec_sa;
}

TurbSaArithmetic::SourceResult SharedSource(
    const Cell & cell, const SaModelConstants & c )
{
    TurbSaArithmetic::SourceInput input;
    input.density = cell.density;
    input.workingVariable = cell.workingVariable;
    input.laminarViscosity = cell.laminarViscosity;
    for ( int component = 0; component < 9; ++ component )
    {
        input.velocityGradient[ component ] = cell.velocityGradient[ component ];
    }
    for ( int component = 0; component < 3; ++ component )
    {
        input.workingVariableGradient[ component ] =
            cell.workingVariableGradient[ component ];
    }
    input.lengthScale = cell.lengthScale;
    input.cellVolume = cell.cellVolume;
    input.constants = c;
    return TurbSaArithmetic::ComputeSource( input );
}

// Shear cell: strong rotation, the branch SaModifiedStrainRate takes when
// sBar >= -c2*omega.
Cell ShearCell()
{
    Cell cell;
    cell.density = 1.0;
    cell.workingVariable = 1.0e-4;
    cell.laminarViscosity = 1.0e-5;
    cell.velocityGradient[ 0 ] = 10.0;
    cell.velocityGradient[ 1 ] = 2.0;
    cell.velocityGradient[ 2 ] = -1.0;
    cell.velocityGradient[ 3 ] = 0.5;
    cell.velocityGradient[ 4 ] = -3.0;
    cell.velocityGradient[ 5 ] = 0.25;
    cell.velocityGradient[ 6 ] = 1.5;
    cell.velocityGradient[ 7 ] = 0.75;
    cell.velocityGradient[ 8 ] = 0.1;
    cell.workingVariableGradient[ 0 ] = 0.01;
    cell.workingVariableGradient[ 1 ] = -0.2;
    cell.workingVariableGradient[ 2 ] = 0.03;
    cell.lengthScale = 0.5;
    cell.cellVolume = 2.0e-3;
    return cell;
}

// Near-wall destruction-dominated cell: nonzero implicit spectrum and a
// negative residual, so both outputs carry information.
Cell NearWallCell()
{
    Cell cell;
    cell.density = 1.0;
    cell.workingVariable = 1.0e-2;
    cell.laminarViscosity = 1.0e-5;
    cell.velocityGradient[ 0 ] = 1.0;
    cell.lengthScale = 1.0e-3;
    cell.cellVolume = 1.0e-9;
    return cell;
}

}

TEST( SharedSaArithmetic, ShearCellMatchesHostCalcSrcSa )
{
    const SaModelConstants c = MakeConstants();
    const Cell cell = ShearCell();

    double hostResidual = 0.0;
    double hostSpectrum = 0.0;
    HostSource( cell, c, hostResidual, hostSpectrum );

    const TurbSaArithmetic::SourceResult shared = SharedSource( cell, c );
    EXPECT_NEAR( shared.residual, hostResidual,
                 kTolerance * std::fabs( hostResidual ) );
}

TEST( SharedSaArithmetic, NearWallCellMatchesHostCalcSrcSa )
{
    const SaModelConstants c = MakeConstants();
    const Cell cell = NearWallCell();

    double hostResidual = 0.0;
    double hostSpectrum = 0.0;
    HostSource( cell, c, hostResidual, hostSpectrum );
    ASSERT_LT( hostResidual, 0.0 );
    ASSERT_GT( hostSpectrum, 0.0 );

    const TurbSaArithmetic::SourceResult shared = SharedSource( cell, c );
    EXPECT_NEAR( shared.residual, hostResidual,
                 kTolerance * std::fabs( hostResidual ) );
    EXPECT_NEAR( shared.implicitSpectrum, hostSpectrum,
                 kTolerance * std::fabs( hostSpectrum ) );
}

TEST( SharedSaArithmetic, ProductionOptionTwoIsVorticityMagnitude )
{
    // iprod_sa = 2, the shipped case value, matches no CalcSaProd branch, so
    // CalcSrcSa's `str = work4` leaves the vorticity magnitude.  A shared
    // implementation that assumed the classical branch would disagree here.
    const SaModelConstants c = MakeConstants();
    const Cell cell = ShearCell();
    const TurbSaArithmetic::StrainWork work =
        TurbSaArithmetic::ComputeStrainWork( cell.velocityGradient );
    EXPECT_DOUBLE_EQ(
        TurbSaArithmetic::StrainProduction( work, 2 ), work.work4 );
}

// Discriminating power of the parity assertion above.
//
// A parity test that passes for any implementation proves nothing.  Each case
// here perturbs one input the source term genuinely consumes on a cell where
// that term carries weight, and requires the comparison to reject it:
//
//   * cb1 and cb2OverSigma on the shear cell (production and diffusion
//     dominate there);
//   * cw1K on the near-wall cell (destruction dominates there);
//   * a 1e-10 relative error, below the 1e-14 assertion but not zero.
//
// The 0.0535 relative difference for cb2OverSigma looks small only because the
// diffusion term is a sixth of the residual on that cell; it is still eleven
// orders of magnitude above the tolerance, which is the point.
TEST( SharedSaArithmetic, PerturbationsAreDetected )
{
    const SaModelConstants c = MakeConstants();
    const Cell shear = ShearCell();
    const Cell nearWall = NearWallCell();

    double hostResidual = 0.0;
    double hostSpectrum = 0.0;

    // Baseline: must agree, or the perturbations below cannot be attributed.
    HostSource( shear, c, hostResidual, hostSpectrum );
    const TurbSaArithmetic::SourceResult baseline = SharedSource( shear, c );
    EXPECT_NEAR( baseline.residual, hostResidual,
                 kTolerance * std::fabs( hostResidual ) );

    // Production coefficient.
    {
        SaModelConstants perturbed = c;
        perturbed.cb1 *= 1.01;
        HostSource( shear, c, hostResidual, hostSpectrum );
        const TurbSaArithmetic::SourceResult shared =
            SharedSource( shear, perturbed );
        EXPECT_GT( std::fabs( shared.residual - hostResidual ),
                   kTolerance * std::fabs( hostResidual ) * 1.0e6 )
            << "a wrong cb1 must not pass the parity assertion";
    }

    // Destruction coefficient, on the cell where destruction dominates.
    {
        SaModelConstants perturbed = c;
        perturbed.cw1K *= 1.01;
        HostSource( nearWall, c, hostResidual, hostSpectrum );
        const TurbSaArithmetic::SourceResult shared =
            SharedSource( nearWall, perturbed );
        EXPECT_GT( std::fabs( shared.residual - hostResidual ),
                   kTolerance * std::fabs( hostResidual ) * 1.0e6 )
            << "a wrong cw1K must not pass the parity assertion";
    }

    // Diffusion coefficient.
    {
        SaModelConstants perturbed = c;
        perturbed.cb2OverSigma *= 1.01;
        HostSource( shear, c, hostResidual, hostSpectrum );
        const TurbSaArithmetic::SourceResult shared =
            SharedSource( shear, perturbed );
        EXPECT_GT( std::fabs( shared.residual - hostResidual ),
                   kTolerance * std::fabs( hostResidual ) * 1.0e6 )
            << "a wrong cb2 sigma weight must not pass the parity assertion";
    }

    // A relative error just under the gate must still be rejected by a stricter
    // comparison; this pins that the assertion is not comparing a value with
    // itself.
    {
        HostSource( shear, c, hostResidual, hostSpectrum );
        const double off = hostResidual * ( 1.0 + 1.0e-10 );
        EXPECT_GT( std::fabs( off - hostResidual ),
                   kTolerance * std::fabs( hostResidual ) );
    }
}

TEST( SharedSaArithmetic, EddyViscosityMatchesHostClosure )
{
    const SaModelConstants c = MakeConstants();

    SaEddyViscosityInput input;
    input.density = 1.2;
    input.workingVariable = 1.0e-3;
    input.laminarViscosity = 1.8e-5;
    input.cv1Cubed = c.cv1Cubed;
    input.maxViscosityRatio = c.maxViscosityRatio;
    const SaEddyViscosityResult host = SaComputeEddyViscosity( input );

    TurbSaArithmetic::EddyViscosityInput sharedInput;
    sharedInput.density = input.density;
    sharedInput.workingVariable = input.workingVariable;
    sharedInput.laminarViscosity = input.laminarViscosity;
    sharedInput.constants = c;
    const TurbSaArithmetic::EddyViscosityResult shared =
        TurbSaArithmetic::ComputeEddyViscosity( sharedInput );

    EXPECT_NEAR( shared.eddyViscosity, host.eddyViscosity,
                 kTolerance * std::fabs( host.eddyViscosity ) );
    EXPECT_NEAR( shared.chi, host.chi, kTolerance * std::fabs( host.chi ) );
    EXPECT_NEAR( shared.fv1, host.fv1, kTolerance * std::fabs( host.fv1 ) );
}

TEST( SharedSaArithmetic, InviscidFaceFluxMatchesHostPolicy )
{
    const double normal[ 3 ] = { 0.8, 0.5, 0.2 };
    const double meshVelocityNormal = 0.0;
    const double area = 1.7;
    const double phiLeft = 1.0e-3;
    const double phiRight = 2.0e-3;
    const double velocityLeft[ 3 ] = { 0.30, -0.10, 0.05 };
    const double velocityRight[ 3 ] = { 0.10, -0.05, 0.02 };

    const double vnl = RelativeFaceNormalVelocity(
        normal[ 0 ], normal[ 1 ], normal[ 2 ], meshVelocityNormal,
        velocityLeft[ 0 ], velocityLeft[ 1 ], velocityLeft[ 2 ] );
    const double vnr = RelativeFaceNormalVelocity(
        normal[ 0 ], normal[ 1 ], normal[ 2 ], meshVelocityNormal,
        velocityRight[ 0 ], velocityRight[ 1 ], velocityRight[ 2 ] );

    double vnlPositive = 0.0;
    double vnlNegative = 0.0;
    SplitRelativeNormalVelocity( vnl, vnlPositive, vnlNegative );
    double vnrPositive = 0.0;
    double vnrNegative = 0.0;
    SplitRelativeNormalVelocity( vnr, vnrPositive, vnrNegative );

    // SA sets rho_coef = 1, so the density weight is exactly one.
    const double hostFlux = area * UpwindScalarConvectiveFlux(
        vnlPositive, TurbConvectiveDensityWeight( 1.0, 1.0 ), phiLeft,
        vnrNegative, TurbConvectiveDensityWeight( 1.0, 0.9 ), phiRight );

    TurbSaArithmetic::InviscidFaceInput input;
    input.workingVariableLeft = phiLeft;
    input.workingVariableRight = phiRight;
    input.densityLeft = 1.0;
    input.densityRight = 0.9;
    for ( int component = 0; component < 3; ++ component )
    {
        input.velocityLeft[ component ] = velocityLeft[ component ];
        input.velocityRight[ component ] = velocityRight[ component ];
        input.normal[ component ] = normal[ component ];
    }
    input.meshVelocityNormal = meshVelocityNormal;
    input.faceArea = area;
    input.densityWeight = 1.0;

    EXPECT_NEAR( TurbSaArithmetic::ComputeInviscidFaceFlux( input ), hostFlux,
                 kTolerance * std::fabs( hostFlux ) );
}

TEST( SharedSaArithmetic, DiffusiveDiffusivityMatchesHostPolicy )
{
    const double vislLeft = 1.8e-5;
    const double vislRight = 2.0e-5;
    const double rhoLeft = 1.0;
    const double rhoRight = 0.9;
    const double phiLeft = 1.0e-3;
    const double phiRight = 1.2e-3;

    // CalcFaceVisFlux1Equ:339-343.
    const double hostDiffusivity = SaFaceEffectiveDiffusivity(
        0.5 * ( vislLeft * ( 1.0 / ( rhoLeft + 1.0e-40 ) )
              + vislRight * ( 1.0 / ( rhoRight + 1.0e-40 ) ) ),
        0.5 * ( phiLeft + phiRight ),
        1.5 );

    TurbSaArithmetic::DiffusiveFaceInput input;
    input.workingVariableLeft = phiLeft;
    input.workingVariableRight = phiRight;
    input.laminarViscosityLeft = vislLeft;
    input.laminarViscosityRight = vislRight;
    input.densityLeft = rhoLeft;
    input.densityRight = rhoRight;
    input.invSigma = 1.5;

    EXPECT_NEAR( TurbSaArithmetic::ComputeEffectiveDiffusivity( input ),
                 hostDiffusivity, kTolerance * std::fabs( hostDiffusivity ) );
}

TEST( SharedSaArithmetic, DiffusiveFluxMatchesHostPolicy )
{
    const double normalGradient = -3.5;
    const double area = 0.7;
    const double invReynolds = 1.0e-6;
    const double diffusivity = 0.01;

    // TurbScalarDiffusiveFlux, TurbVisFluxPolicy.cpp:25.
    const double hostFlux = TurbScalarDiffusiveFlux(
        invReynolds, diffusivity, normalGradient, area );

    TurbSaArithmetic::DiffusiveFaceInput input;
    input.densityLeft = 1.0;
    input.densityRight = 1.0;
    input.laminarViscosityLeft = 0.0;
    input.laminarViscosityRight = 0.0;
    input.workingVariableLeft = diffusivity / 1.5;
    input.workingVariableRight = diffusivity / 1.5;
    input.invSigma = 1.5;
    input.normalGradient = normalGradient;
    input.faceArea = area;
    input.invReynolds = invReynolds;

    EXPECT_NEAR( TurbSaArithmetic::ComputeDiffusiveFaceFlux( input ), hostFlux,
                 kTolerance * std::fabs( hostFlux ) );
}

TEST( SharedSaArithmetic, UpdateMatchesHostUnderRelaxation )
{
    const double qOld = 2.0e-4;
    const double delta = -5.0e-5;
    const double relax = 0.7;
    EXPECT_DOUBLE_EQ(
        TurbSaArithmetic::UnderRelaxedUpdate( qOld, delta, relax ),
        SaUnderRelaxedUpdate( qOld, delta, relax ) );
}

// ---------------------------------------------------------------------------
// Gradient: the shared pieces against UGradPolicy.
//
// The full Green-Gauss sum is assembled in a kernel pass, but every numeric
// operation in it comes from TurbSaArithmetic.h, so the pieces are what a unit
// test can pin.  The face-flux test above already covers the face weights
// indirectly through the normal gradient.
// ---------------------------------------------------------------------------
TEST( SharedSaArithmetic, GradientFaceWeightsMatchHostPolicy )
{
    const double dxl = -0.04;
    const double dyl = 0.01;
    const double dzl = 0.005;
    const double dxr = 0.05;
    const double dyr = -0.01;
    const double dzr = -0.004;

    double hostLeft = 0.0;
    double hostRight = 0.0;
    InverseDistanceFaceWeights(
        dxl, dyl, dzl, dxr, dyr, dzr, hostLeft, hostRight );

    double sharedLeft = 0.0;
    double sharedRight = 0.0;
    TurbSaArithmetic::InverseDistanceWeights(
        dxl, dyl, dzl, dxr, dyr, dzr, sharedLeft, sharedRight );

    EXPECT_DOUBLE_EQ( sharedLeft, hostLeft );
    EXPECT_DOUBLE_EQ( sharedRight, hostRight );
}

TEST( SharedSaArithmetic, GradientScaleMatchesHostPolicy )
{
    const double volume = 2.0e-3;
    const double value = 0.0173;
    double hostX = value;
    double hostY = value;
    double hostZ = value;
    ScaleGradientByVolume( volume, hostX, hostY, hostZ );

    double sharedX = value;
    double sharedY = value;
    double sharedZ = value;
    TurbSaArithmetic::ScaleGradientByVolume(
        volume, sharedX, sharedY, sharedZ );

    EXPECT_DOUBLE_EQ( sharedX, hostX );
    EXPECT_DOUBLE_EQ( sharedY, hostY );
    EXPECT_DOUBLE_EQ( sharedZ, hostZ );
}
