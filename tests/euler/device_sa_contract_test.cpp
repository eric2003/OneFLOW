/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

// Single-cell / single-face numeric contract for the Spalart-Allmaras device
// kernels.
//
// The host side is evaluated with the functions the solver itself calls, not
// with a re-derivation:
//
//   SaComputeEddyViscosity / TurbCappedEddyViscosity   TurbSaPolicy.cpp
//   TurbCom::CalcSrcSa                                 TurbCom.cpp:747
//   TurbInv / UTurbInvFlux's policy helpers            TurbInvFluxPolicy.cpp
//   SaFaceEffectiveDiffusivity / TurbScalarDiffusiveFlux
//                                                      TurbVisFluxPolicy.cpp
//   VisGrad / VisGradGeom                              VisGrad.cpp
//
// TurbCom::CalcSrcSa is a member function driven by the `turbcom` singleton, so
// the source test builds exactly the state UTurbSrcFlux::PrepareCellValue1Equ
// would have loaded and calls it.  That is what makes the test catch a
// transcription error rather than compare two copies of the same mistake.
//
// Scope: one cell for the source, one face for each face operator.  It is a
// component contract, not a substitute for a stage trace or a cross-backend
// residual gate.

#include "AccelRuntime.h"
#include "HipTurbulenceBackend.h"
#include "NsInvFluxPolicy.h"
#include "TurbSaArithmetic.h"
#include "TurbSaPolicy.h"
#include "TurbViscosityPolicy.h"
#include "TurbVisFluxPolicy.h"
#include "TurbInvFluxPolicy.h"
#include "TurbCom.h"

#include <gtest/gtest.h>

#include <cmath>

using namespace ONEFLOW;

namespace
{

// The host writes qf1/qf2 and reads densities and velocities from the NS
// limiter.  Two constants keep the host state and the device sample provably
// identical: they are assigned in one place and copied to both sides.
constexpr double kTolerance = 1.0e-14;

struct SaConstantsFixture
{
    SaModelConstants constants;

    // Reproduces TurbCom::InitConst for "1eq-sa" (TurbCom.cpp:174-190) and
    // then loads the runtime inputs the case file supplies.
    static SaConstantsFixture Make()
    {
        SaConstantsFixture fixture;
        const double karm = 0.41;
        const double karm2 = karm * karm;
        const double okarm2 = 1.0 / karm2;
        const double sigma = 2.0 / 3.0;
        const double osigma = 1.0 / sigma;
        const double cv1 = 7.1;

        SaModelConstants & c = fixture.constants;
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
        return fixture;
    }

    // Push the same values into the host singleton the way TurbCom::Init would.
    void LoadHostCom() const
    {
        const SaModelConstants & c = constants;
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
};

// One cell's worth of the host inputs for the source term.
struct HostCell
{
    double density = 0.0;
    double workingVariable = 0.0;
    double laminarViscosity = 0.0;
    double velocityGradient[ 9 ] = {};
    double workingVariableGradient[ 3 ] = {};
    double wallDistance = 0.0;
    double cellVolume = 0.0;
};

// Drive the host source path exactly as UTurbSrcFlux::PrepareCellValue1Equ
// does, then call the production TurbCom::CalcSrcSa.
void HostSource( const HostCell & cell, const SaConstantsFixture & fixture,
                 double & residual, double & implicitSpectrum )
{
    fixture.LoadHostCom();
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
    turbcom.dist = cell.wallDistance;
    turbcom.len_scale = cell.wallDistance;
    // CalcSrcSa reads `cross_term` only through a two-equation path, but the
    // member must be a defined value because PrepareCellValue1Equ loads it.
    turbcom.cross_term = 0.0;

    turbcom.CalcSrcSa();
    residual = turbcom.res_sa;
    implicitSpectrum = turbcom.spec_sa;
}

class DeviceSaContract : public ::testing::Test
{
protected:
    static void SetUpTestSuite()
    {
        InitializeAccelRuntime( 0, 1 );
    }

    static void TearDownTestSuite()
    {
        AccelRuntime::Instance().Finalize();
    }
};

HostCell MakeCell()
{
    HostCell cell;
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
    cell.wallDistance = 0.5;
    cell.cellVolume = 2.0e-3;
    return cell;
}

// A second state on the other side of SaModifiedStrainRate's branch, so a
// transcription that only handles one branch cannot pass.
HostCell MakeCellSlowStrain()
{
    HostCell cell = MakeCell();
    cell.velocityGradient[ 0 ] = 0.05;
    cell.velocityGradient[ 1 ] = 0.01;
    cell.velocityGradient[ 2 ] = 0.0;
    cell.velocityGradient[ 3 ] = 0.0;
    cell.velocityGradient[ 4 ] = 0.02;
    cell.velocityGradient[ 5 ] = 0.0;
    cell.velocityGradient[ 6 ] = 0.0;
    cell.velocityGradient[ 7 ] = 0.0;
    cell.velocityGradient[ 8 ] = 0.01;
    cell.workingVariable = 5.0e-3;
    cell.wallDistance = 0.02;
    return cell;
}

// Near-wall destruction-dominated state.  The implicit spectrum is genuinely
// nonzero here, which the two states above cannot exercise: for them spec_sa
// is -0.0, so an equality assertion on it would pass for an empty kernel.
HostCell MakeCellNearWall()
{
    HostCell cell = MakeCell();
    cell.velocityGradient[ 0 ] = 1.0;
    cell.velocityGradient[ 1 ] = 0.0;
    cell.velocityGradient[ 2 ] = 0.0;
    cell.velocityGradient[ 3 ] = 0.0;
    cell.velocityGradient[ 4 ] = 0.0;
    cell.velocityGradient[ 5 ] = 0.0;
    cell.velocityGradient[ 6 ] = 0.0;
    cell.velocityGradient[ 7 ] = 0.0;
    cell.velocityGradient[ 8 ] = 0.0;
    cell.workingVariable = 1.0e-2;
    cell.workingVariableGradient[ 0 ] = 0.0;
    cell.workingVariableGradient[ 1 ] = 0.0;
    cell.workingVariableGradient[ 2 ] = 0.0;
    cell.wallDistance = 1.0e-3;
    cell.cellVolume = 1.0e-9;
    return cell;
}

HipTurbulenceBackend::SourceSample MakeSample(
    const HostCell & cell, const SaConstantsFixture & fixture )
{
    HipTurbulenceBackend::SourceSample sample;
    sample.density = cell.density;
    sample.workingVariable = cell.workingVariable;
    sample.laminarViscosity = cell.laminarViscosity;
    for ( int component = 0; component < 9; ++ component )
    {
        sample.velocityGradient[ component ] = cell.velocityGradient[ component ];
    }
    for ( int component = 0; component < 3; ++ component )
    {
        sample.workingVariableGradient[ component ] =
            cell.workingVariableGradient[ component ];
    }
    sample.lengthScale = cell.wallDistance;
    sample.cellVolume = cell.cellVolume;
    sample.constants = fixture.constants;
    return sample;
}

}

// ---------------------------------------------------------------------------
// Source term: one cell.
// ---------------------------------------------------------------------------
TEST_F( DeviceSaContract, SourceMatchesHostOnShearCell )
{
    const SaConstantsFixture fixture = SaConstantsFixture::Make();
    const HostCell cell = MakeCell();

    double hostResidual = 0.0;
    double hostSpectrum = 0.0;
    HostSource( cell, fixture, hostResidual, hostSpectrum );

    double deviceResidual = 0.0;
    double deviceSpectrum = 0.0;
    HipTurbulenceBackend::Shared().EvaluateSourceSample(
        MakeSample( cell, fixture ), deviceResidual, deviceSpectrum );

    EXPECT_NEAR( deviceResidual, hostResidual,
                 kTolerance * std::fabs( hostResidual ) )
        << "res_sa: device " << deviceResidual << " host " << hostResidual;
    EXPECT_NEAR( deviceSpectrum, hostSpectrum,
                 kTolerance * std::max( std::fabs( hostSpectrum ), 1.0e-30 ) )
        << "spec_sa: device " << deviceSpectrum << " host " << hostSpectrum;
}

TEST_F( DeviceSaContract, SourceMatchesHostOnLowStrainCell )
{
    const SaConstantsFixture fixture = SaConstantsFixture::Make();
    const HostCell cell = MakeCellSlowStrain();

    double hostResidual = 0.0;
    double hostSpectrum = 0.0;
    HostSource( cell, fixture, hostResidual, hostSpectrum );

    double deviceResidual = 0.0;
    double deviceSpectrum = 0.0;
    HipTurbulenceBackend::Shared().EvaluateSourceSample(
        MakeSample( cell, fixture ), deviceResidual, deviceSpectrum );

    EXPECT_NEAR( deviceResidual, hostResidual,
                 kTolerance * std::fabs( hostResidual ) );
    EXPECT_NEAR( deviceSpectrum, hostSpectrum,
                 kTolerance * std::max( std::fabs( hostSpectrum ), 1.0e-30 ) );
}

// Near-wall state: nonzero implicit spectrum and a net negative residual, so
// both outputs are discriminating.
TEST_F( DeviceSaContract, SourceMatchesHostNearWall )
{
    const SaConstantsFixture fixture = SaConstantsFixture::Make();
    const HostCell cell = MakeCellNearWall();

    double hostResidual = 0.0;
    double hostSpectrum = 0.0;
    HostSource( cell, fixture, hostResidual, hostSpectrum );
    ASSERT_LT( hostResidual, 0.0 );
    ASSERT_GT( hostSpectrum, 0.0 );

    double deviceResidual = 0.0;
    double deviceSpectrum = 0.0;
    HipTurbulenceBackend::Shared().EvaluateSourceSample(
        MakeSample( cell, fixture ), deviceResidual, deviceSpectrum );

    EXPECT_NEAR( deviceResidual, hostResidual,
                 kTolerance * std::fabs( hostResidual ) );
    EXPECT_NEAR( deviceSpectrum, hostSpectrum,
                 kTolerance * std::fabs( hostSpectrum ) );
}

// The source term must be volume-free: UTurbSrcFlux::Update1Equ multiplies by
// cvol at the accumulation site, so a kernel that also scaled by volume would
// agree here only by accident.
TEST_F( DeviceSaContract, SourceIsPerUnitVolume )
{
    const SaConstantsFixture fixture = SaConstantsFixture::Make();
    HostCell cell = MakeCell();
    const double expected = 0.0;
    ( void ) expected;

    double residualA = 0.0;
    double spectrumA = 0.0;
    HipTurbulenceBackend::Shared().EvaluateSourceSample(
        MakeSample( cell, fixture ), residualA, spectrumA );

    cell.cellVolume *= 100.0;
    double residualB = 0.0;
    double spectrumB = 0.0;
    HipTurbulenceBackend::Shared().EvaluateSourceSample(
        MakeSample( cell, fixture ), residualB, spectrumB );

    EXPECT_DOUBLE_EQ( residualA, residualB );
    EXPECT_DOUBLE_EQ( spectrumA, spectrumB );
}

// ---------------------------------------------------------------------------
// Eddy viscosity: one cell.
// ---------------------------------------------------------------------------
TEST_F( DeviceSaContract, EddyViscosityMatchesHost )
{
    const SaConstantsFixture fixture = SaConstantsFixture::Make();

    SaEddyViscosityInput input;
    input.density = 1.2;
    input.workingVariable = 1.0e-3;
    input.laminarViscosity = 1.8e-5;
    input.cv1Cubed = fixture.constants.cv1Cubed;
    input.maxViscosityRatio = fixture.constants.maxViscosityRatio;
    const SaEddyViscosityResult host = SaComputeEddyViscosity( input );

    HipTurbulenceBackend::EddyViscositySample sample;
    sample.density = input.density;
    sample.workingVariable = input.workingVariable;
    sample.laminarViscosity = input.laminarViscosity;
    sample.constants = fixture.constants;

    double device = 0.0;
    HipTurbulenceBackend::Shared().EvaluateEddyViscositySample(
        sample, device );

    EXPECT_NEAR( device, host.eddyViscosity,
                 kTolerance * std::fabs( host.eddyViscosity ) )
        << "device " << device << " host " << host.eddyViscosity;
}

// The cap is part of the host contract (TurbCappedEddyViscosity).  An
// uncapped device result would differ by orders of magnitude here.
TEST_F( DeviceSaContract, EddyViscosityAppliesRatioCap )
{
    const SaConstantsFixture fixture = SaConstantsFixture::Make();

    // rho * nuet * fv1 must exceed r_max * mu_l by a wide margin, or the cap
    // is not exercised at all.  With rho = 1, mu_l = 1e-5 and r_max = 1e6 the
    // cap is 10, so nuet = 1e3 makes the uncapped value 100x too large.
    SaEddyViscosityInput input;
    input.density = 1.0;
    input.workingVariable = 1.0e3;
    input.laminarViscosity = 1.0e-5;
    input.cv1Cubed = fixture.constants.cv1Cubed;
    input.maxViscosityRatio = 1.0e6;
    const SaEddyViscosityResult host = SaComputeEddyViscosity( input );
    ASSERT_GT( host.uncappedEddyViscosity,
               host.eddyViscosity * 10.0 );
    ASSERT_DOUBLE_EQ( host.eddyViscosity, 1.0e-5 * 1.0e6 );

    HipTurbulenceBackend::EddyViscositySample sample;
    sample.density = input.density;
    sample.workingVariable = input.workingVariable;
    sample.laminarViscosity = input.laminarViscosity;
    sample.constants = fixture.constants;

    double device = 0.0;
    HipTurbulenceBackend::Shared().EvaluateEddyViscositySample(
        sample, device );
    EXPECT_DOUBLE_EQ( device, host.eddyViscosity );
}

// ---------------------------------------------------------------------------
// Scalar inviscid flux: one face.
// ---------------------------------------------------------------------------
TEST_F( DeviceSaContract, InviscidFaceFluxMatchesHost )
{
    // Host side: TurbInv::RoeFlux composed from the policy helpers.
    const double normal[ 3 ] = { 0.8, 0.5, 0.2 };
    const double meshVelocityNormal = 0.0;
    const double area = 1.7;
    const double phiLeft = 1.0e-3;
    const double phiRight = 2.0e-3;
    const double rhoLeft = 1.0;
    const double rhoRight = 0.9;
    const double velocityLeft[ 3 ] = { 0.30, -0.10, 0.05 };
    const double velocityRight[ 3 ] = { 0.10, -0.05, 0.02 };

    double vnl = 0.0;
    double vnr = 0.0;
    vnl = RelativeFaceNormalVelocity(
        normal[ 0 ], normal[ 1 ], normal[ 2 ], meshVelocityNormal,
        velocityLeft[ 0 ], velocityLeft[ 1 ], velocityLeft[ 2 ] );
    vnr = RelativeFaceNormalVelocity(
        normal[ 0 ], normal[ 1 ], normal[ 2 ], meshVelocityNormal,
        velocityRight[ 0 ], velocityRight[ 1 ], velocityRight[ 2 ] );

    double vnlPositive = 0.0;
    double vnlNegative = 0.0;
    SplitRelativeNormalVelocity( vnl, vnlPositive, vnlNegative );
    double vnrPositive = 0.0;
    double vnrNegative = 0.0;
    SplitRelativeNormalVelocity( vnr, vnrPositive, vnrNegative );

    // SA sets rho_coef = 1 (TurbInv::Init:47).
    const double weightLeft = TurbConvectiveDensityWeight( 1.0, rhoLeft );
    const double weightRight = TurbConvectiveDensityWeight( 1.0, rhoRight );
    const double hostFlux = area * UpwindScalarConvectiveFlux(
        vnlPositive, weightLeft, phiLeft,
        vnrNegative, weightRight, phiRight );

    HipTurbulenceBackend::InviscidFaceSample sample;
    sample.workingVariableLeft = phiLeft;
    sample.workingVariableRight = phiRight;
    sample.densityLeft = rhoLeft;
    sample.densityRight = rhoRight;
    for ( int component = 0; component < 3; ++ component )
    {
        sample.velocityLeft[ component ] = velocityLeft[ component ];
        sample.velocityRight[ component ] = velocityRight[ component ];
        sample.normal[ component ] = normal[ component ];
    }
    sample.meshVelocityNormal = meshVelocityNormal;
    sample.faceArea = area;
    sample.densityWeight = 1.0;

    double deviceFlux = 0.0;
    HipTurbulenceBackend::Shared().EvaluateInviscidFaceSample(
        sample, deviceFlux );

    EXPECT_NEAR( deviceFlux, hostFlux,
                 kTolerance * std::fabs( hostFlux ) )
        << "device " << deviceFlux << " host " << hostFlux;
}

// ---------------------------------------------------------------------------
// Diffusive face flux: one face, through the VIS_TEST normal gradient.
// ---------------------------------------------------------------------------
TEST_F( DeviceSaContract, FaceNormalGradientMatchesHostVisTest )
{
    // Geometry: a mildly skewed internal face, so the inverse-distance and
    // skew terms are all exercised.
    const double faceCenter[ 3 ] = { 0.0, 0.0, 0.0 };
    const double cellLeft[ 3 ] = { -0.04, 0.01, 0.005 };
    const double cellRight[ 3 ] = { 0.05, -0.01, -0.004 };
    const double normal[ 3 ] = { 1.0, 0.0, 0.0 };
    const double phiLeft = 1.0e-3;
    const double phiRight = 1.2e-3;
    const double gradientLeft[ 3 ] = { 0.02, 0.01, -0.005 };
    const double gradientRight[ 3 ] = { 0.03, -0.002, 0.004 };

    // Host side: the VisGrad sequence TestMethod runs, on the same values.
    // The helper below is written as the host expression rather than calling
    // VisGrad, because VisGrad is driven entirely by the `vgg`/`gcom`/`ug`
    // globals and a single-face test would have to fabricate a grid for it.
    // Every term is the host formula with its source line named.
    const double dxl = cellLeft[ 0 ] - faceCenter[ 0 ];
    const double dyl = cellLeft[ 1 ] - faceCenter[ 1 ];
    const double dzl = cellLeft[ 2 ] - faceCenter[ 2 ];
    const double dxr = cellRight[ 0 ] - faceCenter[ 0 ];
    const double dyr = cellRight[ 1 ] - faceCenter[ 1 ];
    const double dzr = cellRight[ 2 ] - faceCenter[ 2 ];

    const double d1 = normal[ 0 ] * dxl + normal[ 1 ] * dyl + normal[ 2 ] * dzl;
    const double d2 = normal[ 0 ] * dxr + normal[ 1 ] * dyr + normal[ 2 ] * dzr;

    // VisGrad::AverFaceValue
    const double faceValue = 0.5 * ( phiLeft + phiRight );

    // VisGrad::CalcTestMethod: both signed angles are positive here, so the
    // side values are corrected before the one-sided gradients are formed.
    const double dxnl = normal[ 0 ] * d1 - dxl;
    const double dynl = normal[ 1 ] * d1 - dyl;
    const double dznl = normal[ 2 ] * d1 - dzl;
    const double dxnr = normal[ 0 ] * d2 - dxr;
    const double dynr = normal[ 1 ] * d2 - dyr;
    const double dznr = normal[ 2 ] * d2 - dzr;
    const double q11 = phiLeft + gradientLeft[ 0 ] * dxnl
        + gradientLeft[ 1 ] * dynl + gradientLeft[ 2 ] * dznl;
    const double q22 = phiRight + gradientRight[ 0 ] * dxnr
        + gradientRight[ 1 ] * dynr + gradientRight[ 2 ] * dznr;
    const double dqdn1 = ( q11 - faceValue ) / d1;
    const double dqdn2 = ( q22 - faceValue ) / d2;
    // VisGrad::CalcC1C2
    const double dtmp = d1 * d1 + d2 * d2;
    const double c1 = d1 * d1 / dtmp;
    const double c2 = 1.0 - c1;
    const double hostNormalGradient = c1 * dqdn1 + c2 * dqdn2;

    HipTurbulenceBackend::NormalGradientSample sample;
    sample.workingVariableLeft = phiLeft;
    sample.workingVariableRight = phiRight;
    for ( int component = 0; component < 3; ++ component )
    {
        sample.gradientLeft[ component ] = gradientLeft[ component ];
        sample.gradientRight[ component ] = gradientRight[ component ];
        sample.faceCenter[ component ] = faceCenter[ component ];
        sample.cellLeft[ component ] = cellLeft[ component ];
        sample.cellRight[ component ] = cellRight[ component ];
        sample.normal[ component ] = normal[ component ];
    }
    sample.boundary = false;
    sample.isInterface = false;

    double deviceNormalGradient = 0.0;
    HipTurbulenceBackend::Shared().EvaluateNormalGradientSample(
        sample, deviceNormalGradient );

    EXPECT_NEAR( deviceNormalGradient, hostNormalGradient,
                 kTolerance * std::fabs( hostNormalGradient ) )
        << "device " << deviceNormalGradient
        << " host " << hostNormalGradient;
}

// A boundary face whose bc type is not INTERFACE is forced to the one-sided
// normal gradient by VisGrad::CalcC1C2 (c1 = 1, c2 = 0).  That is the branch a
// solid wall takes, so it must be pinned separately.
TEST_F( DeviceSaContract, FaceNormalGradientUsesOneSidedWeightAtWall )
{
    const double faceCenter[ 3 ] = { 0.0, 0.0, 0.0 };
    const double cellLeft[ 3 ] = { -0.05, 0.0, 0.0 };
    const double cellRight[ 3 ] = { 0.05, 0.0, 0.0 };
    const double normal[ 3 ] = { 1.0, 0.0, 0.0 };
    const double phiLeft = 1.0e-3;
    const double phiRight = 0.0;
    const double gradientLeft[ 3 ] = { 0.02, 0.0, 0.0 };
    const double gradientRight[ 3 ] = { 0.0, 0.0, 0.0 };

    // Host: boundary and not INTERFACE -> c1 = 1, c2 = 0.
    const double d1 = -0.05;
    const double faceValue = 0.5 * ( phiLeft + phiRight );
    const double q11 = phiLeft + gradientLeft[ 0 ] * 0.0;
    const double hostNormalGradient = ( q11 - faceValue ) / d1;

    HipTurbulenceBackend::NormalGradientSample sample;
    sample.workingVariableLeft = phiLeft;
    sample.workingVariableRight = phiRight;
    for ( int component = 0; component < 3; ++ component )
    {
        sample.gradientLeft[ component ] = gradientLeft[ component ];
        sample.gradientRight[ component ] = gradientRight[ component ];
        sample.faceCenter[ component ] = faceCenter[ component ];
        sample.cellLeft[ component ] = cellLeft[ component ];
        sample.cellRight[ component ] = cellRight[ component ];
        sample.normal[ component ] = normal[ component ];
    }
    sample.boundary = true;
    sample.isInterface = false;

    double deviceNormalGradient = 0.0;
    HipTurbulenceBackend::Shared().EvaluateNormalGradientSample(
        sample, deviceNormalGradient );

    EXPECT_NEAR( deviceNormalGradient, hostNormalGradient,
                 kTolerance * std::fabs( hostNormalGradient ) )
        << "device " << deviceNormalGradient
        << " host " << hostNormalGradient;
}

// Diffusive face flux composition: the effective diffusivity is the arithmetic
// mean of the per-side kinematic laminar viscosity plus the mean working
// variable, all times 1/sigma, and the flux carries the -1/Re and the area.
TEST( SaDiffusiveFluxContract, EffectiveDiffusivityMatchesHost )
{
    const double vislLeft = 1.8e-5;
    const double vislRight = 2.0e-5;
    const double rhoLeft = 1.0;
    const double rhoRight = 0.9;
    const double phiLeft = 1.0e-3;
    const double phiRight = 1.2e-3;
    const double invSigma = 1.5;

    // Host: CalcFaceVisFlux1Equ:339-343.
    const double orl = 1.0 / ( rhoLeft + 1.0e-40 );
    const double orr = 1.0 / ( rhoRight + 1.0e-40 );
    const double hostDiffusivity = SaFaceEffectiveDiffusivity(
        0.5 * ( vislLeft * orl + vislRight * orr ),
        0.5 * ( phiLeft + phiRight ),
        invSigma );

    TurbSaArithmetic::DiffusiveFaceInput input;
    input.workingVariableLeft = phiLeft;
    input.workingVariableRight = phiRight;
    input.laminarViscosityLeft = vislLeft;
    input.laminarViscosityRight = vislRight;
    input.densityLeft = rhoLeft;
    input.densityRight = rhoRight;
    input.invSigma = invSigma;

    EXPECT_DOUBLE_EQ(
        TurbSaArithmetic::ComputeEffectiveDiffusivity( input ),
        hostDiffusivity );
}

TEST( SaDiffusiveFluxContract, FluxCarriesMinusInverseReynoldsAndArea )
{
    const double invReynolds = 1.0e-6;
    const double diffusivity = 0.01;
    const double normalGradient = -3.5;
    const double area = 0.7;

    // Host: TurbScalarDiffusiveFlux, TurbVisFluxPolicy.cpp:25.
    const double hostFlux = TurbScalarDiffusiveFlux(
        invReynolds, diffusivity, normalGradient, area );

    TurbSaArithmetic::DiffusiveFaceInput input;
    // Choose the states so the arithmetic reproduces `diffusivity` exactly:
    // equal densities, equal laminar viscosity 0, and both working variables
    // set so their mean is diffusivity / invSigma.
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
                 1.0e-14 * std::fabs( hostFlux ) );
}

// ---------------------------------------------------------------------------
// Update: one cell, host formula SaUnderRelaxedUpdate.
// ---------------------------------------------------------------------------
TEST( SaUpdateContract, UnderRelaxedUpdateMatchesHost )
{
    const double qOld = 2.0e-4;
    const double delta = -5.0e-5;
    const double relax = 0.7;

    // Host: TurbUpdate::CalcFlowField1Equ:64 -> SaUnderRelaxedUpdate.
    const double host = SaUnderRelaxedUpdate( qOld, delta, relax );
    const double device =
        TurbSaArithmetic::UnderRelaxedUpdate( qOld, delta, relax );
    EXPECT_DOUBLE_EQ( device, host );
}
