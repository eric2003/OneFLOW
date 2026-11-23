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

// Single-face numeric contract between the host Roe flux and the device kernel.
//
// The two implementations must agree to round-off on the same left/right state.
// The host side is evaluated with the shared NsInvFluxPolicy helpers, i.e. the
// functions the solver itself uses, so a change to the device kernel cannot
// drift away from the production formula unnoticed.
//
// Scope: one face, inviscid five-equation Roe with the Harten entropy fix.
// This is a component contract, not a substitute for the 3-step trace or the
// 200-step cross-backend gate.

#include "AccelRuntime.h"
#include "HipFluxBackend.h"
#include "NsInvFluxPolicy.h"

#include <gtest/gtest.h>

#include <cmath>

using namespace ONEFLOW;

namespace
{

struct RefState
{
    double left[ 5 ];    // rho, u, v, w, p
    double right[ 5 ];   // rho, u, v, w, p
    double normal[ 3 ];
    double meshVelocityNormal;
    double entropy1;
    double entropy2;
    bool entropyFix;
    const char * name;
};

// Host Roe, mirroring NsInvFlux::Roe but fed from one face's state.  Uses the
// shared policy helpers so the reference is the production formula.
void HostRoeFlux( const RefState & state, double gammaGas, double flux[ 5 ] )
{
    const double xfn = state.normal[ 0 ];
    const double yfn = state.normal[ 1 ];
    const double zfn = state.normal[ 2 ];
    const double vfn = state.meshVelocityNormal;

    const double rl = state.left[ 0 ];
    const double ul = state.left[ 1 ];
    const double vl = state.left[ 2 ];
    const double wl = state.left[ 3 ];
    const double pl = state.left[ 4 ];
    const double rr = state.right[ 0 ];
    const double ur = state.right[ 1 ];
    const double vr = state.right[ 2 ];
    const double wr = state.right[ 3 ];
    const double pr = state.right[ 4 ];

    const double v2l = ul * ul + vl * vl + wl * wl;
    const double v2r = ur * ur + vr * vr + wr * wr;
    const double hint_l = gammaGas / ( gammaGas - 1.0 ) * pl / rl;
    const double hint_r = gammaGas / ( gammaGas - 1.0 ) * pr / rr;
    const double hl = TotalEnthalpy( hint_l, v2l );
    const double hr = TotalEnthalpy( hint_r, v2r );

    const double vnl = RelativeFaceNormalVelocity(
        xfn, yfn, zfn, vfn, ul, vl, wl );
    const double vnr = RelativeFaceNormalVelocity(
        xfn, yfn, zfn, vfn, ur, vr, wr );

    double rm = 0.0;
    double um = 0.0;
    double vm = 0.0;
    double wm = 0.0;
    double pm = 0.0;
    double hm = 0.0;
    double gamaMid = 0.0;
    RoeAverageState(
        rl, ul, vl, wl, pl, hl, gammaGas,
        rr, ur, vr, wr, pr, hr, gammaGas,
        rm, um, vm, wm, pm, hm, gamaMid );

    const double v2 = um * um + vm * vm + wm * wm;
    const double vnflow = xfn * um + yfn * vm + zfn * wm;
    const double vnrel = vnflow - vfn;
    const double gamm1 = gamaMid - 1.0;
    const double c2 = gamm1 * ( hm - 0.5 * v2 );
    const double cm = std::sqrt( std::fabs( c2 ) );

    double aeig1 = 0.0;
    double aeig2 = 0.0;
    double aeig3 = 0.0;
    RoeAbsoluteEigenvalues( vnrel, cm, aeig1, aeig2, aeig3 );
    double meig1 = 0.0;
    double meig2 = 0.0;
    double meig3 = 0.0;
    ApplyHartenEntropyFix(
        aeig1, aeig2, aeig3,
        state.entropy1, state.entropy2, state.entropyFix,
        meig1, meig2, meig3 );

    double flux1[ 5 ] = {};
    double flux2[ 5 ] = {};
    EulerConvectiveFlux(
        rl, ul, vl, wl, pl, hl, vnl, xfn, yfn, zfn, vfn,
        flux1[ 0 ], flux1[ 1 ], flux1[ 2 ], flux1[ 3 ], flux1[ 4 ] );
    EulerConvectiveFlux(
        rr, ur, vr, wr, pr, hr, vnr, xfn, yfn, zfn, vfn,
        flux2[ 0 ], flux2[ 1 ], flux2[ 2 ], flux2[ 3 ], flux2[ 4 ] );

    double total[ 5 ] = {};
    for ( int eq = 0; eq < 5; ++ eq ) total[ eq ] = flux1[ eq ] + flux2[ eq ];

    // Conserved states, in the same order the solver uses (PrimToQ).
    double q1[ 5 ] = {};
    double q2[ 5 ] = {};
    q1[ 0 ] = rl;
    q1[ 1 ] = rl * ul;
    q1[ 2 ] = rl * vl;
    q1[ 3 ] = rl * wl;
    q1[ 4 ] = pl / ( gammaGas - 1.0 ) + 0.5 * rl * v2l;
    q2[ 0 ] = rr;
    q2[ 1 ] = rr * ur;
    q2[ 2 ] = rr * vr;
    q2[ 3 ] = rr * wr;
    q2[ 4 ] = pr / ( gammaGas - 1.0 ) + 0.5 * rr * v2r;

    double dq[ 5 ] = {};
    for ( int eq = 0; eq < 5; ++ eq ) dq[ eq ] = q2[ eq ] - q1[ eq ];

    double xi1 = 0.0;
    double xi2 = 0.0;
    RoeDissipationWaveCoefficients( meig1, meig2, meig3, c2, cm, xi1, xi2 );

    const double dc = RoeContinuityNormalJump(
        vnflow, xfn, yfn, zfn, dq[ 0 ], dq[ 1 ], dq[ 2 ], dq[ 3 ] );
    const double dh = RoeTotalEnthalpyJump(
        um, vm, wm, gamaMid, dq[ 0 ], dq[ 1 ], dq[ 2 ], dq[ 3 ], dq[ 4 ] );

    double term1 = 0.0;
    double term2 = 0.0;
    RoeDissipationCouplingTerms( dh, dc, c2, xi1, xi2, term1, term2 );

    ApplyRoeDissipationFiveEq(
        meig1, term1, term2, um, vm, wm, hm, vnflow,
        xfn, yfn, zfn,
        dq[ 0 ], dq[ 1 ], dq[ 2 ], dq[ 3 ], dq[ 4 ],
        total[ 0 ], total[ 1 ], total[ 2 ], total[ 3 ], total[ 4 ] );

    for ( int eq = 0; eq < 5; ++ eq ) flux[ eq ] = 0.5 * total[ eq ];
}

RefState MakeState( const char * name )
{
    RefState state = {};
    state.name = name;
    state.normal[ 0 ] = 0.8;
    state.normal[ 1 ] = 0.5;
    state.normal[ 2 ] = 0.2;
    state.meshVelocityNormal = 0.0;
    state.entropy1 = 0.05;
    state.entropy2 = 0.05;
    state.entropyFix = true;
    return state;
}

// `gamma` is also a libm symbol; keep the gas constant under a distinct name.
constexpr double kGamma = 1.4;
constexpr double tolerance = 1.0e-14;

// The backend needs an initialized accelerator runtime; the probe and the
// smoke test do the same before touching device state.
class DeviceRoeContract : public ::testing::Test
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

}

TEST_F( DeviceRoeContract, SubsonicFaceMatchesHost )
{
    RefState state = MakeState( "subsonic" );
    const double left[ 5 ] = { 1.0, 0.30, -0.10, 0.05, 1.0 };
    const double right[ 5 ] = { 0.98, 0.33, -0.08, 0.04, 0.97 };
    for ( int i = 0; i < 5; ++ i )
    {
        state.left[ i ] = left[ i ];
        state.right[ i ] = right[ i ];
    }

    HipFluxBackend::DeviceRoeSample sample;
    for ( int i = 0; i < 5; ++ i )
    {
        sample.left[ i ] = left[ i ];
        sample.right[ i ] = right[ i ];
    }
    for ( int i = 0; i < 3; ++ i ) sample.normal[ i ] = state.normal[ i ];
    sample.gamma = kGamma;
    sample.entropy1 = state.entropy1;
    sample.entropy2 = state.entropy2;
    sample.entropyFix = state.entropyFix;

    double deviceFlux[ 5 ] = {};
    HipFluxBackend::Shared().EvaluateDeviceRoe( sample, deviceFlux );

    double hostFlux[ 5 ] = {};
    HostRoeFlux( state, kGamma, hostFlux );

    for ( int eq = 0; eq < 5; ++ eq )
    {
        EXPECT_NEAR( deviceFlux[ eq ], hostFlux[ eq ], tolerance )
            << "equation " << eq;
    }
}

TEST_F( DeviceRoeContract, TransonicFaceMatchesHost )
{
    RefState state = MakeState( "transonic" );
    const double left[ 5 ] = { 1.0, 0.75, 0.20, -0.10, 1.0 };
    const double right[ 5 ] = { 0.72, 0.90, 0.15, -0.05, 0.62 };
    for ( int i = 0; i < 5; ++ i )
    {
        state.left[ i ] = left[ i ];
        state.right[ i ] = right[ i ];
    }

    HipFluxBackend::DeviceRoeSample sample;
    for ( int i = 0; i < 5; ++ i )
    {
        sample.left[ i ] = left[ i ];
        sample.right[ i ] = right[ i ];
    }
    for ( int i = 0; i < 3; ++ i ) sample.normal[ i ] = state.normal[ i ];
    sample.gamma = kGamma;
    sample.entropy1 = state.entropy1;
    sample.entropy2 = state.entropy2;
    sample.entropyFix = state.entropyFix;

    double deviceFlux[ 5 ] = {};
    HipFluxBackend::Shared().EvaluateDeviceRoe( sample, deviceFlux );

    double hostFlux[ 5 ] = {};
    HostRoeFlux( state, kGamma, hostFlux );

    for ( int eq = 0; eq < 5; ++ eq )
    {
        EXPECT_NEAR( deviceFlux[ eq ], hostFlux[ eq ], tolerance )
            << "equation " << eq;
    }
}

TEST_F( DeviceRoeContract, EntropyFixOffMatchesHost )
{
    RefState state = MakeState( "entropy-fix-off" );
    const double left[ 5 ] = { 1.0, 0.30, -0.10, 0.05, 1.0 };
    const double right[ 5 ] = { 0.98, 0.33, -0.08, 0.04, 0.97 };
    state.entropyFix = false;
    for ( int i = 0; i < 5; ++ i )
    {
        state.left[ i ] = left[ i ];
        state.right[ i ] = right[ i ];
    }

    HipFluxBackend::DeviceRoeSample sample;
    for ( int i = 0; i < 5; ++ i )
    {
        sample.left[ i ] = left[ i ];
        sample.right[ i ] = right[ i ];
    }
    for ( int i = 0; i < 3; ++ i ) sample.normal[ i ] = state.normal[ i ];
    sample.gamma = kGamma;
    sample.entropy1 = state.entropy1;
    sample.entropy2 = state.entropy2;
    sample.entropyFix = state.entropyFix;

    double deviceFlux[ 5 ] = {};
    HipFluxBackend::Shared().EvaluateDeviceRoe( sample, deviceFlux );

    double hostFlux[ 5 ] = {};
    HostRoeFlux( state, kGamma, hostFlux );

    for ( int eq = 0; eq < 5; ++ eq )
    {
        EXPECT_NEAR( deviceFlux[ eq ], hostFlux[ eq ], tolerance )
            << "equation " << eq;
    }
}
