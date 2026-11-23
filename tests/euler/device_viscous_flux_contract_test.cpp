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

// Single-face numeric contract between the host viscous flux chain and the
// device kernel.
//
// The host side below is the production chain, not a paraphrase of it:
//   UNsVisFlux::PrepareFaceValue  VIS_STD branch (StdMethod)
//   UNsVisFlux::CalcFaceVisFlux   CalcHeatFlux -> CalcStress -> CalcNsVisFlux
//   UNsVisFlux::UpdateFaceVisFlux -A/Re scaling
// every arithmetic step of which delegates to the shared policy helpers
// NsVisFluxPolicy (StokesStressTensor, ProjectStressMomentumFlux,
// ViscousEnergyFlux, EffectiveThermalConductivity,
// ConductiveHeatFluxNormal, ViscousFaceFluxCoefficient) and UGradPolicy
// (InverseDistanceFaceWeights, CorrectFaceGradient).  Those are the functions
// the solver itself calls, so a change to either side cannot drift unnoticed.
//
// Scope: one interior face, two cells, both sides populated with distinct
// state/gradient/viscosity, and a second case with the face normal reversed so
// a signed-term error cannot hide behind a symmetric setup.  This is a
// component contract, not a substitute for the 3-step trace or the 200-step
// cross-backend gate.

#include "AccelRuntime.h"
#include "HipFluxBackend.h"
// The host reference links the production policy functions (NsVisFluxPolicy,
// UGradPolicy, Sutherland) so the expectation is the solver's own formula, and
// ViscousFluxArithmetic.h is the shared inline decomposition the device kernel
// compiles.  Comparing the two catches a change in either side.
#include "NsVisFluxPolicy.h"
#include "Sutherland.h"
#include "UGradPolicy.h"

#include <gtest/gtest.h>

#include <cmath>

using namespace ONEFLOW;

namespace
{

// The contract tolerance: 1e-14 relative to the face flux magnitude.  The host
// chain and the device kernel evaluate the same expression in the same order,
// so this is a round-off budget, not a model tolerance.
constexpr double kTolerance = 1.0e-14;

// One face of the host viscous chain, transcribed from
// UNsVisFlux::CalcFaceVisFlux with the VIS_STD preparation.
//
// Layout of the inputs mirrors ViscousFluxView: per-cell component arrays
// (primitive state, gradient, temperature, temperature gradient, viscosities)
// and per-face geometry.  Cell 0 is the left (lc) cell and cell 1 is the right
// (rc) cell, exactly as ug.lcf/ug.rcf index them.
struct HostViscousInputs
{
    Real q[ 5 ][ 2 ] = {};       // rho, u, v, w, p at the two cells
    Real dqdx[ 5 ][ 2 ] = {};
    Real dqdy[ 5 ][ 2 ] = {};
    Real dqdz[ 5 ][ 2 ] = {};
    Real temperature[ 2 ] = {};
    Real dTdx[ 2 ] = {};
    Real dTdy[ 2 ] = {};
    Real dTdz[ 2 ] = {};
    Real laminarViscosity[ 2 ] = {};
    Real turbulentViscosity[ 2 ] = {};

    Real xCell[ 2 ] = {};
    Real yCell[ 2 ] = {};
    Real zCell[ 2 ] = {};

    Real xfn = 0.0;
    Real yfn = 0.0;
    Real zfn = 0.0;
    Real faceArea = 0.0;

    Real gamma = 1.4;
    // Reciprocals of Pr = 0.72 and 0.9.  These were previously 0.72 and 0.90,
    // i.e. the Prandtl numbers themselves, which is exactly the mistake the
    // solver made when filling the view: the test and the call site agreed with
    // each other and both disagreed with the host, so the contract stayed green
    // while the viscous residual was off by 6.7e-3.
    Real inversePrandtlLaminar = 1.3888888888888889;
    Real inversePrandtlTurbulent = 1.1111111111111112;
    Real specificHeatCp = 0.0;
    Real inverseReynolds = 0.0;
};

// nscom.oprl/oprt are inverse Prandtl numbers (NsCom::Init sets oprl = 1/prl),
// so the view carries the inverse form directly.  The helper below makes the
// distinction explicit at the call site instead of burying it.
Real HostFaceViscousFlux(
    const HostViscousInputs & in, Real flux[ 5 ], Real & wallHeatFlux,
    Real & qNormalOut )
{
    const int lc = 0;
    const int rc = 1;

    // --- PrepareFaceValue: side reads and the two viscosity averages ---------
    const Real visl1 = in.laminarViscosity[ lc ];
    const Real visl2 = in.laminarViscosity[ rc ];
    const Real vist1 = in.turbulentViscosity[ lc ];
    const Real vist2 = in.turbulentViscosity[ rc ];
    const Real visl = half * ( visl1 + visl2 );
    const Real vist = half * ( vist1 + vist2 );
    const Real vis = visl + vist;

    // --- AverGrad ------------------------------------------------------------
    Real dqdx[ 5 ] = {};
    Real dqdy[ 5 ] = {};
    Real dqdz[ 5 ] = {};
    for ( int eq = 0; eq < 5; ++ eq )
    {
        dqdx[ eq ] = half * ( in.dqdx[ eq ][ lc ] + in.dqdx[ eq ][ rc ] );
        dqdy[ eq ] = half * ( in.dqdy[ eq ][ lc ] + in.dqdy[ eq ][ rc ] );
        dqdz[ eq ] = half * ( in.dqdz[ eq ][ lc ] + in.dqdz[ eq ][ rc ] );
    }

    // --- CalcGradCoef / UnitizeDisplacement ----------------------------------
    Real dx = in.xCell[ rc ] - in.xCell[ lc ];
    Real dy = in.yCell[ rc ] - in.yCell[ lc ];
    Real dz = in.zCell[ rc ] - in.zCell[ lc ];
    Real ods = 0.0;
    UnitizeDisplacement( dx, dy, dz, ods );

    // --- AverFaceValue -------------------------------------------------------
    Real qFace[ 5 ] = {};
    for ( int eq = 0; eq < 5; ++ eq )
    {
        qFace[ eq ] = ArithmeticFaceAverage( in.q[ eq ][ lc ], in.q[ eq ][ rc ] );
    }
    Real dTdx = ArithmeticFaceAverage( in.dTdx[ lc ], in.dTdx[ rc ] );
    Real dTdy = ArithmeticFaceAverage( in.dTdy[ lc ], in.dTdy[ rc ] );
    Real dTdz = ArithmeticFaceAverage( in.dTdz[ lc ], in.dTdz[ rc ] );

    // --- CorrectFaceGrad (visQ then visT) ------------------------------------
    // VisGrad::CorrectFaceGrad calls CorrectGrad, which is a one-line wrapper
    // around CorrectFaceGradient; the wrapper is skipped here because it lives
    // in VisGrad.cpp and pulls the VisGrad singleton into the link for no gain.
    for ( int eq = 0; eq < 5; ++ eq )
    {
        CorrectFaceGradient( in.q[ eq ][ lc ], in.q[ eq ][ rc ],
                             dqdx[ eq ], dqdy[ eq ], dqdz[ eq ],
                             dx, dy, dz, ods );
    }
    CorrectFaceGradient( in.temperature[ lc ], in.temperature[ rc ],
                         dTdx, dTdy, dTdz, dx, dy, dz, ods );

    // --- CalcNormalGrad ------------------------------------------------------
    const Real dTdn =
        in.xfn * dTdx + in.yfn * dTdy + in.zfn * dTdz;

    // --- CalcStress / StokesStressTensor -------------------------------------
    Real txx = 0.0;
    Real tyy = 0.0;
    Real tzz = 0.0;
    Real txy = 0.0;
    Real txz = 0.0;
    Real tyz = 0.0;
    StokesStressTensor(
        vis,
        dqdx[ 1 ], dqdy[ 1 ], dqdz[ 1 ],
        dqdx[ 2 ], dqdy[ 2 ], dqdz[ 2 ],
        dqdx[ 3 ], dqdy[ 3 ], dqdz[ 3 ],
        txx, tyy, tzz, txy, txz, tyz );

    // --- CalcNsVisFlux -------------------------------------------------------
    Real fvis[ 5 ] = {};
    fvis[ 0 ] = 0.0;
    ProjectStressMomentumFlux(
        in.xfn, in.yfn, in.zfn,
        txx, tyy, tzz, txy, txz, tyz,
        fvis[ 1 ], fvis[ 2 ], fvis[ 3 ] );

    // --- AddHeatFlux ---------------------------------------------------------
    // AddChemHeatFlux is empty for chemModel != 1 and FaceNormalDot on the
    // zeroed qx/qy/qz is zero, so this is the whole normal heat flux.
    const Real kcp = EffectiveThermalConductivity(
        visl, vist, in.inversePrandtlLaminar, in.inversePrandtlTurbulent,
        in.specificHeatCp );
    const Real qx = 0.0;
    const Real qy = 0.0;
    const Real qz = 0.0;
    Real qNormal = 0.0;
    qNormal += FaceNormalDot( in.xfn, in.yfn, in.zfn, qx, qy, qz );
    qNormal += ConductiveHeatFluxNormal( kcp, dTdn );
    qNormalOut = qNormal;

    fvis[ 4 ] = ViscousEnergyFlux(
        qFace[ 1 ], qFace[ 2 ], qFace[ 3 ],
        fvis[ 1 ], fvis[ 2 ], fvis[ 3 ], qNormal );

    // --- UpdateFaceVisFlux ---------------------------------------------------
    const Real coeff = ViscousFaceFluxCoefficient(
        in.inverseReynolds, in.faceArea );
    for ( int eq = 0; eq < 5; ++ eq ) flux[ eq ] = coeff * fvis[ eq ];

    // --- SaveHeatFlux --------------------------------------------------------
    wallHeatFlux = - in.inverseReynolds * qNormal;
    return qNormal;
}

// Fill a ViscousFluxView over the same one-face problem.  Component pointers
// point at the two-cell arrays of HostViscousInputs, so the view sees a
// "component array of two cells" per equation, which is the layout
// ViscousFluxView defines.
void FillView( const HostViscousInputs & in, ViscousFluxView & view,
               Real * flux, int stride )
{
    // The face centre is not read by the viscous chain (the host takes it from
    // ug.xfc only for the unused tangential decomposition), but a view has to
    // carry one, so use the midpoint rather than a cell centre.
    static Real xFace = 0.0;
    static Real yFace = 0.0;
    static Real zFace = 0.0;
    xFace = 0.5 * ( in.xCell[ 0 ] + in.xCell[ 1 ] );
    yFace = 0.5 * ( in.yCell[ 0 ] + in.yCell[ 1 ] );
    zFace = 0.5 * ( in.zCell[ 0 ] + in.zCell[ 1 ] );

    view.nCells = 2;
    view.nGhostCells = 0;
    view.nFaces = 1;
    view.nBoundaryFaces = 0;
    view.nEquations = 5;

    view.xFace = & xFace;
    view.yFace = & yFace;
    view.zFace = & zFace;
    view.xNormal = & in.xfn;
    view.yNormal = & in.yfn;
    view.zNormal = & in.zfn;
    view.faceArea = & in.faceArea;
    view.xCell = & in.xCell[ 0 ];
    view.yCell = & in.yCell[ 0 ];
    view.zCell = & in.zCell[ 0 ];
    static const Real volume[ 2 ] = { 1.0, 1.0 };
    view.cellVolume = volume;
    static const int left[ 1 ] = { 0 };
    static const int right[ 1 ] = { 1 };
    view.leftCell = left;
    view.rightCell = right;

    for ( int eq = 0; eq < 5; ++ eq )
    {
        view.q[ eq ] = & in.q[ eq ][ 0 ];
        view.dqdx[ eq ] = & in.dqdx[ eq ][ 0 ];
        view.dqdy[ eq ] = & in.dqdy[ eq ][ 0 ];
        view.dqdz[ eq ] = & in.dqdz[ eq ][ 0 ];
    }

    view.temperature = & in.temperature[ 0 ];
    view.dTdx = & in.dTdx[ 0 ];
    view.dTdy = & in.dTdy[ 0 ];
    view.dTdz = & in.dTdz[ 0 ];
    view.laminarViscosity = & in.laminarViscosity[ 0 ];
    view.turbulentViscosity = & in.turbulentViscosity[ 0 ];

    view.gamma = in.gamma;
    view.inversePrandtlLaminar = in.inversePrandtlLaminar;
    view.inversePrandtlTurbulent = in.inversePrandtlTurbulent;
    view.specificHeatCp = in.specificHeatCp;
    view.inverseReynolds = in.inverseReynolds;
    view.wallTemperature = 0.0;

    view.faceFlux = flux;
    view.faceFluxStride = stride;
}

// Largest |host component| on the face.  The viscous flux components are small
// and partly cancel, so a per-component relative test would divide by a value
// that is itself a cancellation residue.  Normalising the component differences
// by this scale keeps the comparison relative to the face flux magnitude while
// staying meaningful when one component happens to be near zero.
double HostFluxScale( const Real hostFlux[ 5 ] )
{
    double scale = 0.0;
    for ( int eq = 0; eq < 5; ++ eq )
    {
        const double magnitude = std::fabs( hostFlux[ eq ] );
        if ( magnitude > scale ) scale = magnitude;
    }
    return scale > 0.0 ? scale : 1.0;
}

// Compare the five components of one face against the host flux.
//
// deviceFlux is the equation-major output buffer and stride is the view's
// faceFluxStride, so component `eq` of this face lives at eq * stride + face.
// The view validates stride >= nEquations, so a small face count still needs a
// buffer of ( nEquations - 1 ) * stride + 1 entries; reading it as a plain
// five-element array would read past the end.
void ExpectFaceMatchesHost( const Real hostFlux[ 5 ],
                            const Real * deviceFlux,
                            int stride, int face,
                            const char * context )
{
    const double scale = HostFluxScale( hostFlux );
    for ( int eq = 0; eq < 5; ++ eq )
    {
        const Real device = deviceFlux[ eq * stride + face ];
        EXPECT_LE( std::fabs( device - hostFlux[ eq ] ) / scale,
                   kTolerance )
            << context << " equation " << eq
            << " host " << hostFlux[ eq ]
            << " device " << device;
    }
}

// A face whose normal is deliberately not axis-aligned and whose two cells are
// not mirror images, so an omitted sign or a swapped index changes the answer.
HostViscousInputs MakeInputs()
{
    HostViscousInputs in;

    // Primitive state: rho, u, v, w, p.
    const Real qLeft[ 5 ] = { 1.02, 0.31, -0.08, 0.05, 0.98 };
    const Real qRight[ 5 ] = { 0.94, 0.36, 0.02, -0.03, 0.91 };
    for ( int eq = 0; eq < 5; ++ eq )
    {
        in.q[ eq ][ 0 ] = qLeft[ eq ];
        in.q[ eq ][ 1 ] = qRight[ eq ];
    }

    // Distinct asymmetric gradients on the two sides.
    for ( int eq = 0; eq < 5; ++ eq )
    {
        const Real s = 0.1 * static_cast< Real >( eq + 1 );
        in.dqdx[ eq ][ 0 ] =  0.20 * s;
        in.dqdy[ eq ][ 0 ] = -0.07 * s;
        in.dqdz[ eq ][ 0 ] =  0.05 * s;
        in.dqdx[ eq ][ 1 ] = -0.11 * s;
        in.dqdy[ eq ][ 1 ] =  0.13 * s;
        in.dqdz[ eq ][ 1 ] = -0.04 * s;
    }

    in.temperature[ 0 ] = 1.02;
    in.temperature[ 1 ] = 0.96;
    in.dTdx[ 0 ] =  0.09;
    in.dTdy[ 0 ] = -0.03;
    in.dTdz[ 0 ] =  0.02;
    in.dTdx[ 1 ] = -0.05;
    in.dTdy[ 1 ] =  0.07;
    in.dTdz[ 1 ] =  0.01;

    in.laminarViscosity[ 0 ] = 1.31e-3;
    in.laminarViscosity[ 1 ] = 1.27e-3;
    in.turbulentViscosity[ 0 ] = 4.4e-3;
    in.turbulentViscosity[ 1 ] = 4.9e-3;

    in.xCell[ 0 ] = -0.50;
    in.yCell[ 0 ] = -0.25;
    in.zCell[ 0 ] = -0.10;
    in.xCell[ 1 ] =  0.60;
    in.yCell[ 1 ] =  0.35;
    in.zCell[ 1 ] =  0.20;

    // Not normalised on purpose: the host also uses the stored values, and the
    // projection terms are linear in them.
    in.xfn = 0.62;
    in.yfn = 0.55;
    in.zfn = 0.31;
    in.faceArea = 0.0873;

    in.gamma = 1.4;
    in.inversePrandtlLaminar = 1.3888888888888889;
    in.inversePrandtlTurbulent = 1.1111111111111112;
    // nscom.const_cp = 1 / ( ( gama_ref - 1 ) * mach_ref^2 ).
    in.specificHeatCp = 1.0 / ( ( 1.4 - 1.0 ) * 0.8395 * 0.8395 );
    in.inverseReynolds = 1.0 / 1.172e7;

    return in;
}

class DeviceViscousFluxContract : public ::testing::Test
{
protected:
    static void SetUpTestSuite()
    {
        InitializeAccelRuntime( 0, 1 );
        // A CPU backend here would allocate host memory for every DeviceBuffer
        // and hand it to a HIP kernel: the kernel faults, the copy-back reads
        // untouched memory, and the comparison can look like a pass.  Refuse
        // to run in that configuration instead.
        ASSERT_EQ( AccelRuntime::Instance().Backend().Kind(),
                   AccelBackendKind::HIP )
            << "viscous flux contract requires the HIP backend";
    }

    static void TearDownTestSuite()
    {
        AccelRuntime::Instance().Finalize();
    }
};

}

// The headline contract: one face, both sides populated, device flux equal to
// the host chain to 1e-14 relative.
TEST_F( DeviceViscousFluxContract, InteriorFaceMatchesHostChain )
{
    const HostViscousInputs in = MakeInputs();

    Real hostFlux[ 5 ] = {};
    Real hostWallHeatFlux = 0.0;
    Real hostQNormal = 0.0;
    HostFaceViscousFlux( in, hostFlux, hostWallHeatFlux, hostQNormal );

    constexpr int kStride = 5;
    Real deviceFlux[ ( 5 - 1 ) * kStride + 1 ] = {};
    ViscousFluxView view;
    FillView( in, view, deviceFlux, kStride );
    HipFluxBackend::Shared().CalcViscousFaceFlux( view );

    ExpectFaceMatchesHost( hostFlux, deviceFlux, kStride, 0, "interior face" );

    // The continuity row must be identically zero on both sides: a nonzero
    // value would mean momentum leaked into the density equation.
    EXPECT_EQ( hostFlux[ 0 ], 0.0 );
    EXPECT_EQ( deviceFlux[ 0 ], 0.0 );
}

// Reversing the face normal must flip every momentum component's sign and
// leave the energy term's conduction part on the same side of the sign
// convention as the host.  This is the case a symmetric fixture hides.
TEST_F( DeviceViscousFluxContract, ReversedNormalStillMatchesHost )
{
    HostViscousInputs in = MakeInputs();
    in.xfn = - in.xfn;
    in.yfn = - in.yfn;
    in.zfn = - in.zfn;

    Real hostFlux[ 5 ] = {};
    Real hostWallHeatFlux = 0.0;
    Real hostQNormal = 0.0;
    HostFaceViscousFlux( in, hostFlux, hostWallHeatFlux, hostQNormal );

    constexpr int kStride = 5;
    Real deviceFlux[ ( 5 - 1 ) * kStride + 1 ] = {};
    ViscousFluxView view;
    FillView( in, view, deviceFlux, kStride );
    HipFluxBackend::Shared().CalcViscousFaceFlux( view );

    ExpectFaceMatchesHost( hostFlux, deviceFlux, kStride, 0, "reversed normal" );
    EXPECT_EQ( deviceFlux[ 0 ], 0.0 );
}

// Laminar-only configuration (no SA solve): the view's turbulent viscosity is
// null and the device must give exactly the host answer with vist == 0.
TEST_F( DeviceViscousFluxContract, LaminarOnlyMatchesHost )
{
    // A laminar-only view has no turbulent field; the host chain then reads an
    // all-zero vist component, which the device models by a null pointer.
    HostViscousInputs laminarIn = MakeInputs();
    laminarIn.turbulentViscosity[ 0 ] = 0.0;
    laminarIn.turbulentViscosity[ 1 ] = 0.0;
    Real laminarHostFlux[ 5 ] = {};
    Real laminarWallHeatFlux = 0.0;
    Real laminarQNormal = 0.0;
    HostFaceViscousFlux( laminarIn, laminarHostFlux, laminarWallHeatFlux,
                         laminarQNormal );

    constexpr int kStride = 5;
    Real deviceFlux[ ( 5 - 1 ) * kStride + 1 ] = {};
    ViscousFluxView view;
    FillView( laminarIn, view, deviceFlux, kStride );
    view.turbulentViscosity = nullptr;
    HipFluxBackend::Shared().CalcViscousFaceFlux( view );

    ExpectFaceMatchesHost( laminarHostFlux, deviceFlux, kStride, 0,
                           "laminar only" );
}

// Wall heat flux has its own expression in SaveHeatFlux (-oreynolds * qNormal)
// and its own non-uniform layout.  The two-face, one-wall problem below checks
// both the value and that only the classified face is written.
TEST_F( DeviceViscousFluxContract, WallHeatFluxMatchesHost )
{
    const HostViscousInputs in = MakeInputs();

    Real hostFlux[ 5 ] = {};
    Real hostWallHeatFlux = 0.0;
    Real hostQNormal = 0.0;
    HostFaceViscousFlux( in, hostFlux, hostWallHeatFlux, hostQNormal );

    // Two identical faces so the kernel writes slot 0 and leaves slot 1 alone.
    Real q[ 5 ][ 4 ] = {};
    Real dqdx[ 5 ][ 4 ] = {};
    Real dqdy[ 5 ][ 4 ] = {};
    Real dqdz[ 5 ][ 4 ] = {};
    Real temperature[ 4 ] = {};
    Real dTdx[ 4 ] = {};
    Real dTdy[ 4 ] = {};
    Real dTdz[ 4 ] = {};
    Real laminar[ 4 ] = {};
    Real turbulent[ 4 ] = {};
    for ( int eq = 0; eq < 5; ++ eq )
    {
        q[ eq ][ 0 ] = in.q[ eq ][ 0 ]; q[ eq ][ 1 ] = in.q[ eq ][ 1 ];
        q[ eq ][ 2 ] = in.q[ eq ][ 0 ]; q[ eq ][ 3 ] = in.q[ eq ][ 1 ];
        dqdx[ eq ][ 0 ] = in.dqdx[ eq ][ 0 ]; dqdx[ eq ][ 1 ] = in.dqdx[ eq ][ 1 ];
        dqdx[ eq ][ 2 ] = in.dqdx[ eq ][ 0 ]; dqdx[ eq ][ 3 ] = in.dqdx[ eq ][ 1 ];
        dqdy[ eq ][ 0 ] = in.dqdy[ eq ][ 0 ]; dqdy[ eq ][ 1 ] = in.dqdy[ eq ][ 1 ];
        dqdy[ eq ][ 2 ] = in.dqdy[ eq ][ 0 ]; dqdy[ eq ][ 3 ] = in.dqdy[ eq ][ 1 ];
        dqdz[ eq ][ 0 ] = in.dqdz[ eq ][ 0 ]; dqdz[ eq ][ 1 ] = in.dqdz[ eq ][ 1 ];
        dqdz[ eq ][ 2 ] = in.dqdz[ eq ][ 0 ]; dqdz[ eq ][ 3 ] = in.dqdz[ eq ][ 1 ];
    }
    temperature[ 0 ] = in.temperature[ 0 ];
    temperature[ 1 ] = in.temperature[ 1 ];
    temperature[ 2 ] = in.temperature[ 0 ];
    temperature[ 3 ] = in.temperature[ 1 ];
    dTdx[ 0 ] = in.dTdx[ 0 ]; dTdx[ 1 ] = in.dTdx[ 1 ];
    dTdx[ 2 ] = in.dTdx[ 0 ]; dTdx[ 3 ] = in.dTdx[ 1 ];
    dTdy[ 0 ] = in.dTdy[ 0 ]; dTdy[ 1 ] = in.dTdy[ 1 ];
    dTdy[ 2 ] = in.dTdy[ 0 ]; dTdy[ 3 ] = in.dTdy[ 1 ];
    dTdz[ 0 ] = in.dTdz[ 0 ]; dTdz[ 1 ] = in.dTdz[ 1 ];
    dTdz[ 2 ] = in.dTdz[ 0 ]; dTdz[ 3 ] = in.dTdz[ 1 ];
    laminar[ 0 ] = in.laminarViscosity[ 0 ];
    laminar[ 1 ] = in.laminarViscosity[ 1 ];
    laminar[ 2 ] = in.laminarViscosity[ 0 ];
    laminar[ 3 ] = in.laminarViscosity[ 1 ];
    turbulent[ 0 ] = in.turbulentViscosity[ 0 ];
    turbulent[ 1 ] = in.turbulentViscosity[ 1 ];
    turbulent[ 2 ] = in.turbulentViscosity[ 0 ];
    turbulent[ 3 ] = in.turbulentViscosity[ 1 ];

    const Real xCell[ 4 ] = { in.xCell[ 0 ], in.xCell[ 1 ],
                              in.xCell[ 0 ], in.xCell[ 1 ] };
    const Real yCell[ 4 ] = { in.yCell[ 0 ], in.yCell[ 1 ],
                              in.yCell[ 0 ], in.yCell[ 1 ] };
    const Real zCell[ 4 ] = { in.zCell[ 0 ], in.zCell[ 1 ],
                              in.zCell[ 0 ], in.zCell[ 1 ] };
    const Real xNormal[ 2 ] = { in.xfn, in.xfn };
    const Real yNormal[ 2 ] = { in.yfn, in.yfn };
    const Real zNormal[ 2 ] = { in.zfn, in.zfn };
    const Real faceArea[ 2 ] = { in.faceArea, in.faceArea };
    const int left[ 2 ] = { 0, 2 };
    const int right[ 2 ] = { 1, 3 };

    ViscousFluxView view;
    view.nCells = 4;
    view.nGhostCells = 0;
    view.nFaces = 2;
    view.nBoundaryFaces = 1;
    view.nEquations = 5;
    view.xFace = xCell;
    view.yFace = yCell;
    view.zFace = zCell;
    view.xNormal = xNormal;
    view.yNormal = yNormal;
    view.zNormal = zNormal;
    view.faceArea = faceArea;
    view.xCell = xCell;
    view.yCell = yCell;
    view.zCell = zCell;
    static const Real volume[ 4 ] = { 1.0, 1.0, 1.0, 1.0 };
    view.cellVolume = volume;
    view.leftCell = left;
    view.rightCell = right;
    for ( int eq = 0; eq < 5; ++ eq )
    {
        view.q[ eq ] = q[ eq ];
        view.dqdx[ eq ] = dqdx[ eq ];
        view.dqdy[ eq ] = dqdy[ eq ];
        view.dqdz[ eq ] = dqdz[ eq ];
    }
    view.temperature = temperature;
    view.dTdx = dTdx;
    view.dTdy = dTdy;
    view.dTdz = dTdz;
    view.laminarViscosity = laminar;
    view.turbulentViscosity = turbulent;
    view.gamma = in.gamma;
    view.inversePrandtlLaminar = in.inversePrandtlLaminar;
    view.inversePrandtlTurbulent = in.inversePrandtlTurbulent;
    view.specificHeatCp = in.specificHeatCp;
    view.inverseReynolds = in.inverseReynolds;

    // Equation-major with stride 5 and two faces: the highest index the kernel
    // writes is ( 5 - 1 ) * 5 + 1 = 21, so the buffer must hold 22 entries.
    constexpr int kStride = 5;
    Real flux[ ( 5 - 1 ) * kStride + 2 ] = {};
    view.faceFlux = flux;
    view.faceFluxStride = kStride;

    // Only face 0 is a solid face.  The buffer has room for two slots so the
    // kernel has to respect nWallFaces rather than write both.
    Real wallHeatFlux[ 2 ] = { 0.0, 12345.0 };
    view.wallHeatFlux = wallHeatFlux;
    view.nWallFaces = 1;
    const int wallFaceIds[ 1 ] = { 0 };

    HipFluxBackend::Shared().CalcViscousFaceFlux( view, wallFaceIds );

    const double scale =
        std::max( std::fabs( hostWallHeatFlux ), 1.0e-300 );
    EXPECT_LE( std::fabs( wallHeatFlux[ 0 ] - hostWallHeatFlux ) / scale,
               kTolerance )
        << "host " << hostWallHeatFlux << " device " << wallHeatFlux[ 0 ];
    EXPECT_EQ( wallHeatFlux[ 1 ], 12345.0 )
        << "a face that is not a solid face must not be written";

    // The classified face's flux must still match the host as well.
    ExpectFaceMatchesHost( hostFlux, flux, kStride, 0, "wall face" );
}

// SutherlandViscosityNondim: the sqrt( t^3 ) form, compared against the host
// function itself rather than a re-derived formula.
TEST_F( DeviceViscousFluxContract, SutherlandMatchesHostNondimensionalForm )
{
    // ts = ts_dim / tref_dim from the m6 case: 110.4 / 300.
    const Real sutherlandTs = 110.4 / 300.0;
    const Real temperatures[ 4 ] = { 0.72, 0.96, 1.0, 1.35 };
    for ( int i = 0; i < 4; ++ i )
    {
        const Real host =
            SutherlandViscosityNondim( temperatures[ i ], sutherlandTs );
        const Real device = HipFluxBackend::Shared().EvaluateDeviceSutherland(
            temperatures[ i ], sutherlandTs );
        const double scale = std::max( std::fabs( host ), 1.0e-300 );
        EXPECT_LE( std::fabs( device - host ) / scale, kTolerance )
            << "temperature " << temperatures[ i ];
    }
}

// The base entry must fail closed when the caller asks for wall heat flux: the
// backend-neutral view has no boundary classification, and writing nothing
// would leave the solver's heat-flux vector stale.
TEST_F( DeviceViscousFluxContract, WallHeatFluxWithoutFaceIdsFailsClosed )
{
    const HostViscousInputs in = MakeInputs();
    Real deviceFlux[ 5 ] = {};
    ViscousFluxView view;
    FillView( in, view, deviceFlux, 5 );
    Real wallHeatFlux = 0.0;
    view.wallHeatFlux = & wallHeatFlux;
    view.nWallFaces = 1;

    EXPECT_THROW( HipFluxBackend::Shared().CalcViscousFaceFlux( view ),
                  std::invalid_argument );
}
