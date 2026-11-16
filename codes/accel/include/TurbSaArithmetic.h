/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\*---------------------------------------------------------------------------*/

#pragma once

// SaModelConstants lives in AccelViews.h, because the SA views carry one.  The
// dependency therefore runs from this header to that one and never the other
// way: AccelViews.h is reachable from every execution space, and a dependency
// in the opposite direction would put this arithmetic - and every kernel that
// uses it - into CUDA's and Kokkos's include graph for a plain POD.
#include "AccelViews.h"

#include <cmath>

// Host/device callable marker.
//
// These functions are the one implementation that both the host chain and the
// HIP kernels compile.  Without this attribute the device compiler rejects the
// call from a __global__ function ("reference to __host__ function ... in
// __global__ function"), and duplicating the formulas into device-only copies
// is exactly the transcription drift this header exists to prevent.
//
// The guard is written on __HIPCC__ rather than __CUDACC__ so the same header
// also works in a CUDA translation unit if that backend is ever enabled; on a
// plain host compiler the macro is empty.
#if defined( __HIPCC__ ) || defined( __CUDACC__ )
#define ONEFLOW_TURB_SA_HD __host__ __device__
#else
#define ONEFLOW_TURB_SA_HD
#endif

BeginNameSpace( ONEFLOW )

// Spalart-Allmaras one-equation closure arithmetic, shared by every execution
// space.
//
// These are the same formulas as the host chain, written as header-inline
// scalar functions so a device kernel and a host loop compile identical
// expressions instead of two transcriptions that can drift.  The host sources
// stay the reference:
//
//   codes/turb/src/TurbSaPolicy.cpp      closure terms (SaChi .. SaNegativePart)
//   codes/turb/src/TurbCom.cpp:438       TurbCom::CalcWorkVar
//   codes/turb/src/TurbCom.cpp:450       TurbCom::CalcSaProd
//   codes/turb/src/TurbCom.cpp:747       TurbCom::CalcSrcSa
//   codes/turb/src/TurbSrcFlux.cpp:282   TurbSrcFlux::CalcCellVist1Equ
//   codes/turb/src/TurbInvFlux.cpp:80    TurbInv::RoeFlux
//   codes/turb/src/TurbVisFluxPolicy.cpp face diffusivity and diffusive flux
//   codes/turb/src/TurbUpdate.cpp:64     TurbUpdate::CalcFlowField1Equ
//
// Nothing here allocates, reads a global, or touches a container: every input
// is a scalar or a fixed-size array.  No ONEFLOW math helper is called because
// those are ordinary host functions and are not device-callable; the few
// comparisons they stand in for are written out.  `sqrt`, `pow` and `exp` are
// unqualified for the same reason ViscousFluxArithmetic.h leaves them
// unqualified: the device toolchain resolves them from its own math headers.
//
// The SA variable is a transported scalar, not a conserved mean-flow equation.
// Its residual is `production - destruction + diffusion` with a source term and
// its update is a single scalar under-relaxation, so none of this belongs in
// the five-equation kernels or in UpdatePrimitiveStateKernel, which derives
// pressure from total energy.

namespace TurbSaArithmetic
{

// Constant.h defines SMALL = 1.0e-40 and the host closure divides by
// (x + SMALL) throughout.  The same value is required for bit parity.
constexpr Real kSmall = 1.0e-40;

ONEFLOW_TURB_SA_HD inline Real Absolute( Real value )
{
    // Same form as HXMath's ABS so the NaN branch matches.
    return ( value < 0.0 ) ? - value : value;
}

ONEFLOW_TURB_SA_HD inline Real Minimum( Real a, Real b )
{
    // HXMath::MIN(a, b) = (a < b) ? a : b.  Written out because the helper is
    // host-only.
    return ( a < b ) ? a : b;
}

ONEFLOW_TURB_SA_HD inline Real Maximum( Real a, Real b )
{
    return ( a > b ) ? a : b;
}

// SQR(a, b, c) evaluates (a*a + b*b) + c*c.
ONEFLOW_TURB_SA_HD inline Real Square3( Real a, Real b, Real c )
{
    return a * a + b * b + c * c;
}

// DIST(a, b, c) = sqrt(SQR(a, b, c)).
ONEFLOW_TURB_SA_HD inline Real Distance3( Real a, Real b, Real c )
{
    return sqrt( Square3( a, b, c ) );
}

// ---------------------------------------------------------------------------
// Eddy viscosity: TurbSrcFlux::CalcCellVist1Equ:282 -> SaComputeEddyViscosity
// (TurbSaPolicy.cpp:113) -> TurbCappedEddyViscosity (TurbViscosityPolicy.cpp:6).
// ---------------------------------------------------------------------------
struct EddyViscosityInput
{
    Real density = 0.0;
    Real workingVariable = 0.0;
    Real laminarViscosity = 0.0;
    SaModelConstants constants;
};

struct EddyViscosityResult
{
    Real chi = 0.0;
    Real fv1 = 0.0;
    Real uncappedEddyViscosity = 0.0;
    Real eddyViscosity = 0.0;
};

ONEFLOW_TURB_SA_HD inline EddyViscosityResult ComputeEddyViscosity( const EddyViscosityInput & in )
{
    const SaModelConstants & c = in.constants;
    // TurbSaPolicy.cpp:8  SaChi(nuet, rho/(visl + SMALL))
    const Real inverseKinematicViscosity =
        in.density / ( in.laminarViscosity + kSmall );
    EddyViscosityResult result;
    result.chi = in.workingVariable * inverseKinematicViscosity + kSmall;

    // TurbSaPolicy.cpp:13  fv1 = chi^3 / (chi^3 + cv1^3)
    const Real chi3 = result.chi * result.chi * result.chi;
    result.fv1 = chi3 / ( chi3 + c.cv1Cubed );

    // TurbSaPolicy.cpp:108  mu_t = rho * nuTilde * fv1
    result.uncappedEddyViscosity =
        in.density * in.workingVariable * result.fv1;

    // TurbViscosityPolicy.cpp:6  mu_t = min(mu_t, r_max * mu_l)
    result.eddyViscosity = Minimum(
        in.laminarViscosity * c.maxViscosityRatio,
        result.uncappedEddyViscosity );
    return result;
}

// ---------------------------------------------------------------------------
// Strain measure: TurbCom::CalcWorkVar:438.
//
// velocityGradient is the host `turbcom.dudx, dudy, dudz, dvdx, dvdy, dvdz,
// dwdx, dwdy, dwdz` order, i.e. dqdx/dqdy/dqdz of the mean-flow u, v, w.
// ---------------------------------------------------------------------------
struct StrainWork
{
    Real work1 = 0.0;   // ui,i
    Real work2 = 0.0;   // u1,1^2 + u2,2^2 + u3,3^2
    Real work3 = 0.0;   // off-diagonal strain
    Real work4 = 0.0;   // vorticity magnitude
};

ONEFLOW_TURB_SA_HD inline StrainWork ComputeStrainWork( const Real velocityGradient[ 9 ] )
{
    const Real dudx = velocityGradient[ 0 ];
    const Real dudy = velocityGradient[ 1 ];
    const Real dudz = velocityGradient[ 2 ];
    const Real dvdx = velocityGradient[ 3 ];
    const Real dvdy = velocityGradient[ 4 ];
    const Real dvdz = velocityGradient[ 5 ];
    const Real dwdx = velocityGradient[ 6 ];
    const Real dwdy = velocityGradient[ 7 ];
    const Real dwdz = velocityGradient[ 8 ];

    StrainWork work;
    work.work1 = dudx + dvdy + dwdz;
    work.work2 = Square3( dudx, dudy, dudz );
    work.work3 = Square3( dudy + dvdx, dvdz + dwdy, dwdx + dudz );
    work.work4 = Distance3( dwdy - dvdz, dudz - dwdx, dvdx - dudy );
    return work;
}

// TurbCom::CalcSaProd:450, reduced to the value CalcSrcSa actually reads.
//
// CalcSrcSa applies CalcSaProd and then assigns `str = work4`, so the strain
// rate the source uses is the *post-branch* work4 in every option.  For options
// 4 and 5 CalcSaProd stores the pre-branch work4 in `str` and overwrites
// work4 with the rotational correction, so `str = work4` yields the corrected
// value.  Option 1 stores the classical production in work4.  Option 2 (the
// flagship case) hits no branch, so work4 stays the vorticity magnitude.
ONEFLOW_TURB_SA_HD inline Real StrainProduction(
    const StrainWork & work, int productionOption )
{
    if ( productionOption == 1 )
    {
        // SaStrainProductionClassical: sqrt(|2 work2 + work3 - 2/3 work1^2|)
        const Real stp = 2.0 * work.work2 + work.work3
            - ( 2.0 / 3.0 ) * work.work1 * work.work1;
        return sqrt( Absolute( stp ) );
    }
    if ( productionOption == 4 )
    {
        // SaStrainProductionRotationalV4: omega - 2 max(0, omega - sqrt(|stp|))
        const Real stp = 2.0 * work.work2 + work.work3
            - ( 2.0 / 3.0 ) * work.work1 * work.work1;
        const Real omega = work.work4;
        const Real root = sqrt( Absolute( stp ) );
        return omega - 2.0 * Maximum( 0.0, omega - root );
    }
    if ( productionOption == 5 )
    {
        // SaStrainProductionRotationalV5 with |stp| fed in.
        const Real stp = Absolute(
            2.0 * work.work2 + work.work3
            - ( 2.0 / 3.0 ) * work.work1 * work.work1 );
        const Real omega = work.work4;
        return omega * ( 1.0 - Maximum( 0.0,
            ( omega * omega - stp ) / ( omega * omega + stp + kSmall ) ) );
    }
    return work.work4;
}

// ---------------------------------------------------------------------------
// Per-cell source term: TurbCom::CalcSrcSa:747.
//
// The host reads `rho`, `nuet`, `visl`, the nine mean-flow velocity gradients,
// the three gradients of the turbulence variable, `len_scale` (which is the
// wall distance unless a DES/DDES length scale is active) and `cvol`, and
// writes `turbcom.res_sa` and `turbcom.spec_sa`.
//
// Two host reads are not carried here because CalcSrcSa never uses them:
// `turbcom.vist` (loaded by PrepareCellValue1Equ:327) and `cross_term`
// (:329, a two-equation term).  `dft2dnu` is computed at TurbCom.cpp:810 but
// never read, so it is omitted; that changes no value.
// ---------------------------------------------------------------------------
struct SourceInput
{
    Real density = 0.0;
    Real workingVariable = 0.0;
    Real laminarViscosity = 0.0;
    Real velocityGradient[ 9 ] = {};
    Real workingVariableGradient[ 3 ] = {};
    Real lengthScale = 0.0;         // len_scale = wall distance for SA
    Real cellVolume = 0.0;
    SaModelConstants constants;
};

struct SourceResult
{
    Real residual = 0.0;        // res_sa
    Real implicitSpectrum = 0.0; // spec_sa
    // Diagnostic intermediates, in case a caller wants them; not part of the
    // host residual contract.
    Real fv1 = 0.0;
    Real fv2 = 0.0;
    Real fw = 0.0;
    Real modifiedStrain = 0.0;
};

ONEFLOW_TURB_SA_HD inline Real NegativePart( Real value )
{
    // TurbSaPolicy.cpp:103: -1/2 (x - |x|)
    return - 0.5 * ( value - Absolute( value ) );
}

ONEFLOW_TURB_SA_HD inline SourceResult ComputeSource( const SourceInput & in )
{
    const SaModelConstants & c = in.constants;

    // TurbCom.cpp:749-751
    const Real olam = in.density / ( in.laminarViscosity + kSmall );
    const Real d2 = in.lengthScale * in.lengthScale;
    const Real od2 = 1.0 / d2;

    // TurbCom.cpp:753-754
    const StrainWork work = ComputeStrainWork( in.velocityGradient );
    const Real omega = StrainProduction( work, c.productionOption );

    // TurbCom.cpp:759-762.  `compress` is set to 0.0 at :757 and never read.
    const Real xsi = in.workingVariable * olam + kSmall;
    const Real xsi3 = xsi * xsi * xsi;
    const Real fv1 = xsi3 / ( xsi3 + c.cv1Cubed );
    const Real fv2 = 1.0 - xsi / ( 1.0 + xsi * fv1 );

    // TurbCom.cpp:765-771, absorbing the Reynolds number into rs.
    const Real rs = in.workingVariable * c.okarm2 * od2 * c.invReynolds;
    const Real nuetRs = in.workingVariable * rs;
    const Real sBar = rs * fv2;

    // SaModifiedStrainRate, TurbSaPolicy.cpp:29.
    Real std = 0.0;
    if ( sBar >= - c.sac2 * omega )
    {
        std = omega + sBar;
    }
    else
    {
        const Real term1 = c.sac2 * c.sac2 * omega + c.sac3 * sBar;
        const Real term2 = ( c.sac3 - 2.0 * c.sac2 ) * omega - sBar;
        std = omega + omega * term1 / term2;
    }

    // TurbCom.cpp:773-781.  SaG/SaFw were rewritten to multiplicative integer
    // powers (ccb901b78); the intermediates they share with the Jacobian terms
    // are formed once here exactly as the host does.
    const Real ostd = 1.0 / ( std + kSmall );
    const Real r = Minimum( rs / ( std + kSmall ), 10.0 );
    const Real r2 = r * r;
    const Real r6 = r2 * r2 * r2;
    const Real r5 = r6 / r;
    const Real g = r + c.cw2 * ( r6 - r );
    const Real g2 = g * g;
    const Real g6 = g2 * g2 * g2;
    const Real fw = g * pow(
        ( 1.0 + c.cw3ToSix ) / ( g6 + c.cw3ToSix ), 1.0 / 6.0 );

    // TurbCom.cpp:783-790.  SaFt2, TurbSaPolicy.cpp:69.
    const Real xsi2 = xsi * xsi;
    const Real ft2 = c.ft2Enabled ? c.ct3 * exp( - c.ct4 * xsi2 ) : 0.0;

    const Real grd2 = Square3( in.workingVariableGradient[ 0 ],
                               in.workingVariableGradient[ 1 ],
                               in.workingVariableGradient[ 2 ] );

    // TurbCom.cpp:792-798: SaProduction, SaDiffusion, SaDestruction,
    // SaLimitedDiffusion, SaSourceResidual.
    const Real prod = c.cb1 * ( 1.0 - ft2 ) * std * in.workingVariable;
    Real diff = c.cb2OverSigma * grd2 * c.invReynolds;
    const Real dest = ( c.cw1K * fw - c.cb1 * ft2 ) * nuetRs;
    diff = Minimum( c.productionRatio * prod, diff );
    const Real srcTerm = prod - dest + diff;

    // TurbCom.cpp:800-823: the implicit part of the Jacobian, which the host
    // folds into spec_sa.
    const Real fv1Ratio = fv1 / ( xsi + kSmall );
    const Real fv1Denom = 1.0 + xsi * fv1;
    const Real dfv2dk =
        ( 3.0 * c.cv1Cubed * fv1Ratio * fv1Ratio - 1.0 )
        / ( fv1Denom * fv1Denom );
    const Real dsdnu =
        c.invReynolds * c.okarm2 * od2 * ( fv2 + xsi * dfv2dk );

    const Real dfwdg = fw / ( g + kSmall ) * ( 1.0 - g6 / ( g6 + c.cw3ToSix ) );
    const Real dgdr = 1.0 + c.cw2 * ( 6.0 * r5 - 1.0 );
    const Real drdnu = c.invReynolds * c.okarm2 * od2 * ostd
        * ( 1.0 - in.workingVariable * ostd * dsdnu );
    const Real dfwdnu = dfwdg * dgdr * drdnu;

    const Real prod0 = c.cb1 * ( 1.0 - ft2 ) * std;
    const Real dest0 = ( c.cw1 * fw - c.cb1 * ft2 * c.okarm2 )
        * in.workingVariable * od2 * c.invReynolds;

    const Real prodp = c.cb1 * ( ( 1.0 - ft2 ) * dsdnu );
    const Real destp = ( c.cw1 * ( fw + dfwdnu * in.workingVariable )
        - c.cb1 * ft2 * c.okarm2 ) * od2 * c.invReynolds;

    const Real prde0 = prod0 - dest0;
    const Real prdep = prodp - destp;

    SourceResult result;
    result.residual = srcTerm;
    result.implicitSpectrum = NegativePart( prde0 )
        + NegativePart( prdep ) * in.workingVariable;
    result.fv1 = fv1;
    result.fv2 = fv2;
    result.fw = fw;
    result.modifiedStrain = std;
    return result;
}

// ---------------------------------------------------------------------------
// Scalar upwind face flux: TurbInv::RoeFlux:80 and
// UTurbInvFlux::UpdateFaceInvFlux:164.
//
// For SA (vismodel = 3) the host sets rho_coef = 1, so
// TurbConvectiveDensityWeight returns exactly 1 and the flux is the plain
// upwind scalar transport.  rhoCoef is still passed so the same kernel serves
// the two-equation case (rho_coef = 0).
// ---------------------------------------------------------------------------
struct InviscidFaceInput
{
    Real workingVariableLeft = 0.0;
    Real workingVariableRight = 0.0;
    Real densityLeft = 0.0;
    Real densityRight = 0.0;
    Real velocityLeft[ 3 ] = {};    // reconstructed u, v, w
    Real velocityRight[ 3 ] = {};
    Real normal[ 3 ] = {};
    Real meshVelocityNormal = 0.0;
    Real faceArea = 0.0;
    Real densityWeight = 1.0;       // rho_coef: 1 for SA, 0 for the 2eq pair
};

ONEFLOW_TURB_SA_HD inline Real SplitPositive( Real relativeNormalVelocity )
{
    // TurbInvFluxPolicy.cpp:6, vn+ = 1/2 (vn + |vn|)
    return 0.5 * ( relativeNormalVelocity + Absolute( relativeNormalVelocity ) );
}

ONEFLOW_TURB_SA_HD inline Real SplitNegative( Real relativeNormalVelocity )
{
    return 0.5 * ( relativeNormalVelocity - Absolute( relativeNormalVelocity ) );
}

ONEFLOW_TURB_SA_HD inline Real ConvectiveDensityWeight( Real densityWeight, Real density )
{
    // TurbInvFluxPolicy.cpp:13, w = alpha + (1 - alpha) rho.
    return densityWeight * 1.0 + ( 1.0 - densityWeight ) * density;
}

ONEFLOW_TURB_SA_HD inline Real ComputeInviscidFaceFlux( const InviscidFaceInput & in )
{
    // NsInvFluxPolicy.cpp:7, vn = n . V - v_mesh.
    const Real vnl = in.normal[ 0 ] * in.velocityLeft[ 0 ]
        + in.normal[ 1 ] * in.velocityLeft[ 1 ]
        + in.normal[ 2 ] * in.velocityLeft[ 2 ]
        - in.meshVelocityNormal;
    const Real vnr = in.normal[ 0 ] * in.velocityRight[ 0 ]
        + in.normal[ 1 ] * in.velocityRight[ 1 ]
        + in.normal[ 2 ] * in.velocityRight[ 2 ]
        - in.meshVelocityNormal;

    const Real vnlPositive = SplitPositive( vnl );
    const Real vnrNegative = SplitNegative( vnr );

    const Real weightLeft =
        ConvectiveDensityWeight( in.densityWeight, in.densityLeft );
    const Real weightRight =
        ConvectiveDensityWeight( in.densityWeight, in.densityRight );

    // TurbInvFluxPolicy.cpp:18, F = vn+ w_L phi_L + vn- w_R phi_R.
    const Real flux = vnlPositive * weightLeft * in.workingVariableLeft
        + vnrNegative * weightRight * in.workingVariableRight;

    return in.faceArea * flux;
}

// ---------------------------------------------------------------------------
// Scalar diffusive face flux: UTurbVisFlux::CalcFaceVisFlux1Equ:338.
//
// Two host details matter for parity and are easy to get wrong:
//
//  1. The turbulent diffusivity is `ArithmeticMean(q1[IKE], q2[IKE])`, and for
//     the one-equation model IKE == ISA == 0, so this is the mean of the
//     *working variable*, not of the eddy viscosity.  `turbcom.vist1/vist2` are
//     loaded at :325-326 but are unused on this branch.  Mirroring the host
//     means using the mean working variable.
//  2. The laminar part is the mean of `visl/rho` per side, not `visl` over the
//     mean density.
// ---------------------------------------------------------------------------
struct DiffusiveFaceInput
{
    // Cell-centred working variable at each side.  The host fills visTurb.q1/
    // q2 from `(*uturbf.q)[IKE][lc/rc]` (UTurbVisFlux.cpp:303-311), and for the
    // one-equation model IKE == ISA, so these are the cell values, not the
    // reconstructed face values.
    Real workingVariableLeft = 0.0;
    Real workingVariableRight = 0.0;
    // Mean-flow laminar viscosity and density per side.  The host forms the
    // kinematic value as visl * (1/rho) in CalcFaceVisFlux1Equ:339.
    Real laminarViscosityLeft = 0.0;
    Real laminarViscosityRight = 0.0;
    Real densityLeft = 0.0;
    Real densityRight = 0.0;
    Real normalGradient = 0.0;      // d(nuet)/dn at the face
    Real faceArea = 0.0;
    Real invSigma = 1.5;
    Real invReynolds = 0.0;
};

ONEFLOW_TURB_SA_HD inline Real ComputeEffectiveDiffusivity( const DiffusiveFaceInput & in )
{
    const Real orl = 1.0 / ( in.densityLeft + kSmall );
    const Real orr = 1.0 / ( in.densityRight + kSmall );
    // ArithmeticMean of the per-side kinematic laminar viscosity, exactly as
    // CalcFaceVisFlux1Equ:339 writes `ArithmeticMean(visl1 * orl, visl2 * orr)`.
    const Real kinematicLaminar =
        0.5 * ( in.laminarViscosityLeft * orl + in.laminarViscosityRight * orr );
    const Real workingVariableMean =
        0.5 * ( in.workingVariableLeft + in.workingVariableRight );
    // SaFaceEffectiveDiffusivity, TurbVisFluxPolicy.cpp:11.
    return ( kinematicLaminar + workingVariableMean ) * in.invSigma;
}

ONEFLOW_TURB_SA_HD inline Real ComputeDiffusiveFaceFlux( const DiffusiveFaceInput & in )
{
    // TurbScalarDiffusiveFlux, TurbVisFluxPolicy.cpp:25:
    // -1/Re * effectiveDiffusivity * d(nuet)/dn * area.
    return - in.invReynolds * ComputeEffectiveDiffusivity( in )
        * in.normalGradient * in.faceArea;
}

// ---------------------------------------------------------------------------
// Explicit update: TurbUpdate::CalcFlowField1Equ:64 through
// SaUnderRelaxedUpdate (TurbSaPolicy.cpp:169).
//
// TurbUpdate::ModifyValue1Equ:73 (the non-negative floor and the
// max-viscosity-ratio cap) has no caller anywhere in the tree, so the host
// update that actually runs is the under-relaxed sum alone.  It is not mirrored
// here; adding a live clamp on the device would make the device and host
// diverge by exactly the clamp.
// ---------------------------------------------------------------------------
ONEFLOW_TURB_SA_HD inline Real UnderRelaxedUpdate( Real qOld, Real deltaQ, Real relaxCoeff )
{
    const Real qCandidate = qOld + deltaQ;
    return ( 1.0 - relaxCoeff ) * qOld + relaxCoeff * qCandidate;
}

// Upper bound on the working variable, from SaWorkingVariableUpperBound
// (TurbSaPolicy.cpp:151), kept for the future boundary/limit stage.
ONEFLOW_TURB_SA_HD inline Real WorkingVariableUpperBound(
    Real maxViscosityRatio, Real laminarViscosity,
    Real density, Real fv1 )
{
    return maxViscosityRatio * laminarViscosity / ( density * fv1 );
}

// ---------------------------------------------------------------------------
// Cell gradient: Green-Gauss with inverse-distance face weights.
//
// Mirrors CalcGradGGCellWeight (codes/usolver/src/UGrad.cpp:88), which is what
// UTurbGrad::CalcGrad reaches through Grad::CalcGrad for the turbulence
// variable.  The three host passes are kept separate here because a device
// kernel cannot do them at different parallelism:
//
//   1. FaceAccumulateGradient -- left cell gets +S*value, the right cell of an
//      *internal* face gets -S*value.  Boundary faces only touch the left cell
//      (UGrad.cpp:121 `continue`), so a boundary ghost never accumulates.
//   2. ScaleGradientByVolume -- divide by the owner cell volume.
//   3. CopyGradient -- boundary ghost inherits the left cell gradient.
// ---------------------------------------------------------------------------
ONEFLOW_TURB_SA_HD inline void InverseDistanceWeights(
    Real dxl, Real dyl, Real dzl,
    Real dxr, Real dyr, Real dzr,
    Real & weightLeft, Real & weightRight )
{
    // UGradPolicy.cpp:11
    const Real delt1 = Distance3( dxl, dyl, dzl );
    const Real delt2 = Distance3( dxr, dyr, dzr );
    const Real delta = 1.0 / ( delt1 + delt2 + kSmall );
    weightLeft = delt2 * delta;
    weightRight = delt1 * delta;
}

// One face's contribution to the gradient.  `sign` is +1 for the left cell and
// -1 for an internal face's right cell; UGradPolicy.cpp:67.
ONEFLOW_TURB_SA_HD inline void AccumulateFaceGradient(
    Real snx, Real sny, Real snz, Real faceValue, Real sign,
    Real & dqdx, Real & dqdy, Real & dqdz )
{
    const Real w = sign * faceValue;
    dqdx += snx * w;
    dqdy += sny * w;
    dqdz += snz * w;
}

ONEFLOW_TURB_SA_HD inline void ScaleGradientByVolume(
    Real volume, Real & dqdx, Real & dqdy, Real & dqdz )
{
    // UGradPolicy.cpp:47
    const Real ovol = 1.0 / volume;
    dqdx *= ovol;
    dqdy *= ovol;
    dqdz *= ovol;
}

// ---------------------------------------------------------------------------
// Face tangents: GCom::CalcTangent, codes/global/src/Com.cpp:66.
//
// Mirrors the host branch order and its degenerate case exactly; the second
// tangent is the cross product of the normal with the first.
// ---------------------------------------------------------------------------
ONEFLOW_TURB_SA_HD inline bool ComputeFaceTangents(
    const Real normal[ 3 ], Real tangent1[ 3 ], Real tangent2[ 3 ] )
{
    const Real xfn = normal[ 0 ];
    const Real yfn = normal[ 1 ];
    const Real zfn = normal[ 2 ];

    if ( xfn != 0.0 )
    {
        tangent1[ 0 ] = yfn;
        tangent1[ 1 ] = - xfn;
        tangent1[ 2 ] = 0.0;
    }
    else if ( yfn != 0.0 )
    {
        tangent1[ 0 ] = - yfn;
        tangent1[ 1 ] = xfn;
        tangent1[ 2 ] = 0.0;
    }
    else if ( zfn != 0.0 )
    {
        tangent1[ 0 ] = 0.0;
        tangent1[ 1 ] = - zfn;
        tangent1[ 2 ] = yfn;
    }
    else
    {
        tangent1[ 0 ] = 0.0;
        tangent1[ 1 ] = 0.0;
        tangent1[ 2 ] = 0.0;
        tangent2[ 0 ] = 0.0;
        tangent2[ 1 ] = 0.0;
        tangent2[ 2 ] = 0.0;
        return false;
    }

    const Real inverseLength = 1.0 / Distance3(
        tangent1[ 0 ], tangent1[ 1 ], tangent1[ 2 ] );
    tangent1[ 0 ] *= inverseLength;
    tangent1[ 1 ] *= inverseLength;
    tangent1[ 2 ] *= inverseLength;

    tangent2[ 0 ] = yfn * tangent1[ 2 ] - zfn * tangent1[ 1 ];
    tangent2[ 1 ] = zfn * tangent1[ 0 ] - xfn * tangent1[ 2 ];
    tangent2[ 2 ] = xfn * tangent1[ 1 ] - yfn * tangent1[ 0 ];
    return true;
}

// ---------------------------------------------------------------------------
// Face normal gradient of the turbulence variable.
//
// The flagship case sets iturb_visflux = 2, which the host maps to
// UTurbVisFlux::TestMethod (UTurbVisFlux.cpp:88), and that variant is the one
// mirrored here:
//
//   VisGradGeom::PrepareCellGeom (VisGrad.cpp:57)  d1, d2, skew angles
//   VisGrad::CalcTestMethod       (VisGrad.cpp:264) face-value correction
//   VisGrad::CalcC1C2             (VisGrad.cpp:206) weight, with the
//                                                   non-interface boundary
//                                                   face forced to one side
//   VisGrad::ModifyFaceGrad       (VisGrad.cpp:241) tangential blend
//
// The two gradient-ratio weights fw1/fw2 come from VisGradGeom::CalcFaceWeight
// (VisGrad.cpp:35), which the host calls for every variant through
// UTurbVisFlux::PrepareFaceValue, and the tangents come from
// GCom::CalcTangent.  `side1`/`side2` are the two cell gradients at the face's
// left and right cell.
//
// `isInterface` selects the host's CalcC1C2 exception: for a boundary face
// whose bc type is not INTERFACE the host sets c1 = 1, c2 = 0, i.e. it uses the
// left side's one-sided normal gradient.  That is the value the flagship's
// solid and farfield faces take.
// ---------------------------------------------------------------------------
struct FaceNormalGradientInput
{
    Real workingVariableLeft = 0.0;
    Real workingVariableRight = 0.0;
    Real gradientLeft[ 3 ] = {};
    Real gradientRight[ 3 ] = {};
    Real xFace = 0.0;
    Real yFace = 0.0;
    Real zFace = 0.0;
    Real xCellLeft = 0.0;
    Real yCellLeft = 0.0;
    Real zCellLeft = 0.0;
    Real xCellRight = 0.0;
    Real yCellRight = 0.0;
    Real zCellRight = 0.0;
    Real normal[ 3 ] = {};
    bool boundary = false;
    bool isInterface = false;
};

struct FaceNormalGradientResult
{
    Real normalGradient = 0.0;
    Real gradient[ 3 ] = {};   // the reconstructed face gradient
};

ONEFLOW_TURB_SA_HD inline FaceNormalGradientResult ComputeFaceNormalGradient(
    const FaceNormalGradientInput & in )
{
    const Real xfn = in.normal[ 0 ];
    const Real yfn = in.normal[ 1 ];
    const Real zfn = in.normal[ 2 ];

    // VisGradGeom::CalcFaceWeight: inverse-distance weights on the same
    // centre-to-face displacements the Green-Gauss weights use.
    Real fw1 = 0.0;
    Real fw2 = 0.0;
    InverseDistanceWeights(
        in.xFace - in.xCellLeft, in.yFace - in.yCellLeft,
        in.zFace - in.zCellLeft,
        in.xFace - in.xCellRight, in.yFace - in.yCellRight,
        in.zFace - in.zCellRight,
        fw1, fw2 );

    // VisGradGeom::PrepareCellGeom: signed distance from each cell centre to
    // the face along the normal, plus the skew displacements.
    const Real dxl = in.xCellLeft - in.xFace;
    const Real dyl = in.yCellLeft - in.yFace;
    const Real dzl = in.zCellLeft - in.zFace;
    const Real dxr = in.xCellRight - in.xFace;
    const Real dyr = in.yCellRight - in.yFace;
    const Real dzr = in.zCellRight - in.zFace;

    const Real d1 = xfn * dxl + yfn * dyl + zfn * dzl;
    const Real d2 = xfn * dxr + yfn * dyr + zfn * dzr;

    const Real dxnl = xfn * d1 - dxl;
    const Real dynl = yfn * d1 - dyl;
    const Real dznl = zfn * d1 - dzl;
    const Real dxnr = xfn * d2 - dxr;
    const Real dynr = yfn * d2 - dyr;
    const Real dznr = zfn * d2 - dzr;

    // The host stores the two skew angles through FaceSkewAngleDegrees
    // (UGradPolicy.cpp:84), i.e. asin(clamp(proj/|d|, -1, 1)) * 180/pi, with
    // proj = -d1 on the left and +d2 on the right.  Two distinct host tests
    // read them:
    //
    //   FaceAngleIsValid  angle1 > skewAngle (20 degrees) and angle2 > 20
    //   TestSatisfied     angle1 > 0 and angle2 > 0 and |d1|,|d2| > SMALL
    //
    // Because asin is monotone, comparing angle against 20 degrees is the same
    // as comparing proj/|d| against sin(20 degrees), so the asin and the
    // degrees conversion are not needed and the test stays exact.
    const Real skewSine = 0.34202014332566871;   // sin( 20 degrees )
    const Real lengthLeft = Distance3( dxl, dyl, dzl );
    const Real lengthRight = Distance3( dxr, dyr, dzr );
    const Real ratioLeft = - d1 / ( lengthLeft + kSmall );
    const Real ratioRight = d2 / ( lengthRight + kSmall );
    const bool angleValid = ratioLeft > skewSine && ratioRight > skewSine;
    const bool signedAnglesPositive = ratioLeft > 0.0 && ratioRight > 0.0;

    // VisGrad::AverFaceValue (VisGrad.cpp:144).
    const Real leftValue = in.workingVariableLeft;
    const Real rightValue = in.workingVariableRight;
    Real faceValue = 0.5 * ( leftValue + rightValue );

    // VisGrad::AccurateSideValue + AccurateFaceValue (VisGrad.cpp:224).
    Real q11 = leftValue;
    Real q22 = rightValue;
    if ( angleValid )
    {
        q11 += in.gradientLeft[ 0 ] * dxnl + in.gradientLeft[ 1 ] * dynl
            + in.gradientLeft[ 2 ] * dznl;
        q22 += in.gradientRight[ 0 ] * dxnr + in.gradientRight[ 1 ] * dynr
            + in.gradientRight[ 2 ] * dznr;
        // The host recomputes `q` from the untouched q1/q2, so faceValue is
        // unchanged; the assignment is kept for clarity only.
        faceValue = 0.5 * ( leftValue + rightValue );
    }

    // VisGrad::TestSatisfied -> CalcC1C2 (VisGrad.cpp:178, :206).  The host
    // requires both signed skew angles to be positive (not the 20-degree test
    // above) plus |d1|,|d2| > SMALL.
    Real c1 = 0.0;
    Real c2 = 0.0;
    const bool oneSided = signedAnglesPositive
        && Absolute( d1 ) > kSmall && Absolute( d2 ) > kSmall;
    if ( oneSided )
    {
        const Real dtmp = d1 * d1 + d2 * d2;
        c1 = d1 * d1 / dtmp;
        c2 = 1.0 - c1;
        if ( in.boundary && ! in.isInterface )
        {
            c1 = 1.0;
            c2 = 0.0;
        }
    }

    // VisGrad::CalcTestMethod: the one-sided normal gradients, then the blend.
    Real normalGradient = 0.0;
    if ( oneSided )
    {
        const Real dqdn1 = ( q11 - faceValue ) / d1;
        const Real dqdn2 = ( q22 - faceValue ) / d2;
        normalGradient = c1 * dqdn1 + c2 * dqdn2;
    }

    // VisGrad::ModifyFaceGrad: inverse-distance gradient blend, then rebuild
    // the gradient from the normal and the two tangential components.
    Real blended[ 3 ];
    for ( int component = 0; component < 3; ++ component )
    {
        blended[ component ] = fw1 * in.gradientLeft[ component ]
            + fw2 * in.gradientRight[ component ];
    }

    Real tangent1[ 3 ];
    Real tangent2[ 3 ];
    const bool tangents = ComputeFaceTangents( in.normal, tangent1, tangent2 );

    FaceNormalGradientResult result;
    result.normalGradient = normalGradient;
    if ( ! tangents )
    {
        for ( int component = 0; component < 3; ++ component )
        {
            result.gradient[ component ] = blended[ component ];
        }
        return result;
    }

    const Real dqdt1 = tangent1[ 0 ] * blended[ 0 ]
        + tangent1[ 1 ] * blended[ 1 ] + tangent1[ 2 ] * blended[ 2 ];
    const Real dqdt2 = tangent2[ 0 ] * blended[ 0 ]
        + tangent2[ 1 ] * blended[ 1 ] + tangent2[ 2 ] * blended[ 2 ];

    result.gradient[ 0 ] = xfn * normalGradient
        + tangent1[ 0 ] * dqdt1 + tangent2[ 0 ] * dqdt2;
    result.gradient[ 1 ] = yfn * normalGradient
        + tangent1[ 1 ] * dqdt1 + tangent2[ 1 ] * dqdt2;
    result.gradient[ 2 ] = zfn * normalGradient
        + tangent1[ 2 ] * dqdt1 + tangent2[ 2 ] * dqdt2;
    return result;
}

}

EndNameSpace
