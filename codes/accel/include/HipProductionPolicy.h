/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\*---------------------------------------------------------------------------*/

#pragma once

#include "AccelRuntime.h"

#include <cstdlib>
#include <cstring>
#include <iostream>

BeginNameSpace( ONEFLOW )

// Policy flags are queried from the hottest paths in the solver: every stage
// asks HipProductionFeatureEnabled / HipDeviceAuthoritativeStateEnabled /
// HipSkipHostInterfaceExchange, and the interface field loop asks once per
// field.  std::getenv scans the whole environment block on each call, so a
// single RK step was paying hundreds of scans.  Resolve each queried name once
// and return the cached pointer; a flag is a process-lifetime decision here.
inline const char * CachedEnv( const char * name )
{
    struct Entry
    {
        const char * name;
        const char * value;
    };
    static Entry entries[ 32 ];
    static int count = 0;
    for ( int i = 0; i < count; ++ i )
    {
        if ( std::strcmp( entries[ i ].name, name ) == 0 )
            return entries[ i ].value;
    }
    const char * value = std::getenv( name );
    if ( count < 32 )
    {
        entries[ count ].name = name;
        entries[ count ].value = value;
        ++ count;
    }
    return value;
}

// Production HIP follows the runtime backend selection. Legacy per-stage
// variables remain overrides: unset selects automatically, "0" disables the
// stage, and "1" requests it subject to capability checks.
inline bool HipProductionBackendSelected()
{
#ifdef ONEFLOW_ENABLE_HIP
    // Asked by several policy helpers on every stage.  The runtime is created
    // once per process and its kind never changes afterwards, so cache the
    // positive answer and keep re-checking until the runtime exists.
    static bool confirmedHip = false;
    if ( confirmedHip ) return true;
    AccelRuntime & runtime = AccelRuntime::Instance();
    const bool hip = runtime.IsInitialized()
        && runtime.Backend().Kind() == AccelBackendKind::HIP;
    if ( hip ) confirmedHip = true;
    return hip;
#else
    return false;
#endif
}

inline bool HipProductionFeatureExplicitlyEnabled( const char * name )
{
    const char * value = CachedEnv(  name  );
    return value != nullptr && value[ 0 ] == '1';
}

inline bool HipProductionFeatureEnabled( const char * legacyName )
{
    const char * value = CachedEnv(  legacyName  );
    if ( value != nullptr )
    {
        // An explicit =1 still has to have a HIP backend behind it.  Without
        // this check the hook launches HIP kernels over host pointers, which is
        // the same failure as a missing upload but harder to attribute: the
        // buffers come from CpuBackend's plain operator new.
        return value[ 0 ] == '1' && HipProductionBackendSelected();
    }
    if ( ! HipProductionBackendSelected() ) return false;

    // Stages enabled by `ONEFLOW_ACCEL_BACKEND=HIP` alone.  Every one of them
    // has passed the same target-node gates (3-step 36-record trace against
    // the legacy oracle plus fixed-CFL 50-step stability, single- and
    // multi-zone, and the resident-state contract at 0 difference):
    //   * fused primitive flux/residual,
    //   * the state-owned interface/halo exchange,
    //   * device Green-Gauss gradient,
    //   * device reconstruction (cross-zone ghost gradient still travels
    //     through the host-staged exchange contract),
    //   * device primitive state update plus the resident residual chain,
    //   * device local boundary model,
    //   * resident local time-step.
    // The per-stage UseHip*/Supported guards still fall back to the host path
    // when the configuration is outside the validated scope (viscous,
    // limited, multi-grid, non-Lax-Friedrichs, more than one local grid, ...),
    // and an explicitly requested but unsupported stage still fails fast.
    return std::strcmp( legacyName, "ONEFLOW_ENABLE_UNS_HIP_BATCH" ) == 0
        || std::strcmp( legacyName, "ONEFLOW_ENABLE_UNS_HIP_GRADIENT" ) == 0
        || std::strcmp( legacyName, "ONEFLOW_ENABLE_UNS_HIP_INTERFACE" ) == 0
        || std::strcmp( legacyName, "ONEFLOW_ENABLE_UNS_HIP_RECONSTRUCTION" ) == 0
        || std::strcmp( legacyName, "ONEFLOW_ENABLE_UNS_HIP_STATE_UPDATE" ) == 0
        || std::strcmp( legacyName, "ONEFLOW_ENABLE_UNS_HIP_BOUNDARY" ) == 0
        || std::strcmp( legacyName, "ONEFLOW_ENABLE_UNS_HIP_TIME_STEP" ) == 0;
}

// Device-packed gradient interface exchange.  The send side packs the
// resident device gradient directly into the interface storage order; the
// receive side unpacks the interface storage directly into the device ghost
// gradient.  This removes the host cell gather/scatter and the separate
// ghost-gradient upload while preserving the legacy one-round lag.  The host
// mirror of boundary-adjacent owned cells is still produced for aero forces.
// Skip the legacy host dataBook interface exchange when the resident device
// exchange already owns the halo.  The device path packs/sends/unpacks the
// same cells, and the per-iteration post-process round trip still refreshes
// the host mirror, so the extra host MPI pass is pure duplication.  Set
// ONEFLOW_HIP_SKIP_HOST_INTERFACE=0 to restore the legacy double exchange for
// bisection; the 200-step cross-backend gate covers the default.
inline bool HipSkipHostInterfaceExchange()
{
    const char * value = CachedEnv(  "ONEFLOW_HIP_SKIP_HOST_INTERFACE"  );
    if ( value != nullptr ) return value[ 0 ] != '0';
    return HipProductionBackendSelected();
}

// The host mirror of the boundary-adjacent gradient only feeds aerodynamic
// forces, the optional boundary-gradient store and field dumps.  Those
// readers run after the Runge-Kutta loop, so an intermediate stage's mirror
// is overwritten before anyone reads it.  Mirroring every stage costs a
// device pack plus one D2H of the sampled cells.  Request it only on steps
// whose host consumers actually run, and only on the last stage of those
// steps.  ONEFLOW_HIP_GRADIENT_MIRROR_ON_DEMAND=0 restores the per-stage
// mirror for bisection.
inline bool HipGradientHostMirrorOnDemand()
{
    const char * value = CachedEnv(  "ONEFLOW_HIP_GRADIENT_MIRROR_ON_DEMAND"  );
    return value == nullptr || value[ 0 ] != '0';
}

// rkStage < 0 means the caller is not inside a staged loop.  A single-stage
// or unstaged evaluation still mirrors when a post-step reader exists.
inline bool HipGradientMirrorOnThisStage(
    bool mirrorOnDemand,
    bool forceEveryStage,
    bool postStepReader,
    int rkStage,
    int stageCount )
{
    if ( ! mirrorOnDemand || forceEveryStage ) return true;
    if ( ! postStepReader ) return false;
    if ( rkStage < 0 || stageCount <= 1 ) return true;
    return rkStage + 1 >= stageCount;
}

inline bool HipPackedGradientExchangeEnabled()
{
    const char * value = CachedEnv(  "ONEFLOW_HIP_GRADIENT_PACKED"  );
    if ( value != nullptr ) return value[ 0 ] == '1';
    return HipProductionBackendSelected();
}

// Fields the legacy host exchange still carries but that no consumer reads on
// the validated inviscid device path: the limiter is meaningless with ilim=0,
// the time-step gradient and heat flux belong to the viscous/heat-flux models,
// visl/vist are laminar/turbulent viscosities, and dq is the LUSGS increment.
// Skipping both their pack and their download leaves the message layout intact
// (the storage is still exchanged) and lets the limiter be filled once instead
// of once per RK stage.
inline bool HipUnusedInterfaceFieldsSkipped()
{
    const char * value = CachedEnv(  "ONEFLOW_HIP_SKIP_UNUSED_INTERFACE"  );
    if ( value != nullptr && value[ 0 ] == '0' ) return false;
    return HipProductionBackendSelected();
}

inline bool HipSkipHostFieldPackRequested()
{
    const char * value = CachedEnv(  "ONEFLOW_HIP_SKIP_HOST_FIELD_PACK"  );
    return value == nullptr || value[ 0 ] != '0';
}

// Device-authoritative primitive state.
//
// The explicit RK device chain (gradient, reconstruction, residual, state
// update, boundary, local timestep) already keeps q on the device for the whole
// stage; only the host mirror forced a 5-equation D2H per stage.  With this
// switch the device owns the newest q: the state update no longer downloads
// the field, the boundary stage no longer copies ghost q/bc_q back, and the
// host-only gamma/temperature loops are skipped.  Host readers (aerodynamic
// force, pressure coefficient, restart, visualization, diagnostics) pull back
// exactly what they need: the boundary-adjacent cells per step, the complete
// field on output steps.  The host mirror must never be written while it is
// stale; the backend rejects such an upload instead of silently restoring an
// old field.
inline bool HipDeviceAuthoritativeStateEnabled()
{
    const char * value =
        CachedEnv(  "ONEFLOW_HIP_DEVICE_AUTHORITATIVE_Q"  );
    if ( value != nullptr ) return value[ 0 ] != '0';
    return HipProductionBackendSelected();
}

// Diagnostic bisection switch: keep the per-stage ghost q D2H even when the
// device-authoritative state is on (restores the pre-round host ghost mirror
// at the cost of the transfer).
inline bool HipGhostQMirrorForced()
{
    const char * value = CachedEnv(  "ONEFLOW_HIP_GHOST_Q_MIRROR"  );
    return value != nullptr && value[ 0 ] == '1';
}

// Host RK slot copies. FIELD_FLOW saves the step's base primitive state for
// the host explicit update and for multigrid restriction. FIELD_RHS stores
// the negated residual used by the multigrid defect. The device explicit
// chain does not read either slot. A later device Advance should keep that
// base state on the device; do not bring these full-field host copies back
// as the resident implementation. ONEFLOW_HIP_SKIP_HOST_RK_SLOT_COPY=0
// restores both copies for bisection. Restart and visualization read the
// named q mirror, not FIELD_FLOW, so this skip must not replace
// NsRefreshHostStateMirror.
inline bool HipSkipHostRkSlotCopyRequested()
{
    const char * value = CachedEnv(  "ONEFLOW_HIP_SKIP_HOST_RK_SLOT_COPY"  );
    return value == nullptr || value[ 0 ] != '0';
}

inline bool HipSkipHostRkSlotCopy(
    bool skipRequested,
    bool deviceAuthoritative,
    bool deviceStateUpdate,
    bool finestGrid,
    bool singleGridLevel,
    bool explicitMultiStage,
    bool steadyStrategy )
{
    return skipRequested
        && deviceAuthoritative
        && deviceStateUpdate
        && finestGrid
        && singleGridLevel
        && explicitMultiStage
        && steadyStrategy;
}

// Measurement-only probe: skip the inter-zone interface exchange entirely so
// the per-step cost of communication can be bounded. The flow is then PHYSICALLY
// WRONG (each zone keeps stale ghosts), so it must never be used for accuracy or
// performance claims; it only answers "how much would removing communication
// save". Enabled solely by ONEFLOW_DIAG_NO_COMM=1.
// Keep the time-step / state-update status flags on the device and check them
// once per macro step instead of after every stage call. Enabled by default;
// set ONEFLOW_HIP_DEFER_STATUS=0 to restore per-call error detection.
inline bool HipStatusCheckDeferred()
{
    const char * value = CachedEnv( "ONEFLOW_HIP_DEFER_STATUS" );
    return value == nullptr || value[ 0 ] != '0';
}

// Issue the residual scale from the state-update entry instead of its own call.
// Both kernels touch the same resident arrays and no host reader sits between
// them; the default keeps the two separate entries for bisection.
inline bool HipFuseScaleUpdate()
{
    const char * value = CachedEnv( "ONEFLOW_HIP_FUSE_SCALE_UPDATE" );
    return value == nullptr || value[ 0 ] != '0';
}

inline bool DiagNoCommunication()
{
    const char * value = CachedEnv( "ONEFLOW_DIAG_NO_COMM" );
    return value != nullptr && value[ 0 ] == '1';
}

// A same-rank interface pair lives in one process and one device context, so
// its halo can travel device-to-device instead of through host staging plus an
// MPI self-send. Enabled by default; =0 restores the host-staged local path
// for A/B measurement.
inline bool HipLocalDeviceExchangeEnabled()
{
    const char * value = CachedEnv( "ONEFLOW_HIP_LOCAL_DEVICE_EXCHANGE" );
    return value == nullptr || value[ 0 ] != '0';
}

// Host history levels (q1/q2 and their residual counterparts) exist for the
// dual-time scheme and for the restart stream. With the device-authoritative
// explicit chain neither reads them during stepping: dual time is excluded by
// the capability check, and the restart dump rebuilds them from the resident
// field. The per-step shift copies the whole field once per level.
// ONEFLOW_HIP_SKIP_HOST_HISTORY_SHIFT=0 restores it for bisection.
inline bool HipSkipHostHistoryShiftRequested()
{
    const char * value =
        CachedEnv(  "ONEFLOW_HIP_SKIP_HOST_HISTORY_SHIFT"  );
    return value == nullptr || value[ 0 ] != '0';
}

// Diagnostic bisection switch: skip the per-step boundary-adjacent q mirror
// the aerodynamic force and pressure coefficient read.
inline bool HipBoundaryStateMirrorEnabled()
{
    const char * value =
        CachedEnv(  "ONEFLOW_HIP_BOUNDARY_STATE_MIRROR"  );
    return value == nullptr || value[ 0 ] != '0';
}

// The device gradient interface path owns dqdx/dqdy/dqdz on both sides: the
// send side is packed once per exchange round, and the receive side is
// unpacked on the device from a buffer uploaded once per round.  Outside this
// configuration the host pack/download keeps running so viscous or non-packed
// runs still see the host mirror.
inline bool HipDeviceGradientInterfaceEnabled( bool viscousModelActive )
{
    if ( viscousModelActive ) return false;
    const char * unpack = CachedEnv(  "ONEFLOW_HIP_DEVICE_GRADIENT_UNPACK"  );
    if ( unpack != nullptr && unpack[ 0 ] == '0' ) return false;
    return HipPackedGradientExchangeEnabled()
        && HipSkipHostFieldPackRequested()
        && HipProductionBackendSelected();
}

// One line whenever a device hook is asked for but cannot run.
//
// A decline is not an error - the host loop still computes the right answer -
// but it is indistinguishable from success in the output otherwise, and a run
// that looks accelerated while running entirely on the host is worse than no
// acceleration at all: every measurement taken on it is wrong.  Printing costs
// one line per declined call, and the message names the hook and the reason so
// a log answers "did the device path run" without a debugger.
template < typename Reason >
inline void ReportHipHookDeclined( const char * hook, const Reason & reason )
{
    std::cout << "HIP_HOOK_DECLINED " << hook << " reason=" << reason
              << std::endl;
}

EndNameSpace
