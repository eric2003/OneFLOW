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

BeginNameSpace( ONEFLOW )

// Production HIP follows the runtime backend selection. Legacy per-stage
// variables remain overrides: unset selects automatically, "0" disables the
// stage, and "1" requests it subject to capability checks.
inline bool HipProductionBackendSelected()
{
#ifdef ONEFLOW_ENABLE_HIP
    AccelRuntime & runtime = AccelRuntime::Instance();
    return runtime.IsInitialized()
        && runtime.Backend().Kind() == AccelBackendKind::HIP;
#else
    return false;
#endif
}

inline bool HipProductionFeatureExplicitlyEnabled( const char * name )
{
    const char * value = std::getenv( name );
    return value != nullptr && value[ 0 ] == '1';
}

inline bool HipProductionFeatureEnabled( const char * legacyName )
{
    const char * value = std::getenv( legacyName );
    if ( value != nullptr ) return value[ 0 ] == '1';
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
    const char * value = std::getenv( "ONEFLOW_HIP_SKIP_HOST_INTERFACE" );
    if ( value != nullptr ) return value[ 0 ] != '0';
    return HipProductionBackendSelected();
}

// The host mirror of the boundary-adjacent gradient only feeds aerodynamic
// forces, the optional boundary-gradient store and field dumps.  Mirroring it
// every RK stage costs a device pack plus one H2D of the sampled cells;
// request it only on the steps whose host consumers actually run.
inline bool HipGradientHostMirrorOnDemand()
{
    const char * value = std::getenv( "ONEFLOW_HIP_GRADIENT_MIRROR_ON_DEMAND" );
    return value == nullptr || value[ 0 ] != '0';
}

inline bool HipPackedGradientExchangeEnabled()
{
    const char * value = std::getenv( "ONEFLOW_HIP_GRADIENT_PACKED" );
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
    const char * value = std::getenv( "ONEFLOW_HIP_SKIP_UNUSED_INTERFACE" );
    if ( value != nullptr && value[ 0 ] == '0' ) return false;
    return HipProductionBackendSelected();
}

inline bool HipSkipHostFieldPackRequested()
{
    const char * value = std::getenv( "ONEFLOW_HIP_SKIP_HOST_FIELD_PACK" );
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
        std::getenv( "ONEFLOW_HIP_DEVICE_AUTHORITATIVE_Q" );
    if ( value != nullptr ) return value[ 0 ] != '0';
    return HipProductionBackendSelected();
}

// Diagnostic bisection switch: keep the per-stage ghost q D2H even when the
// device-authoritative state is on (restores the pre-round host ghost mirror
// at the cost of the transfer).
inline bool HipGhostQMirrorForced()
{
    const char * value = std::getenv( "ONEFLOW_HIP_GHOST_Q_MIRROR" );
    return value != nullptr && value[ 0 ] == '1';
}

// Diagnostic bisection switch: skip the per-step boundary-adjacent q mirror
// the aerodynamic force and pressure coefficient read.
inline bool HipBoundaryStateMirrorEnabled()
{
    const char * value =
        std::getenv( "ONEFLOW_HIP_BOUNDARY_STATE_MIRROR" );
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
    const char * unpack = std::getenv( "ONEFLOW_HIP_DEVICE_GRADIENT_UNPACK" );
    if ( unpack != nullptr && unpack[ 0 ] == '0' ) return false;
    return HipPackedGradientExchangeEnabled()
        && HipSkipHostFieldPackRequested()
        && HipProductionBackendSelected();
}

EndNameSpace
