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

#pragma once

#include "NamespaceMacros.h"

#include <cstdint>

BeginNameSpace( ONEFLOW )

// Work the next device consumer must do before it may read the resident
// primitive state.
enum class ResidentStateSync
{
    None,
    GhostCells,
    FullField
};

// Explicit ownership/generation contract for a backend-resident primitive
// state.
//
// The 3D main solver keeps a host mirror of the conserved/primitive field for
// the CPU oracle and for the host boundary/interface path.  "The device copy
// is current" is therefore a property that has to be tracked, not assumed.
// This class is the single place that answers three questions:
//
//   1. Does the resident copy still belong to this solver/zone/grid state
//      (owner token, topology generation, extents, equation count)?
//   2. Is the resident copy still the same field generation as the host
//      mirror (owned-cell updates)?
//   3. Were host ghost cells rewritten in a way the device cannot observe
//      (physical boundary, interface exchange)?
//
// Ghost mutations are tracked with an epoch: the host increments
// `hostGhostEpoch` when it rewrites ghost cells, the device records the epoch
// it currently reflects in `deviceGhostEpoch`.  A ghost mismatch only requires
// an O(nGhostCells) refresh instead of a full field upload.
//
// The contract also tracks *which side owns the newest field*.  A device RK
// stage that does not download its result leaves the device current and the
// host mirror stale; `HostMirrorCurrent()` is false from then until a host
// reader pulls the field back with `NoteHostMirrorRefresh()`.  Host code must
// not write the mirror while it is stale: an upload of stale host data would
// silently overwrite the device state.
class ResidentStateContract
{
public:
    void Invalidate() noexcept
    {
        ownerToken_ = nullptr;
        topologyGeneration_ = 0;
        fieldGeneration_ = 0;
        hostGhostEpoch_ = 0;
        deviceGhostEpoch_ = 0;
        nCells_ = 0;
        nGhostCells_ = 0;
        nEquations_ = 0;
        bound_ = false;
        resident_ = false;
        hostMirrorCurrent_ = false;
    }

    // Bind the contract to one state owner.  A different owner, topology
    // generation or extent always restarts from an empty resident state.
    void Bind( const void * ownerToken,
               std::uint64_t topologyGeneration,
               int nCells,
               int nGhostCells,
               int nEquations ) noexcept
    {
        if ( ownerToken == nullptr || topologyGeneration == 0 || nCells <= 0
             || nGhostCells < 0 || nEquations <= 0 )
        {
            Invalidate();
            return;
        }
        const bool sameBinding =
            bound_
            && ownerToken_ == ownerToken
            && topologyGeneration_ == topologyGeneration
            && nCells_ == nCells
            && nGhostCells_ == nGhostCells
            && nEquations_ == nEquations;
        if ( ! sameBinding )
        {
            fieldGeneration_ = 0;
            hostGhostEpoch_ = 0;
            deviceGhostEpoch_ = 0;
            resident_ = false;
            hostMirrorCurrent_ = false;
        }
        ownerToken_ = ownerToken;
        topologyGeneration_ = topologyGeneration;
        nCells_ = nCells;
        nGhostCells_ = nGhostCells;
        nEquations_ = nEquations;
        bound_ = true;
    }

    bool Bound() const noexcept { return bound_; }
    bool Resident() const noexcept { return resident_; }

    // True when the host mirror holds the same field as the resident device
    // state.  False after a device-only advance (device authoritative), after
    // an invalidation, or before the first upload.
    bool HostMirrorCurrent() const noexcept
    {
        return resident_ && hostMirrorCurrent_;
    }

    bool Matches( const void * ownerToken,
                  std::uint64_t topologyGeneration,
                  int nCells,
                  int nGhostCells,
                  int nEquations ) const noexcept
    {
        return bound_
            && ownerToken_ == ownerToken
            && topologyGeneration_ == topologyGeneration
            && nCells_ == nCells
            && nGhostCells_ == nGhostCells
            && nEquations_ == nEquations;
    }

    const void * OwnerToken() const noexcept { return ownerToken_; }
    std::uint64_t TopologyGeneration() const noexcept
    {
        return topologyGeneration_;
    }
    std::uint64_t FieldGeneration() const noexcept { return fieldGeneration_; }
    std::uint64_t HostGhostEpoch() const noexcept { return hostGhostEpoch_; }
    std::uint64_t DeviceGhostEpoch() const noexcept
    {
        return deviceGhostEpoch_;
    }
    int Cells() const noexcept { return nCells_; }
    int GhostCells() const noexcept { return nGhostCells_; }
    int Equations() const noexcept { return nEquations_; }

    bool HostGhostPending() const noexcept
    {
        return resident_ && hostGhostEpoch_ != deviceGhostEpoch_;
    }

    // Decide what has to happen before the next device read of `q`.
    ResidentStateSync RequiredSync(
        std::uint64_t hostFieldGeneration ) const noexcept
    {
        if ( ! resident_ ) return ResidentStateSync::FullField;
        if ( hostFieldGeneration != fieldGeneration_ )
            return ResidentStateSync::FullField;
        if ( hostGhostEpoch_ != deviceGhostEpoch_ )
            return ResidentStateSync::GhostCells;
        return ResidentStateSync::None;
    }

    // The device received the complete host field (initialization, restart,
    // CPU oracle hand-off, or a host-owned state update).
    void NoteFullFieldSync( std::uint64_t hostFieldGeneration ) noexcept
    {
        fieldGeneration_ = hostFieldGeneration;
        deviceGhostEpoch_ = hostGhostEpoch_;
        resident_ = true;
        hostMirrorCurrent_ = true;
    }

    // The device received only the ghost region.  The owned cells are
    // unchanged, so the mirror keeps whatever currency it already had.
    void NoteGhostSync() noexcept
    {
        deviceGhostEpoch_ = hostGhostEpoch_;
        resident_ = true;
    }

    // The host rewrote ghost cells (host boundary model or host interface
    // exchange).  The owned cells are untouched.
    void NoteHostGhostMutation() noexcept
    {
        ++hostGhostEpoch_;
    }

    // A device kernel advanced the owned cells and the host mirror was
    // refreshed from the device.  Ghost cells are unchanged.
    void NoteDeviceFieldUpdate( std::uint64_t hostFieldGeneration ) noexcept
    {
        fieldGeneration_ = hostFieldGeneration;
        resident_ = true;
        hostMirrorCurrent_ = true;
    }

    // A device kernel writes ghost cells itself (interface unpack or device
    // boundary model).  `hostGhostsMirrored` says whether the host mirror was
    // updated with the same values; when it was not, the mirror is stale.
    void NoteDeviceGhostWrite( bool hostGhostsMirrored ) noexcept
    {
        deviceGhostEpoch_ = hostGhostEpoch_;
        resident_ = true;
        if ( ! hostGhostsMirrored ) hostMirrorCurrent_ = false;
    }

    void NoteDeviceGhostUpdate() noexcept
    {
        NoteDeviceGhostWrite( true );
    }

    // The device holds cells the host mirror does not reflect (for example a
    // ghost write that was not copied back).
    void NoteHostMirrorStale() noexcept
    {
        if ( resident_ ) hostMirrorCurrent_ = false;
    }

    // A device kernel advanced the owned cells and the result was *not*
    // downloaded: the device is now authoritative and the host mirror is
    // stale until a host reader refreshes it.
    void NoteDeviceFieldAdvance( std::uint64_t hostFieldGeneration ) noexcept
    {
        fieldGeneration_ = hostFieldGeneration;
        resident_ = true;
        hostMirrorCurrent_ = false;
    }

    // A host reader pulled the complete resident field back into the host
    // mirror (owned + ghost cells).  Both sides agree again.
    void NoteHostMirrorRefresh( std::uint64_t hostFieldGeneration ) noexcept
    {
        fieldGeneration_ = hostFieldGeneration;
        deviceGhostEpoch_ = hostGhostEpoch_;
        resident_ = true;
        hostMirrorCurrent_ = true;
    }

private:
    const void * ownerToken_ = nullptr;
    std::uint64_t topologyGeneration_ = 0;
    std::uint64_t fieldGeneration_ = 0;
    std::uint64_t hostGhostEpoch_ = 0;
    std::uint64_t deviceGhostEpoch_ = 0;
    int nCells_ = 0;
    int nGhostCells_ = 0;
    int nEquations_ = 0;
    bool bound_ = false;
    bool resident_ = false;
    bool hostMirrorCurrent_ = false;
};

EndNameSpace
