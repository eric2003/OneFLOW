/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\*---------------------------------------------------------------------------*/

// Host-only contract test for the resident per-zone device state generation /
// ownership rules shared by the HIP production path.  It runs without a DCU:
// the generation bookkeeping, not the kernel numerics, is under test.

#include "ResidentStateContract.h"

#include <gtest/gtest.h>

#include <cstdint>

using namespace ONEFLOW;

namespace
{

constexpr int nCells = 100;
constexpr int nGhostCells = 12;
constexpr int nEquations = 5;

int ownerA = 0;
int ownerB = 0;

ResidentStateContract MakeBound( std::uint64_t topologyGeneration = 1,
                                 int cells = nCells,
                                 int ghosts = nGhostCells )
{
    ResidentStateContract contract;
    contract.Bind( & ownerA, topologyGeneration, cells, ghosts, nEquations );
    return contract;
}

}

TEST( ResidentStateContract, UnboundStateRequiresFullUpload )
{
    ResidentStateContract contract;
    EXPECT_FALSE( contract.Bound() );
    EXPECT_FALSE( contract.Resident() );
    EXPECT_EQ( contract.RequiredSync( 1 ),
               ResidentStateSync::FullField );
    EXPECT_EQ( contract.Cells(), 0 );
}

TEST( ResidentStateContract, BindingRejectsInvalidIdentity )
{
    ResidentStateContract contract;
    contract.Bind( nullptr, 1, nCells, nGhostCells, nEquations );
    EXPECT_FALSE( contract.Bound() );
    contract.Bind( & ownerA, 0, nCells, nGhostCells, nEquations );
    EXPECT_FALSE( contract.Bound() );
    contract.Bind( & ownerA, 1, 0, nGhostCells, nEquations );
    EXPECT_FALSE( contract.Bound() );
    contract.Bind( & ownerA, 1, nCells, -1, nEquations );
    EXPECT_FALSE( contract.Bound() );
    contract.Bind( & ownerA, 1, nCells, nGhostCells, 0 );
    EXPECT_FALSE( contract.Bound() );
}

TEST( ResidentStateContract, FullSyncMakesResidentStateCurrent )
{
    ResidentStateContract contract = MakeBound();
    EXPECT_EQ( contract.RequiredSync( 7 ), ResidentStateSync::FullField );

    contract.NoteFullFieldSync( 7 );
    EXPECT_TRUE( contract.Resident() );
    EXPECT_EQ( contract.FieldGeneration(), 7u );
    EXPECT_EQ( contract.RequiredSync( 7 ), ResidentStateSync::None );
    EXPECT_FALSE( contract.HostGhostPending() );
}

TEST( ResidentStateContract, HostFieldMutationInvalidatesOwnedCells )
{
    ResidentStateContract contract = MakeBound();
    contract.NoteFullFieldSync( 7 );
    // A host RK update bumps the owner field generation.
    EXPECT_EQ( contract.RequiredSync( 8 ), ResidentStateSync::FullField );
}

TEST( ResidentStateContract, HostGhostMutationOnlyNeedsGhostRefresh )
{
    ResidentStateContract contract = MakeBound();
    contract.NoteFullFieldSync( 7 );
    contract.NoteHostGhostMutation();
    EXPECT_TRUE( contract.HostGhostPending() );
    EXPECT_EQ( contract.RequiredSync( 7 ), ResidentStateSync::GhostCells );

    contract.NoteGhostSync();
    EXPECT_FALSE( contract.HostGhostPending() );
    EXPECT_EQ( contract.RequiredSync( 7 ), ResidentStateSync::None );
}

TEST( ResidentStateContract, DeviceGhostWriteClearsHostGhostEpoch )
{
    ResidentStateContract contract = MakeBound();
    contract.NoteFullFieldSync( 3 );
    contract.NoteHostGhostMutation();
    EXPECT_TRUE( contract.HostGhostPending() );

    // Interface unpack or the device boundary kernel supplied the ghost cells.
    contract.NoteDeviceGhostUpdate();
    EXPECT_FALSE( contract.HostGhostPending() );
    EXPECT_EQ( contract.RequiredSync( 3 ), ResidentStateSync::None );
}

TEST( ResidentStateContract, DeviceFieldUpdateTracksNewGeneration )
{
    ResidentStateContract contract = MakeBound();
    contract.NoteFullFieldSync( 3 );
    // Device state update writes both device q and the host mirror.
    contract.NoteDeviceFieldUpdate( 4 );
    EXPECT_EQ( contract.FieldGeneration(), 4u );
    EXPECT_EQ( contract.RequiredSync( 4 ), ResidentStateSync::None );
    EXPECT_EQ( contract.RequiredSync( 3 ), ResidentStateSync::FullField );
    // Owned-cell generation change must not silently imply ghost staleness.
    EXPECT_FALSE( contract.HostGhostPending() );
}

TEST( ResidentStateContract, RebindToNewTopologyDropsResidentState )
{
    ResidentStateContract contract = MakeBound();
    contract.NoteFullFieldSync( 3 );
    contract.Bind( & ownerA, 2, nCells, nGhostCells, nEquations );
    EXPECT_FALSE( contract.Resident() );
    EXPECT_EQ( contract.RequiredSync( 3 ), ResidentStateSync::FullField );

    contract.NoteFullFieldSync( 3 );
    EXPECT_TRUE( contract.Resident() );
    // A different owner token never inherits the previous resident state.
    contract.Bind( & ownerB, 2, nCells, nGhostCells, nEquations );
    EXPECT_FALSE( contract.Resident() );
}

TEST( ResidentStateContract, RebindToSameOwnerKeepsResidentState )
{
    ResidentStateContract contract = MakeBound();
    contract.NoteFullFieldSync( 5 );
    contract.NoteHostGhostMutation();
    contract.Bind( & ownerA, 1, nCells, nGhostCells, nEquations );
    EXPECT_TRUE( contract.Resident() );
    EXPECT_TRUE( contract.HostGhostPending() );
    EXPECT_EQ( contract.RequiredSync( 5 ), ResidentStateSync::GhostCells );
}

TEST( ResidentStateContract, ExtentChangeForcesFullUpload )
{
    ResidentStateContract contract = MakeBound();
    contract.NoteFullFieldSync( 5 );
    EXPECT_TRUE( contract.Matches(
        & ownerA, 1, nCells, nGhostCells, nEquations ) );
    EXPECT_FALSE( contract.Matches(
        & ownerA, 1, nCells + 1, nGhostCells, nEquations ) );
    EXPECT_FALSE( contract.Matches(
        & ownerA, 1, nCells, nGhostCells, 3 ) );
    EXPECT_FALSE( contract.Matches(
        & ownerA, 2, nCells, nGhostCells, nEquations ) );

    contract.Bind( & ownerA, 1, nCells + 1, nGhostCells, nEquations );
    EXPECT_FALSE( contract.Resident() );
    EXPECT_EQ( contract.RequiredSync( 5 ), ResidentStateSync::FullField );
}

TEST( ResidentStateContract, ExtraHostGhostMutationReopensGhostSync )
{
    ResidentStateContract contract = MakeBound();
    contract.NoteFullFieldSync( 1 );
    contract.NoteHostGhostMutation();
    contract.NoteDeviceGhostUpdate();
    EXPECT_FALSE( contract.HostGhostPending() );

    // A later host boundary update must not be mistaken for the already
    // synchronized interface round.
    contract.NoteHostGhostMutation();
    EXPECT_TRUE( contract.HostGhostPending() );
    EXPECT_EQ( contract.RequiredSync( 1 ), ResidentStateSync::GhostCells );
}

TEST( ResidentStateContract, InvalidateDropsEverything )
{
    ResidentStateContract contract = MakeBound();
    contract.NoteFullFieldSync( 9 );
    contract.NoteHostGhostMutation();
    contract.Invalidate();
    EXPECT_FALSE( contract.Bound() );
    EXPECT_FALSE( contract.Resident() );
    EXPECT_EQ( contract.HostGhostEpoch(), 0u );
    EXPECT_EQ( contract.DeviceGhostEpoch(), 0u );
    EXPECT_FALSE( contract.HostMirrorCurrent() );
    EXPECT_EQ( contract.RequiredSync( 9 ), ResidentStateSync::FullField );
}

TEST( ResidentStateContract, DeviceAdvanceMakesTheMirrorStale )
{
    ResidentStateContract contract = MakeBound();
    contract.NoteFullFieldSync( 1 );
    EXPECT_TRUE( contract.HostMirrorCurrent() );

    // A device RK stage advanced the owned cells without a download.
    contract.NoteDeviceFieldAdvance( 2 );
    EXPECT_TRUE( contract.Resident() );
    EXPECT_FALSE( contract.HostMirrorCurrent() );
    // The device is current, so no upload is requested.
    EXPECT_EQ( contract.RequiredSync( 2 ), ResidentStateSync::None );
}

TEST( ResidentStateContract, HostMirrorRefreshRestoresBothSides )
{
    ResidentStateContract contract = MakeBound();
    contract.NoteFullFieldSync( 1 );
    contract.NoteDeviceFieldAdvance( 2 );
    contract.NoteDeviceGhostUpdate();

    contract.NoteHostMirrorRefresh( 2 );
    EXPECT_TRUE( contract.HostMirrorCurrent() );
    EXPECT_FALSE( contract.HostGhostPending() );
    EXPECT_EQ( contract.RequiredSync( 2 ), ResidentStateSync::None );
}

TEST( ResidentStateContract, GhostSyncKeepsStaleOwnedCellsStale )
{
    ResidentStateContract contract = MakeBound();
    contract.NoteFullFieldSync( 1 );
    contract.NoteDeviceFieldAdvance( 2 );
    contract.NoteHostGhostMutation();
    contract.NoteGhostSync();
    // Uploading ghosts does not make the old owned-cell mirror current.
    EXPECT_FALSE( contract.HostMirrorCurrent() );
}

TEST( ResidentStateContract, DeviceFieldUpdateKeepsMirrorCurrent )
{
    ResidentStateContract contract = MakeBound();
    contract.NoteFullFieldSync( 1 );
    contract.NoteDeviceFieldUpdate( 2 );
    EXPECT_TRUE( contract.HostMirrorCurrent() );
    EXPECT_EQ( contract.RequiredSync( 2 ), ResidentStateSync::None );
}
