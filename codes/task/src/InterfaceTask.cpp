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

#include "InterfaceTask.h"
#include "InterfaceEndpoint.h"
#include "Parallel.h"
#include "Zone.h"
#include "ZoneState.h"
#include "GridState.h"
#include "PIO.h"
#include "ActionState.h"
#include "DataBook.h"
#include "InterFace.h"
#include "Grid.h"
#include "UnsGrid.h"
#include "FaceTopo.h"
#include "DataBase.h"
#include "SolverDef.h"
#include "FieldWrap.h"
#include "SimuContext.h"
#include "StageProfiler.h"
#include "AccelRuntime.h"
#include "HipProductionPolicy.h"
#ifdef ONEFLOW_ENABLE_HIP
#include "HipFluxBackend.h"
#endif
#include <nlohmann/json.hpp>
#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <array>
#include <map>

BeginNameSpace( ONEFLOW )

namespace
{

bool HipInterfaceRequested()
{
    return HipProductionFeatureEnabled(
        "ONEFLOW_ENABLE_UNS_HIP_INTERFACE" );
}

bool GpuAwareMpiRequested()
{
    const char * enabled =
        std::getenv( "ONEFLOW_ENABLE_UNS_GPU_AWARE_MPI" );
    return enabled != nullptr && enabled[ 0 ] == '1';
}

bool GpuAwareMpiConfirmed()
{
    const char * confirmed =
        std::getenv( "ONEFLOW_GPU_AWARE_MPI_CONFIRMED" );
    return confirmed != nullptr && confirmed[ 0 ] == '1';
}

bool ExportFixtureMetadataRequested()
{
    const char * directory =
        std::getenv( "ONEFLOW_MZ_FIXTURE_EXPORT_DIR" );
    return directory != nullptr && directory[ 0 ] != '\0';
}

void ExportFixtureMetadataShard()
{
    static bool exported = false;
    if ( exported || ! ExportFixtureMetadataRequested() ) return;

    const char * directory =
        std::getenv( "ONEFLOW_MZ_FIXTURE_EXPORT_DIR" );
    const int nEquations = [] {
        const char * value = std::getenv( "ONEFLOW_MZ_FIXTURE_EQUATIONS" );
        if ( value == nullptr || value[ 0 ] == '\0' ) return 5;
        char * end = nullptr;
        const long parsed = std::strtol( value, & end, 10 );
        if ( end == value || *end != '\0' || parsed <= 0 || parsed > 1024 )
            throw std::invalid_argument(
                "ONEFLOW_MZ_FIXTURE_EQUATIONS must be in 1..1024" );
        return static_cast< int >( parsed );
    }();
    const char * memorySpaceValue =
        std::getenv( "ONEFLOW_MZ_FIXTURE_MEMORY_SPACE" );
    const std::string memorySpace =
        memorySpaceValue == nullptr || memorySpaceValue[ 0 ] == '\0'
            ? "HostStaged" : memorySpaceValue;
    if ( memorySpace != "HostStaged" && memorySpace != "DeviceDirect" )
        throw std::invalid_argument(
            "ONEFLOW_MZ_FIXTURE_MEMORY_SPACE must be HostStaged or DeviceDirect" );

    namespace fs = std::filesystem;
    const fs::path outputDirectory( directory );
    fs::create_directories( outputDirectory );
    nlohmann::ordered_json shard;
    shard[ "schema" ] = "oneflow.multizone.fixture/shard/v1";
    shard[ "rank" ] = Parallel::pid;
    shard[ "ranks" ] = Parallel::GetNProc();
    shard[ "zones" ] = nlohmann::ordered_json::array();
    shard[ "endpoints" ] = nlohmann::ordered_json::array();

    for ( int zone = 0; zone < ZoneState::nZones; ++ zone )
    {
        if ( ZoneState::pid[ zone ] != Parallel::pid ) continue;
        Grid * gridBase = Zone::GetGrid( zone, GridState::gridLevel );
        UnsGrid * grid = gridBase != nullptr && IsUnsGrid( gridBase->type )
            ? UnsGridCast( gridBase ) : nullptr;
        if ( grid == nullptr )
            throw std::runtime_error(
                "fixture metadata export could not resolve local unstructured zone" );

        shard[ "zones" ].push_back( {
            { "zone", zone },
            { "rank", Parallel::pid },
            { "owned_cells", grid->nCells },
            { "ghost_cells", grid->nBFaces }
        } );
        if ( grid->interFace == nullptr ) continue;

        InterFace & interface = *grid->interFace;
        for ( int neighborIndex = 0;
              neighborIndex < interface.nNeighbor; ++ neighborIndex )
        {
            const int neighborZone =
                interface.GetInterfacePair( neighborIndex ).nzid;
            if ( neighborZone < 0 ||
                 neighborZone >= static_cast< int >( ZoneState::pid.size() ) )
                throw std::runtime_error(
                    "fixture metadata export found an invalid neighbor zone" );
            const int neighborRank = ZoneState::pid[ neighborZone ];
            const int nInterfaces =
                interface.GetInterfacePair( neighborIndex ).nIFaces;
            const IntField & sendFaces =
                interface.GetInterfaceId( neighborIndex, GREAT_SEND );
            const IntField & receiveFaces =
                interface.GetInterfaceId( neighborIndex, GREAT_RECV );
            if ( sendFaces.size() != static_cast< HXSize_t >( nInterfaces )
                 || receiveFaces.size() != static_cast< HXSize_t >( nInterfaces ) )
                throw std::runtime_error(
                    "fixture metadata export found inconsistent interface face counts" );

            for ( int ordinal = 0; ordinal < nInterfaces; ++ ordinal )
            {
                int sendCell = -1;
                int receiveCell = -1;
                grid->GetFaceTopo().GetSId( sendFaces[ ordinal ], 1, sendCell );
                grid->GetFaceTopo().GetTId(
                    receiveFaces[ ordinal ], 1, receiveCell );
                BoundaryEndpointMap sourceEndpoint{
                    grid->nCells + grid->nBFaces,
                    { sendCell },
                    grid->nCells,
                    BoundaryEndpointRole::OwnedCell
                };
                BoundaryEndpointMap targetEndpoint{
                    grid->nCells + grid->nBFaces,
                    { receiveCell },
                    grid->nCells,
                    BoundaryEndpointRole::GhostCell
                };
                sourceEndpoint.Validate();
                targetEndpoint.Validate();
                shard[ "endpoints" ].push_back( {
                    { "zone", zone },
                    { "rank", Parallel::pid },
                    { "neighbor_zone", neighborZone },
                    { "neighbor_rank", neighborRank },
                    { "ordinal", ordinal },
                    { "direction", "send" },
                    { "face", sendFaces[ ordinal ] },
                    { "cell", sendCell },
                    { "role", "OwnedCell" },
                    { "equations", nEquations },
                    { "memory_space", memorySpace }
                } );
                shard[ "endpoints" ].push_back( {
                    { "zone", zone },
                    { "rank", Parallel::pid },
                    { "neighbor_zone", neighborZone },
                    { "neighbor_rank", neighborRank },
                    { "ordinal", ordinal },
                    { "direction", "receive" },
                    { "face", receiveFaces[ ordinal ] },
                    { "cell", receiveCell },
                    { "role", "GhostCell" },
                    { "equations", nEquations },
                    { "memory_space", memorySpace }
                } );
            }
        }
    }

    const fs::path finalPath =
        outputDirectory / ( "rank-" + std::to_string( Parallel::pid ) + ".json" );
    const fs::path temporaryPath = finalPath.string() + ".tmp";
    {
        std::ofstream output( temporaryPath );
        if ( ! output )
            throw std::runtime_error(
                "fixture metadata export could not open temporary shard" );
        output << shard.dump( 2 ) << '\n';
        if ( ! output )
            throw std::runtime_error(
                "fixture metadata export could not write shard" );
    }
    fs::rename( temporaryPath, finalPath );
    exported = true;
}

#ifdef ONEFLOW_ENABLE_HIP
bool HostInterfaceSkipped()
{
    return HipSkipHostInterfaceExchange() && ZoneState::nZones > 1;
}

namespace {

// The readiness and value-count handshakes in ExchangeHipInterfaceQ are two
// blocking round trips per interface exchange.  They validate the endpoint
// contract, which is a property of the partition, not of the time step: once a
// pair has been proven, the steady state only needs the payload message.
// Caching them removes the per-step waiting that dominated the q halo.  Set
// ONEFLOW_HIP_CACHE_INTERFACE_HANDSHAKE=0 to restore the per-step handshakes.
bool ExchangeHandshakeCacheEnabled()
{
    const char * value =
        std::getenv( "ONEFLOW_HIP_CACHE_INTERFACE_HANDSHAKE" );
    return value == nullptr || value[ 0 ] != '0';
}

std::map< std::array< int, 5 >, int > & ExchangeCounts()
{
    static std::map< std::array< int, 5 >, int > counts;
    return counts;
}

// The handshakes validate two things: the endpoints exist and both sides agree
// on the payload length.  Both are properties of the partition, so they are
// cached per (zone, zone, rank, rank, equations) pair.  The cached value is the
// *verified* length, used symmetrically by both ranks; caching each side's own
// length independently would let a mismatch slip through as silent truncation.
const int * CachedExchangeCount( const std::array< int, 5 > & key )
{
    if ( ! ExchangeHandshakeCacheEnabled() ) return nullptr;
    const auto found = ExchangeCounts().find( key );
    return found == ExchangeCounts().end() ? nullptr : & found->second;
}

void RecordExchangeCount( const std::array< int, 5 > & key, int count )
{
    if ( ! ExchangeHandshakeCacheEnabled() ) return;
    ExchangeCounts()[ key ] = count;
}

}

void ExchangeHipInterfaceQ(
    SimuContext & context,
    int sourceZone,
    int targetZone,
    int sourceRank,
    int targetRank )
{
    AccelRuntime & runtime = AccelRuntime::Instance();
    if ( ! runtime.IsInitialized()
         || runtime.Backend().Kind() != AccelBackendKind::HIP )
        return;

    const bool isSource = Parallel::pid == sourceRank;
    const bool isTarget = Parallel::pid == targetRank;
    if ( ! isSource && ! isTarget ) return;

    const int previousZone = ZoneState::zid;
    const int tagBase =
        700000 + sourceZone * ZoneState::nZones + targetZone;
    UnsGrid * sourceGrid = nullptr;
    UnsGrid * targetGrid = nullptr;
    BoundaryEndpointMap sourceEndpoint;
    BoundaryEndpointMap targetEndpoint;
    HipFluxBackend & backend = HipFluxBackend::Shared();
    const void * sourceKey = nullptr;
    const void * targetKey = nullptr;
    MRField * sourceQ = nullptr;
    MRField * targetQ = nullptr;
    int nEquations = 5;
    const Real * sourceComponents[ 5 ] = {};
    const Real * targetComponents[ 5 ] = {};
    const bool useDeviceDirect =
        sourceRank != targetRank
        && GpuAwareMpiRequested()
        && GpuAwareMpiConfirmed();

    if ( sourceRank != targetRank
         && GpuAwareMpiRequested()
         && ! GpuAwareMpiConfirmed() )
    {
        throw std::runtime_error(
            "GPU-aware MPI was requested without a confirmed device-buffer "
            "MPI probe; set ONEFLOW_GPU_AWARE_MPI_CONFIRMED=1 only after "
            "the target MPI has passed that probe" );
    }

    // Establish readiness symmetrically before any count or payload
    // communication.  A missing/incompatible local endpoint must not leave
    // the peer blocked in HXSwapData.
    int sourceReady = 0;
    int targetReady = 0;
    if ( isSource )
    {
        try
        {
            Grid * gridBase =
                Zone::GetGrid( sourceZone, GridState::gridLevel );
            sourceGrid = gridBase != nullptr && IsUnsGrid( gridBase->type )
                ? UnsGridCast( gridBase ) : nullptr;
            sourceQ = sourceGrid == nullptr
                ? nullptr : GetFieldPointer< MRField >( sourceGrid, "q" );
            sourceReady =
                sourceGrid != nullptr
                && sourceQ != nullptr
                && sourceQ->GetNEqu() == nEquations
                && BuildInterfaceEndpoint(
                    sourceGrid, targetZone, GREAT_SEND,
                    BoundaryEndpointRole::OwnedCell, sourceEndpoint );
        }
        catch ( const std::exception & )
        {
            sourceReady = 0;
        }
    }

    if ( isTarget )
    {
        try
        {
            Grid * gridBase =
                Zone::GetGrid( targetZone, GridState::gridLevel );
            targetGrid = gridBase != nullptr && IsUnsGrid( gridBase->type )
                ? UnsGridCast( gridBase ) : nullptr;
            targetQ = targetGrid == nullptr
                ? nullptr : GetFieldPointer< MRField >( targetGrid, "q" );
            targetReady =
                targetGrid != nullptr
                && targetQ != nullptr
                && targetQ->GetNEqu() == nEquations
                && BuildInterfaceEndpoint(
                    targetGrid, sourceZone, GREAT_RECV,
                    BoundaryEndpointRole::GhostCell, targetEndpoint );
        }
        catch ( const std::exception & )
        {
            targetReady = 0;
        }
    }

    const std::array< int, 5 > exchangeKey{
        sourceZone, targetZone, sourceRank, targetRank, nEquations };
    const int * cachedCount = CachedExchangeCount( exchangeKey );
    const bool validated = cachedCount != nullptr;
    if ( ! validated )
    {
        if ( ! ( isSource && isTarget ) )
        {
            ScopedStageTimer readyTimer( "interface_ready" );
            HXSwapData(
                & sourceReady, 1, sourceRank, targetRank, tagBase );
            HXSwapData(
                & targetReady, 1, targetRank, sourceRank, tagBase + 1 );
        }
        if ( sourceReady == 0 || targetReady == 0 )
        {
            ZoneState::zid = previousZone;
            return;
        }
    }
    else
    {
        // Cached: the endpoint contract was already proven for this pair.
        if ( ( isSource && sourceReady == 0 )
             || ( isTarget && targetReady == 0 ) )
        {
            ZoneState::zid = previousZone;
            return;
        }
    }

    if ( isSource )
    {
        sourceKey = sourceGrid;
        for ( int equation = 0; equation < nEquations; ++ equation )
            sourceComponents[ equation ] =
                ( *sourceQ )[ equation ].data();
    }
    if ( isTarget )
    {
        targetKey = targetGrid;
        for ( int equation = 0; equation < nEquations; ++ equation )
            targetComponents[ equation ] =
                ( *targetQ )[ equation ].data();
    }

    BoundaryBufferView send;
    BoundaryBufferView receive;
    const HipBoundaryEndpointKey sourceEndpointKey{
        sourceZone, targetZone, 0 };
    const HipBoundaryEndpointKey targetEndpointKey{
        targetZone, sourceZone, 0 };
    if ( isSource )
    {
        ZoneState::zid = sourceZone;
        {
            ScopedStageTimer uploadTimer( "interface_upload" );
            backend.UploadBoundaryState(
                sourceKey, sourceGrid->nCells, sourceGrid->nBFaces,
                nEquations, sourceComponents );
        }
        ScopedStageTimer packTimer( "interface_pack" );
        send = backend.PackBoundary(
            sourceKey, sourceEndpoint, nEquations,
            useDeviceDirect
                ? BoundaryTransferMode::DeviceDirect
                : BoundaryTransferMode::HostStaged,
            sourceEndpointKey );
    }
    if ( isTarget )
    {
        ZoneState::zid = targetZone;
        {
            ScopedStageTimer uploadTimer( "interface_upload" );
            backend.UploadBoundaryState(
                targetKey, targetGrid->nCells, targetGrid->nBFaces,
                nEquations, targetComponents );
        }
        ScopedStageTimer packTimer( "interface_pack" );
        receive = backend.PrepareBoundaryReceive(
            targetKey, targetEndpoint, nEquations,
            useDeviceDirect
                ? BoundaryTransferMode::DeviceDirect
                : BoundaryTransferMode::HostStaged,
            targetEndpointKey );
    }

    const BoundaryTransferMode transferMode = useDeviceDirect
        ? BoundaryTransferMode::DeviceDirect
        : BoundaryTransferMode::HostStaged;
    if ( isSource )
    {
        ValidateBoundaryBuffer(
            sourceEndpoint, nEquations, send, transferMode );
    }
    if ( isTarget )
    {
        ValidateBoundaryBuffer(
            targetEndpoint, nEquations, receive, transferMode );
    }
    if ( isSource && isTarget )
    {
        BoundaryExchangeContract contract;
        contract.sourceZone = sourceZone;
        contract.targetZone = targetZone;
        contract.sourceRank = sourceRank;
        contract.targetRank = targetRank;
        contract.nEquations = nEquations;
        contract.source = sourceEndpoint;
        contract.target = targetEndpoint;
        ValidateBoundaryExchangeBuffers(
            contract, send, receive, transferMode );
    }

    {
        ScopedStageTimer mpiTimer( "interface_mpi" );
        int sendCount = isSource
            ? static_cast< int >( send.valueCount ) : 0;
        int payloadCount = 0;
        if ( ! validated )
        {
            HXSwapData( &sendCount, 1, sourceRank, targetRank, tagBase );
            if ( isTarget )
            {
                const int receiveCount =
                    static_cast< int >( receive.valueCount );
                if ( sendCount != receiveCount )
                    throw std::runtime_error(
                        "HIP interface exchange value count mismatch" );
            }
            // The count exchange is itself the validation: after it, both
            // ranks hold the same verified length, so caching each side's own
            // copy is symmetric without a second message.
            payloadCount = sendCount;
            RecordExchangeCount( exchangeKey, payloadCount );
        }
        else
        {
            payloadCount = * cachedCount;
            if ( isSource
                 && static_cast< int >( send.valueCount ) != payloadCount )
            {
                throw std::runtime_error(
                    "HIP interface exchange send length differs from the "
                    "validated contract" );
            }
            if ( isTarget
                 && static_cast< int >( receive.valueCount ) != payloadCount )
            {
                throw std::runtime_error(
                    "HIP interface receive length differs from the validated "
                    "contract" );
            }
        }

        if ( isSource && isTarget )
        {
            // Same-rank exchanges remain host-staged.  Device pointers are not
            // valid inputs to std::copy, and the host path is deterministic for
            // the local multi-zone case.
            if ( send.valueCount != receive.valueCount )
            {
                throw std::runtime_error(
                    "local HIP interface exchange value count mismatch" );
            }
            std::copy(
                send.values,
                send.values + send.valueCount,
                receive.values );
        }
        else
        {
            if ( useDeviceDirect )
            {
                Real * data = isSource
                    ? send.values : receive.values;
                HXSwapDeviceData(
                    data,
                    static_cast< std::size_t >( payloadCount ) * sizeof( Real ),
                    sourceRank,
                    targetRank,
                    tagBase + 1 );
            }
            else
            {
                Real * data = isSource
                    ? send.values : receive.values;
                HXSwapData(
                    data, payloadCount, sourceRank, targetRank, tagBase + 1 );
            }
        }
    }

    if ( isTarget )
    {
        ScopedStageTimer unpackTimer( "interface_unpack" );
        backend.UnpackBoundary(
            targetKey, targetEndpoint, nEquations,
            useDeviceDirect
                ? BoundaryTransferMode::DeviceDirect
                : BoundaryTransferMode::HostStaged,
            targetEndpointKey );
        // The device ghosts are newer than the host mirror unless the paired
        // host exchange rewrote the same cells.  Publish that so the next host
        // reader pulls the field back instead of trusting a stale mirror.
        if ( HostInterfaceSkipped() )
        {
            backend.MarkResidentHostMirrorStale( targetKey );
        }
    }
    ZoneState::zid = previousZone;
    (void) context;
}

// The host interface exchange is the authoritative host-staged halo path.  A
// target zone whose ghosts were rewritten there must publish that mutation to
// the resident-state contract, because the paired HIP unpack is the only
// device-side writer for those cells.
void MarkLocalHostGhostsModified( int zone )
{
    Grid * gridBase = Zone::GetGrid( zone, GridState::gridLevel );
    if ( gridBase == nullptr ) return;
    HipFluxBackend::Shared().NotifyHostGhostStateModified( gridBase );
}

// Safety net: any ghost mutation the HIP interface path did not cover (for
// example an exchange skipped by a readiness failure) is repaired here with a
// single ghost-region upload per zone.  When the device already wrote the
// exchanged cells this is a no-op.
void SyncLocalResidentGhostStates()
{
    for ( int localZoneIndex = 0;
          localZoneIndex < ZoneState::nLocal; ++ localZoneIndex )
    {
        const int zone = ZoneState::localZid[ localZoneIndex ];
        Grid * gridBase = Zone::GetGrid( zone, GridState::gridLevel );
        UnsGrid * grid = gridBase != nullptr && IsUnsGrid( gridBase->type )
            ? UnsGridCast( gridBase ) : nullptr;
        if ( grid == nullptr || grid->nBFaces <= 0 ) continue;
        MRField * q = GetFieldPointer< MRField >( grid, "q" );
        if ( q == nullptr ) continue;
        const int nEquations = static_cast< int >( q->GetNEqu() );
        if ( nEquations < 3 || nEquations > 5 ) continue;
        const Real * components[ 5 ] = {};
        for ( int equation = 0; equation < nEquations; ++ equation )
            components[ equation ] = ( * q )[ equation ].data();
        HipFluxBackend::Shared().SyncResidentGhostStateFromHost(
            grid, grid->nCells, grid->nBFaces, nEquations, components );
    }
}
#endif

}


CUpdateInterface::CUpdateInterface()
{
    ;
}

CUpdateInterface::~CUpdateInterface()
{
    ;
}

void CUpdateInterface::Run()
{
    ActionState::dataBook = this->dataBook.get();
    ExportFixtureMetadataShard();

    std::vector< int > ownerRanks(
        ZoneState::pid.begin(),
        ZoneState::pid.end() );
    std::vector< std::vector< int > > zoneNeighbors(
        ZoneState::nZones );
    for ( int zId = 0; zId < ZoneState::nZones; ++ zId )
    {
        int nNei = Zone::GetNumberOfZoneNeighbors( zId );
        for ( int iNei = 0; iNei < nNei; ++ iNei )
        {
            zoneNeighbors[ zId ].push_back(
                Zone::GetNeighborZoneId( zId, iNei ) );
        }
    }

    const DistributedDomainPlan plan = DistributedDomainPlan::Build(
        ownerRanks,
        zoneNeighbors,
        Parallel::GetNProc() );
#ifdef ONEFLOW_ENABLE_HIP
    // When the legacy host exchange is skipped the host ghost cells are stale,
    // so the resident repair pass must not upload them over the device ghosts.
    const bool trackResidentHalo =
        HipProductionBackendSelected() && ! HostInterfaceSkipped();
#endif
    for ( const BoundaryExchange & exchange : plan.BoundaryExchanges() )
    {
#ifdef ONEFLOW_ENABLE_HIP
        if ( trackResidentHalo
             && Parallel::pid == plan.OwnerRank( exchange.targetZone ) )
        {
            // The host exchange below rewrites target ghost cells.  Mark the
            // mutation per exchange so a successful device unpack clears
            // exactly that exchange; a skipped exchange stays pending and is
            // repaired by the safety-net sync at the end of the round.
            MarkLocalHostGhostsModified( exchange.targetZone );
        }
#endif
        this->SwapInterfaceData( plan, exchange );
    }
#ifdef ONEFLOW_ENABLE_HIP
    if ( trackResidentHalo )
    {
        ScopedStageTimer ghostTimer( "interface_ghost_sync" );
        SyncLocalResidentGhostStates();
    }
#endif
}

void CUpdateInterface::SwapInterfaceData(
    const DistributedDomainPlan & plan,
    const BoundaryExchange & exchange )
{
    const int iZone = exchange.sourceZone;
    const int jZone = exchange.targetZone;
    const int sPid = plan.OwnerRank( iZone );
    const int rPid = plan.OwnerRank( jZone );

#ifdef ONEFLOW_ENABLE_HIP
    const bool skipHostExchange = HostInterfaceSkipped();
#else
    const bool skipHostExchange = false;
#endif
    if ( ! skipHostExchange )
    {
        if ( Parallel::pid == sPid )
        {
            ZoneState::zid  = iZone;
            ZoneState::rzid = jZone;

            ScopedStageTimer sendTimer( "interface_host_send" );
            this->sendAction();
        }

        {
            ScopedStageTimer hostMpiTimer( "interface_host_mpi" );
            HXSwapData( ActionState::dataBook, sPid, rPid );
        }

        if ( Parallel::pid == rPid )
        {
            ZoneState::zid  = jZone;
            ZoneState::szid = iZone;

            ScopedStageTimer recvTimer( "interface_host_recv" );
            this->recvAction();
        }
    }

#ifdef ONEFLOW_ENABLE_HIP
    if ( this->Context() != nullptr && HipInterfaceRequested() )
    {
        ExchangeHipInterfaceQ(
            *this->Context(), iZone, jZone, sPid, rPid );
    }
#endif
}

EndNameSpace
