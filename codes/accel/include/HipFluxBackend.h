/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\---------------------------------------------------------------------------*/

#pragma once

#include "FluxBackend.h"
#include "BoundaryExchangeContract.h"
#include "DeviceBuffer.h"
#include "ResidentStateContract.h"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <unordered_map>
#include <vector>

namespace ONEFLOW
{

struct PrimitiveFaceStateView;

// The NS 3D state is only named here through a reference, an owning factory
// return type and a pointer map, so no domain header is needed. Callers that
// construct or dereference these types include Ns3DDomain.h themselves.
class Ns3DBackendState;
class Ns3DDeviceState;

// Per-zone residual diagnostics computed from the state-owned device residual.
// The solver consumes the same quantities the host residual task used to
// derive by looping the MRField: sum of squares, max |residual| with its cell
// index, and the primitive positivity extremes.  Only these values cross the
// PCIe boundary; the residual field itself stays on the device.
struct HipResidualMetrics
{
    std::vector< Real > sumOfSquares;
    std::vector< Real > maxAbsolute;
    std::vector< int > argmaxCell;
    Real minDensity = 0.0;
    Real minPressure = 0.0;
    int nonPhysicalCells = 0;
    bool valid = false;
};

struct HipBoundaryEndpointKey
{
    int localZone = -1;
    int neighborZone = -1;
    int ghostLevel = 0;

    bool operator==( const HipBoundaryEndpointKey & other ) const noexcept
    {
        return localZone == other.localZone
            && neighborZone == other.neighborZone
            && ghostLevel == other.ghostLevel;
    }
};

struct HipBoundaryEndpointKeyHash
{
    std::size_t operator()(
        const HipBoundaryEndpointKey & key ) const noexcept;
};

struct HipBoundaryEndpointStorage
{
    DeviceBuffer< int > cellIds;
    DeviceBuffer< Real > values;
    std::vector< int > cachedCellIds;
    std::vector< Real > hostValues;
    std::uint64_t topologyGeneration = 0;
    int cellCount = 0;
    int nEquations = 0;
    BoundaryEndpointRole role = BoundaryEndpointRole::Unspecified;
};

// HIP-only storage owned either by a solver/grid Ns3DBackendState or by the
// compatibility fallback in HipFluxBackend. Keeping this type in the HIP
// header avoids leaking DeviceBuffer through backend-neutral contracts.
struct HipGradientStorage
{
    DeviceBuffer< Real > cellState;
    DeviceBuffer< Real > cellGradient;
    // Derived temperature field and its cell gradient.  Temperature is not one
    // of the five conserved variables, so it has its own device storage next to
    // the mean-flow state instead of a sixth plane in cellState/cellGradient.
    // The viscous and heat-conduction kernels read these; nothing allocates
    // them until a viscous configuration is active.
    DeviceBuffer< Real > temperature;
    DeviceBuffer< Real > temperatureGradient;
    DeviceBuffer< Real > xFaceCenter;
    DeviceBuffer< Real > yFaceCenter;
    DeviceBuffer< Real > zFaceCenter;
    DeviceBuffer< Real > xNormal;
    DeviceBuffer< Real > yNormal;
    DeviceBuffer< Real > zNormal;
    DeviceBuffer< Real > meshVelocityNormal;
    DeviceBuffer< Real > faceArea;
    DeviceBuffer< Real > xCellCenter;
    DeviceBuffer< Real > yCellCenter;
    DeviceBuffer< Real > zCellCenter;
    DeviceBuffer< Real > cellVolume;
    DeviceBuffer< int > leftCell;
    DeviceBuffer< int > rightCell;
    // Cell-to-face CSR built once per topology.  The cell-based gradient and
    // residual kernels consume this list; it removes the face-parallel
    // atomicAdd scatter from the production path.
    DeviceBuffer< int > cellFaceOffsets;
    DeviceBuffer< int > cellFaceIds;
    const int * cachedCellFaceLeft = nullptr;
    const int * cachedCellFaceRight = nullptr;
    int cachedCellFaceCells = -1;
    int cachedCellFaceFaces = -1;
    int cachedCellFaceBoundaryFaces = -1;
    bool cellFaceValid = false;
    const void * cachedOwnerToken = nullptr;
    std::uint64_t cachedTopologyGeneration = 0;
    const Real * cachedXFace = nullptr;
    const Real * cachedYFace = nullptr;
    const Real * cachedZFace = nullptr;
    const Real * cachedXNormal = nullptr;
    const Real * cachedYNormal = nullptr;
    const Real * cachedZNormal = nullptr;
    const Real * cachedFaceArea = nullptr;
    const Real * cachedXCell = nullptr;
    const Real * cachedYCell = nullptr;
    const Real * cachedZCell = nullptr;
    const Real * cachedVolume = nullptr;
    const int * cachedLeftCell = nullptr;
    const int * cachedRightCell = nullptr;
    const Real * cachedMeshVelocityNormal = nullptr;
    int cachedCells = 0;
    int cachedGhostCells = -1;
    int cachedFaces = 0;
    int cachedBoundaryFaces = -1;
    bool geometryCacheValid = false;
    DeviceBuffer< Real > qf1;
    DeviceBuffer< Real > qf2;
    // Face flux is part of the state-owned reconstructed production chain.
    // The backend-level deviceFlux_ remains a compatibility buffer for the
    // legacy host-staged path only.
    DeviceBuffer< Real > faceFlux;
    DeviceBuffer< Real > residual;
    DeviceBuffer< Real > inviscidSpectrum;
    DeviceBuffer< Real > timeStep;
    DeviceBuffer< Real > minTimeStep;
    DeviceBuffer< int > timeStepStatus;
    DeviceBuffer< int > stateUpdateStatus;
    DeviceBuffer< unsigned char > boundaryMask;
    DeviceBuffer< ReconstructionBoundaryOperation > boundaryOperation;
    DeviceBuffer< Real > bcQ;
    // Per-cell/per-equation reconstruction limiter factor (host-computed).
    DeviceBuffer< Real > limiter;
    // True while boundaryFaceState holds the newest wall state produced by the
    // device boundary stage, so reconstruction can read it directly instead of
    // mirroring it to the host and uploading it back every stage.
    bool wallStateValid = false;
    int wallStateFaces = 0;
    int wallStateEquations = 0;
    DeviceBuffer< Real > boundaryFaceState;
    // Packed [equation][ghost cell] slice used for a single D2H per stage.
    DeviceBuffer< Real > boundaryGhostPack;
    // Packed [equation][owned cell] slice for on-demand host readers.
    DeviceBuffer< Real > stateMirrorPack;
    DeviceBuffer< int > stateMirrorIds;
    DeviceBuffer< int > boundaryType;
    DeviceBuffer< Real > inflow;
    DeviceBuffer< Real > boundaryXNormal;
    DeviceBuffer< Real > boundaryYNormal;
    DeviceBuffer< Real > boundaryZNormal;
    DeviceBuffer< Real > boundaryMeshVelocityNormal;
    DeviceBuffer< int > boundaryLeftCell;
    DeviceBuffer< int > boundaryRightCell;
    // Boundary geometry/topology is immutable for a bound grid across RK
    // stages. Keep its host-view identity so UpdateBoundaryState does not
    // re-upload the same metadata on every stage. Dynamic inflow values are
    // intentionally not part of this cache.
    const void * cachedBoundaryOwnerToken = nullptr;
    std::uint64_t cachedBoundaryGeometryGeneration = 0;
    const Real * cachedBoundaryXNormal = nullptr;
    const Real * cachedBoundaryYNormal = nullptr;
    const Real * cachedBoundaryZNormal = nullptr;
    const Real * cachedBoundaryMeshVelocityNormal = nullptr;
    const int * cachedBoundaryLeftCell = nullptr;
    const int * cachedBoundaryRightCell = nullptr;
    const int * cachedBoundaryType = nullptr;
    bool boundaryGeometryCacheValid = false;
    // Inflow values are tiny but the upload is a blocking H2D, so it became the
    // drain point for every deferred kernel in the stage. Compare the values
    // (not just the pointer) and upload only when they actually change.
    Real cachedInflowValues[ 8 ] = {};
    int cachedInflowCount = 0;
    std::unordered_map<
        HipBoundaryEndpointKey,
        HipBoundaryEndpointStorage,
        HipBoundaryEndpointKeyHash > boundarySend;
    std::unordered_map<
        HipBoundaryEndpointKey,
        HipBoundaryEndpointStorage,
        HipBoundaryEndpointKeyHash > boundaryReceive;
    const ReconstructionBoundaryOperation * cachedBoundaryOperationHost =
        nullptr;
    std::size_t cachedBoundaryOperationCount = 0;
    std::uint64_t cachedBoundaryTopologyGeneration = 0;
    int cachedEquations = 0;
    int residualCells = 0;
    int residualEquations = 0;
    const Real * cachedTimeStepXNormal = nullptr;
    const Real * cachedTimeStepYNormal = nullptr;
    const Real * cachedTimeStepZNormal = nullptr;
    const Real * cachedTimeStepMeshVelocityNormal = nullptr;
    const Real * cachedTimeStepFaceArea = nullptr;
    const Real * cachedTimeStepVolume = nullptr;
    const int * cachedTimeStepLeftCell = nullptr;
    const int * cachedTimeStepRightCell = nullptr;
    const Real * cachedTimeStepHost = nullptr;
    int cachedTimeStepCells = 0;
    int cachedTimeStepGhostCells = -1;
    int cachedTimeStepFaces = 0;
    int cachedTimeStepBoundaryFaces = -1;
    bool timeStepGeometryCacheValid = false;
    bool residentTimeStepReady = false;
    // Generation of the reconstruction metadata/field view.
    std::uint64_t cachedFieldGeneration = 0;
    // Single source of truth for "is the resident q current, and for which
    // owner/generation".  Ghost-only host mutations are tracked separately so
    // an interface round trip does not force a full field upload.
    ResidentStateContract residentContract;
    bool boundaryMaskValid = false;
    bool residualValid = false;
    // True between CALC_LHS (device residual scaling) and the next flux write.
    // The state update consumes residualValid, but the diagnostics still need
    // the scaled residual, so the two flags are tracked separately.
    bool residualMetricsValid = false;
    // CALC_LHS and UPDATE_FLOWFIELD are two element-wise passes over the same
    // resident arrays with nothing host-visible in between.  With
    // ONEFLOW_HIP_FUSE_SCALE_UPDATE=1 the scale kernel is deferred and issued
    // from the update entry, so a stage pays one host entry and one drain
    // instead of two.
    bool residualScalePending = false;
    Real residualScaleCoefficient = 1.0;
    const Real * cachedStateQ[ 5 ] = {};
    // Cells whose gradient is mirrored to the host for the host-side interface
    // exchange.  Kept resident so the mirror does not re-upload ids per stage.
    DeviceBuffer< int > gradientMirrorCellIds;
    DeviceBuffer< int > gradientPackIds;
    // The unpack caller now passes a cached, stable id list for the partition;
    // re-uploading it was a blocking H2D per ghost level per stage.
    const int * cachedGradientPackIdsHost = nullptr;
    int cachedGradientPackIdCount = -1;
    DeviceBuffer< Real > gradientPackValues;
    std::vector< int > cachedGradientMirrorCellIds;
    DeviceBuffer< Real > gradientMirrorValues;
    int cachedGradientMirrorEquations = 0;
    // Residual diagnostic reduction buffers (small, reused every outer step).
    DeviceBuffer< Real > residualMetricChunkSum;
    DeviceBuffer< Real > residualMetricChunkMax;
    DeviceBuffer< int > residualMetricChunkIndex;
    DeviceBuffer< Real > residualMetricChunkMin;
    int residualMetricChunks = 0;
    int residualMetricEquations = 0;
    // Cache metadata for the bound-state host-staged primitive path. The
    // unbound direct/smoke path keeps using HipFluxBackend's fallback cache.
    const void * cachedFluxOwnerToken = nullptr;
    const Real * cachedFluxXNormal = nullptr;
    const Real * cachedFluxYNormal = nullptr;
    const Real * cachedFluxZNormal = nullptr;
    const Real * cachedFluxMeshVelocityNormal = nullptr;
    const Real * cachedFluxFaceArea = nullptr;
    const int * cachedFluxLeftCell = nullptr;
    const int * cachedFluxRightCell = nullptr;
    const unsigned char * cachedFluxBoundaryMask = nullptr;
    int cachedFluxFaces = 0;
    int cachedFluxBoundaryFaces = -1;
    bool fluxGeometryCacheValid = false;
};

std::unique_ptr< Ns3DBackendState > CreateHipNs3DBackendState(
    int deviceId );

class HipFluxBackend final : public FluxBackend
{
public:
    static HipFluxBackend & Shared();
    static void ReleaseShared();

    void CalcInvFlux(
        const FaceStateView & state,
        FaceFluxView & flux,
        int scheme ) override;

    void AddFaceFlux(
        const FaceFluxView & flux,
        const FaceConnectivityView & connectivity,
        ResidualView & residual ) override;

    // Add the flux produced by the immediately preceding CalcInvFlux call
    // without uploading the host flux array again. This is the production
    // seam used by the 3D main solver HIP batch path.
    void AddCurrentFaceFlux(
        const FaceConnectivityView & connectivity,
        ResidualView & residual );

    // Fused production path for primitive face data. The kernel converts
    // primitive states, evaluates the numerical flux, and scatters directly
    // into the device residual. hostFlux may be null in NoTrace mode.
    void CalcAndAddPrimitiveFaceFlux(
        const PrimitiveFaceStateView & state,
        const FaceConnectivityView & connectivity,
        ResidualView & residual,
        int scheme,
        FaceFluxView * hostFlux = nullptr,
        const void * cacheKey = nullptr,
        Real entropy1 = 0.0,
        Real entropy2 = 0.0,
        bool entropyFix = false );

    // Single-face evaluation of the device Roe flux, for the numeric contract
    // test that pins host and device to the same formula.  Not used by the
    // solver: the production path runs the whole face array.
    struct DeviceRoeSample
    {
        Real left[ 5 ] = {};      // rho, u, v, w, p
        Real right[ 5 ] = {};     // rho, u, v, w, p
        Real normal[ 3 ] = {};
        Real meshVelocityNormal = 0.0;
        Real gamma = 1.4;
        Real entropy1 = 0.0;
        Real entropy2 = 0.0;
        bool entropyFix = false;
    };
    void EvaluateDeviceRoe(
        const DeviceRoeSample & sample, Real flux[ 5 ] ) const;

    // Diagnostic only: copy the resident per-cell limiter factor and the
    // reconstructed face values back to the host so a CPU run and a HIP run can
    // be compared point by point.  Enabled by ONEFLOW_UNS_LIMITER_DUMP.
    void CopyLimiterStateToHost(
        const void * cacheKey, int nCells, int nFaces, int nEquations,
        Real * phi, Real * faceLeft, Real * faceRight );

    // Single-cell-side Barth limiter factor, for the limiter contract test.
    // One face side: linear increment dq = grad . (xFace - xCell), face normal,
    // and the cell-local min/max range of the limited variable.
    Real EvaluateDeviceBarth(
        Real dqdx, Real dqdy, Real dqdz,
        Real dx, Real dy, Real dz,
        Real nx, Real ny, Real nz,
        Real minv, Real maxv,
        Real orientationSign ) const;

    void CalcGradient( const CellGradientView & view );

    // Derived temperature field and its cell gradient for the viscous and
    // heat-conduction slice.  Overrides the fail-closed FluxBackend defaults.
    bool SupportsCellTemperature() const override { return true; }
    void CalcCellTemperature( const CellTemperatureView & view ) override;

    // Viscous stress, energy and conduction flux for one stage.
    //
    // Mirrors the production host chain term for term:
    //   UNsVisFlux::CalcFaceVisFlux (CalcHeatFlux -> CalcStress ->
    //   CalcNsVisFlux), the VIS_STD face preparation in PrepareFaceValue
    //   (average gradient, CalcGradCoef, AverFaceValue, CorrectFaceGrad,
    //   CalcNormalGrad), NsVisFluxPolicy's stress/projection/energy/conductivity
    //   helpers, and UpdateFaceVisFlux's -A/Re scaling.
    //
    // VIS_STD (ivischeme == 0) is the only face scheme implemented: it is the
    // one every production case sets, and the other four differ only in how the
    // face normal gradient is built.  A caller cannot select another scheme
    // through ViscousFluxView, so the kernel cannot silently substitute one.
    //
    // The face tangents in the view are not consumed: GCom::CalcTangent computes
    // them for every face, but their only readers are VisGrad::ModifyFaceGrad
    // and the VIS_TEST/NEW1/NEW2 paths, which VIS_STD never enters.
    bool SupportsViscousFlux() const override { return true; }
    void CalcViscousFaceFlux( const ViscousFluxView & view ) override;

    // Wall heat flux needs the solid faces themselves.  ViscousFluxView carries
    // the output buffer and the count but deliberately no boundary
    // classification, because the accel layer does not know the solver's BC
    // ids.  The caller passes the face ids in the same order it would push the
    // host heat-flux vector, which is exactly the order
    // UNsVisFlux::SaveHeatFlux writes (increasing face id), so the result can be
    // compared element by element with the host field.
    //
    // The override above passes no ids, so it produces no wall heat flux; a view
    // that asks for one through that entry fails closed instead of writing
    // values for faces nobody classified.
    void CalcViscousFaceFlux(
        const ViscousFluxView & view, const int * wallFaceIds );

    // Single-value device mirror of SutherlandViscosityNondim
    // (codes/model/src/Sutherland.cpp), for the numeric contract test that
    // pins it to the host expression.  The nondimensional form is
    // sqrt( t * t * t ) * ( 1 + ts ) / ( t + ts ), not pow( t, 1.5 ).
    Real EvaluateDeviceSutherland(
        Real temperatureNondim, Real sutherlandTsNondim ) const;

    // Compute the local inviscid spectral radius and timestep from resident
    // primitive state. The resulting timestep can remain device-resident for
    // ScaleCurrentResidual.
    void CalcLocalInviscidTimeStep( const CellTimeStepView & view );

    // The time-step and state-update kernels report a non-physical or
    // non-finite result through a device flag that used to be copied back to
    // the host after every call, i.e. three blocking 4-byte round trips per
    // macro step on top of the launches. The default keeps the flag on the
    // device and checks once per macro step through CheckResidentStatus(); set
    // ONEFLOW_HIP_DEFER_STATUS=0 to check per call.
    void CheckResidentStatus();

    // Deferred variant: the macro step only *requests* the check.  The read is
    // performed at the next point that already drains the queue (the
    // host-staged interface pack), so a four-byte status read no longer costs
    // an extra full device drain per step.  Callers outside the interface path
    // fall back to CheckResidentStatus() through FlushRequestedResidentStatus().
    void RequestResidentStatusCheck();
    void FlushRequestedResidentStatus();

    // Update supported local inviscid boundary models on the resident state.
    void UpdateBoundaryState( const BoundaryStateView & view );

    // Reconstruct face primitive values into state-owned device buffers. The
    // host view supplies metadata and optional diagnostic destinations; the
    // numerical path consumes q/gradient/geometry already resident on device.
    // When hasSolidBoundary is set, bcQ must be one contiguous equation-major
    // buffer with stride nBoundaryFaces. Solid-override faces are a prefix,
    // so the kernel indexes that buffer by face.
    void ReconstructFaceValues(
        const CellFaceReconstructionView & view,
        bool downloadToHost = false );

    // Consume the device qf1/qf2 produced by ReconstructFaceValues and scatter
    // the resulting flux directly into the device residual. hostFlux is only
    // populated for explicit trace/diagnostic requests.  residualBase selects
    // whether the host residual field is uploaded as the accumulation base or
    // the device residual is zeroed instead.
    void CalcAndAddReconstructedFaceFlux(
        const FaceConnectivityView & connectivity,
        ResidualView & residual,
        int scheme,
        Real gamma,
        FaceFluxView * hostFlux = nullptr,
        const void * cacheKey = nullptr,
        bool downloadResidual = true,
        ResidualBase residualBase = ResidualBase::HostUpload,
        Real entropy1 = 0.0,
        Real entropy2 = 0.0,
        bool entropyFix = false );

    // Continue the reconstructed-flux path without a residual round trip.
    // ScaleCurrentResidual implements CALC_LHS. UpdatePrimitiveState implements
    // the five-equation primitive -> conserved -> primitive flow update and
    // mirrors the internal-cell primitive state to host when requested.
    void ScaleCurrentResidual( const CellStateUpdateView & view );
    // One entry for the two element-wise per-cell stages of the explicit
    // chain.  ScaleCurrentResidual and UpdatePrimitiveState each launched a
    // kernel, drained, and repeated the same binding/extent checks; the two
    // kernels have no host-visible state between them, so they can be issued
    // back to back from a single entry with a single check at the end.
    void ScaleAndUpdatePrimitiveState(
        const CellStateUpdateView & view,
        bool downloadToHost );

    void UpdatePrimitiveState(
        const CellStateUpdateView & view,
        bool downloadToHost = true );

    BoundaryBufferView PackBoundary(
        const void * cacheKey,
        const BoundaryFieldExchangePlan & plan,
        BoundaryTransferMode mode );
    BoundaryBufferView PackBoundary(
        const void * cacheKey,
        const BoundaryEndpointMap & endpoint,
        int nEquations,
        BoundaryTransferMode mode );
    BoundaryBufferView PackBoundary(
        const void * cacheKey,
        const BoundaryEndpointMap & endpoint,
        int nEquations,
        BoundaryTransferMode mode,
        const HipBoundaryEndpointKey & endpointKey );
    void UploadBoundaryState(
        const void * cacheKey,
        int nCells,
        int nGhostCells,
        int nEquations,
        const Real * const * components );
    BoundaryBufferView PrepareBoundaryReceive(
        const void * cacheKey,
        const BoundaryFieldExchangePlan & plan,
        BoundaryTransferMode mode );
    BoundaryBufferView PrepareBoundaryReceive(
        const void * cacheKey,
        const BoundaryEndpointMap & endpoint,
        int nEquations,
        BoundaryTransferMode mode );
    BoundaryBufferView PrepareBoundaryReceive(
        const void * cacheKey,
        const BoundaryEndpointMap & endpoint,
        int nEquations,
        BoundaryTransferMode mode,
        const HipBoundaryEndpointKey & endpointKey );
    void UnpackBoundary(
        const void * cacheKey,
        const BoundaryFieldExchangePlan & plan,
        BoundaryTransferMode mode );
    void UnpackBoundary(
        const void * cacheKey,
        const BoundaryEndpointMap & endpoint,
        int nEquations,
        BoundaryTransferMode mode );
    void UnpackBoundary(
        const void * cacheKey,
        const BoundaryEndpointMap & endpoint,
        int nEquations,
        BoundaryTransferMode mode,
        const HipBoundaryEndpointKey & endpointKey );

    void BindState( const void * cacheKey, Ns3DDeviceState & state );
    void UnbindState( const void * cacheKey ) noexcept;

    // --- Resident state ownership contract -------------------------------
    // The solver calls these when it writes the host mirror outside the
    // device path.  Owned-cell writes invalidate the whole resident field;
    // ghost writes only invalidate the ghost region.
    void NotifyHostFieldStateModified( const void * cacheKey );
    void NotifyHostGhostStateModified( const void * cacheKey );

    // Refresh the ghost region of the state-owned device q from the host
    // mirror when (and only when) the contract says it is stale.  This is the
    // host-staged halo binding: the exchanged host ghost remains the transfer
    // source of truth, while the device state owns the resident copy.
    bool SyncResidentGhostStateFromHost(
        const void * cacheKey,
        int nCells,
        int nGhostCells,
        int nEquations,
        const Real * const * components );

    // Mirror the resident gradient at the given cells into the solver-owned
    // host components.  A device reconstruction only needs the host mirror at
    // the boundary-adjacent owned cells (the interface send values and the
    // face-indexed boundary gradient), not the whole field: the interface
    // exchange, including its one-round lag, still runs through
    // Grad::SwapBcGrad on that mirror.
    void MirrorGradientCellsToHost(
        const void * cacheKey,
        int nEquations,
        const std::vector< int > & cellIds,
        Real * const * dqdx,
        Real * const * dqdy,
        Real * const * dqdz );

    // Pack the resident gradient at the interface send cells into a
    // [component][equation][face] host buffer.  The host then copies it into
    // the interface send storage in exchange order, so the endpoint cell
    // gather runs on the device instead of the host.  Returns false when this
    // zone has no resident gradient yet (first step, or a zone this rank has
    // not advanced); the caller then packs the host mirror, which is the
    // authoritative state at that point.
    bool PackGradientFacesToHost(
        const void * cacheKey,
        const int * cellIds,
        int nFaces,
        int nEquations,
        Real * hostPacked );

    // Read the interface receive storage into a packed host buffer and scatter
    // it directly into the resident ghost gradient.  This replaces the host
    // ghost scatter plus the separate ghost-gradient H2D sync.
    void UnpackGradientFacesFromHost(
        const void * cacheKey,
        const int * ghostCellIds,
        int nFaces,
        int nEquations,
        const Real * hostPacked );

    // Upload the interface receive payload for the device gradient into the
    // resident pack buffer, once per exchange round.  Every RK stage then
    // unpacks it on the device instead of repeating the same host->device copy.
    void UploadGradientFacesToDevice(
        const void * cacheKey,
        int nFaces,
        int nEquations,
        const Real * hostPacked );

    // Scatter the uploaded receive payload into the resident ghost gradient.
    // The payload buffer is the one UploadGradientFacesToDevice filled.
    void UnpackGradientFacesFromDevice(
        const void * cacheKey,
        const int * ghostCellIds,
        int nFaces,
        int nEquations );

    // Reduce the resident (already LHS-scaled) device residual into the
    // per-zone diagnostic scalars the solver's residual task reports.  Returns
    // false when the grid has no resident residual, so callers keep the host
    // loop as the fallback.
    bool ComputeResidualMetrics(
        const void * cacheKey,
        int nCells,
        int nEquations,
        Real gamma,
        HipResidualMetrics & metrics );

    // Diagnostic contract check: compare the complete resident primitive
    // state (owned + ghost cells) against the host mirror.  Returns false when
    // the grid has no bound resident state.
    bool VerifyResidentStateAgainstHost(
        const void * cacheKey,
        int nCells,
        int nGhostCells,
        int nEquations,
        const Real * const * components,
        Real tolerance,
        Real * maxAbsoluteDifference,
        long * mismatchedValues );

    // True when the resident device state is newer than the host mirror, i.e.
    // the device is authoritative for q and host readers must pull the values
    // back before using them.
    bool HostMirrorStale( const void * cacheKey ) const;

    // Pull the complete resident field (owned + ghost cells) back into the
    // host q mirror.  No-op when the mirror is already current.  Returns true
    // when values were transferred.
    bool RefreshHostPrimitiveState(
        const void * cacheKey,
        int nCells,
        int nGhostCells,
        int nEquations,
        Real * const * components );

    // Pull the resident field at the given owned cells back into the host
    // mirror.  Used by the per-step host readers that only touch
    // boundary-adjacent cells (aerodynamic force, pressure coefficient).
    void MirrorPrimitiveCellsToHost(
        const void * cacheKey,
        int nEquations,
        const int * cellIds,
        int cellCount,
        Real * const * components );

    // The device wrote cells the host mirror does not reflect (an interface
    // unpack whose paired host exchange was skipped).
    void MarkResidentHostMirrorStale( const void * cacheKey );

private:
    // Set by RequestResidentStatusCheck(); consumed by
    // FlushRequestedResidentStatus() at the next point that already drains.
    bool statusCheckPending_ = false;
    DeviceBuffer< Real > qLeft_;
    DeviceBuffer< Real > qRight_;
    DeviceBuffer< Real > xNormal_;
    DeviceBuffer< Real > yNormal_;
    DeviceBuffer< Real > zNormal_;
    DeviceBuffer< Real > meshVelocityNormal_;
    DeviceBuffer< Real > area_;
    DeviceBuffer< Real > deviceFlux_;
    DeviceBuffer< int > leftCell_;
    DeviceBuffer< int > rightCell_;
    DeviceBuffer< unsigned char > boundaryMask_;
    DeviceBuffer< Real > deviceResidual_;
    // Viscous-flux staging.  ViscousFluxView hands over host pointers (that is
    // the shape the backend-neutral contract defines, and the shape a CPU
    // implementation needs), so the HIP implementation uploads one stage's
    // inputs, runs one kernel and copies the face flux back.  These buffers are
    // grown and reused across stages; they are not resident solver state.
    DeviceBuffer< Real > viscousCellState_;    // nEquations * totalCells
    DeviceBuffer< Real > viscousCellGradient_; // nEquations*3 * totalCells
    DeviceBuffer< Real > viscousTemperature_;  // totalCells
    DeviceBuffer< Real > viscousTemperatureGradient_; // 3 * totalCells
    DeviceBuffer< Real > viscousLaminar_;      // totalCells
    DeviceBuffer< Real > viscousTurbulent_;    // totalCells
    DeviceBuffer< Real > viscousCellCenter_;   // 3 * totalCells
    DeviceBuffer< Real > viscousFaceNormals_;  // 3 * nFaces
    DeviceBuffer< Real > viscousFaceArea_;     // nFaces
    DeviceBuffer< int > viscousFaceCells_;     // 2 * nFaces
    // Equation-major face output with the caller's stride.  The highest index
    // written is ( nEquations - 1 ) * stride + ( nFaces - 1 ), so the buffer is
    // not stride * nFaces; see the allocation comment in the implementation.
    DeviceBuffer< Real > viscousFaceFlux_;
    // Wall heat flux, one value per classified solid face.  The view's buffer
    // is a host pointer, so the kernel writes this device buffer and the entry
    // copies it back; handing the host pointer to the kernel is a VMFault.
    DeviceBuffer< Real > viscousWallHeatFlux_; // nWallFaces
    HipGradientStorage fallbackGradientStorage_;
    std::unordered_map< const void *, Ns3DDeviceState * > boundStates_;
    int deviceFluxFaces_ = 0;
    int deviceFluxEquations_ = 0;
    const Real * cachedXNormal_ = nullptr;
    const Real * cachedYNormal_ = nullptr;
    const Real * cachedZNormal_ = nullptr;
    const Real * cachedMeshVelocityNormal_ = nullptr;
    const Real * cachedArea_ = nullptr;
    const int * cachedLeftCell_ = nullptr;
    const int * cachedRightCell_ = nullptr;
    const unsigned char * cachedBoundaryMask_ = nullptr;
    const void * cachedOwner_ = nullptr;
    int cachedFaces_ = 0;
    int cachedBoundaryFaces_ = -1;
    bool geometryCacheValid_ = false;
};

}
