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

#include "UNsInvFlux.h"
#include "UNsGrad.h"
#include "Zone.h"
#include "Atmosphere.h"
#include "UnsGrid.h"
#include "DataBase.h"
#include "FieldWrap.h"
#include "SolverDef.h"
#include "UCom.h"
#include "UNsCom.h"
#include "NsCom.h"
#include "NsIdx.h"
#include "HXMath.h"
#include "Fatal.h"
#include "Boundary.h"
#include "BcRecord.h"
#include "UNsLimiter.h"
#include "ULimiter.h"
#include "FieldManager.h"
#include "Iteration.h"
#include "TurbCom.h"
#include "UTurbCom.h"
#include "SolverDef.h"
#include "SolverState.h"
#include "Ctrl.h"
#include "ResidualBasePolicy.h"
#include "ZoneState.h"
#include "GridState.h"
#include "AccelRuntime.h"
#include "StageProfiler.h"
#include "CpuFluxBackend.h"
#include "EulerCpuAdapter.h"
#include "EulerInvFluxCapability.h"
#include "HipProductionPolicy.h"
#include "Parallel.h"
#ifdef ONEFLOW_ENABLE_HIP
#include "HipFluxBackend.h"
#endif
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <vector>


BeginNameSpace( ONEFLOW )

namespace
{

std::string RankAwareTracePath( const char * configured )
{
    std::string path( configured == nullptr ? "" : configured );
    const std::string token( "%r" );
    const char * rank = nullptr;
    const char * rankNames[] = {
        "OMPI_COMM_WORLD_RANK", "MV2_COMM_WORLD_RANK", "SLURM_PROCID" };
    for ( const char * name : rankNames )
    {
        const char * value = std::getenv( name );
        if ( value != nullptr && value[ 0 ] != '\0' )
        {
            rank = value;
            break;
        }
    }
    const std::string fallback = std::to_string( Parallel::GetPid() );
    if ( rank == nullptr ) rank = fallback.c_str();
    const std::string rankValue( rank );
    std::size_t tokenPos = path.find( token );
    while ( tokenPos != std::string::npos )
    {
        path.replace( tokenPos, token.size(), rankValue );
        tokenPos = path.find( token, tokenPos + rankValue.size() );
    }
    return path;
}

void AppendStageTraceRecord(
    const char * path,
    std::uint32_t kind,
    std::uint32_t sequence,
    std::uint64_t nItems,
    std::uint32_t nEquations,
    const std::vector< const MRField * > & fields )
{
    std::ofstream output( path, std::ios::binary | std::ios::app );
    if ( ! output )
    {
        throw std::runtime_error( "cannot open UNs stage trace file" );
    }

    const char magic[ 8 ] = "OFSTG01";
    const std::int32_t outerStep = Iteration::outerSteps;
    const std::int32_t gridLevel = GridState::gridLevel;
    const std::uint32_t nArrays =
        static_cast< std::uint32_t >( fields.size() );
    output.write( magic, sizeof( magic ) );
    output.write(
        reinterpret_cast< const char * >( & kind ), sizeof( kind ) );
    output.write(
        reinterpret_cast< const char * >( & sequence ), sizeof( sequence ) );
    output.write(
        reinterpret_cast< const char * >( & outerStep ), sizeof( outerStep ) );
    output.write(
        reinterpret_cast< const char * >( & gridLevel ), sizeof( gridLevel ) );
    output.write(
        reinterpret_cast< const char * >( & nEquations ),
        sizeof( nEquations ) );
    output.write(
        reinterpret_cast< const char * >( & nArrays ), sizeof( nArrays ) );
    output.write(
        reinterpret_cast< const char * >( & nItems ), sizeof( nItems ) );

    for ( const MRField * field : fields )
    {
        if ( field == nullptr || field->GetNEqu() < nEquations )
        {
            throw std::runtime_error(
                "UNs stage trace field has an invalid equation extent" );
        }
        for ( std::uint32_t equation = 0;
              equation < nEquations; ++ equation )
        {
            const auto & values = ( * field )[ equation ];
            if ( values.size() < nItems )
            {
                throw std::runtime_error(
                    "UNs stage trace field has an invalid item extent" );
            }
            output.write(
                reinterpret_cast< const char * >( values.data() ),
                static_cast< std::streamsize >(
                    nItems * sizeof( Real ) ) );
        }
    }
    if ( ! output )
    {
        throw std::runtime_error( "failed while writing UNs stage trace" );
    }
}

void AppendStageTraceValueRecord(
    const char * path,
    std::uint32_t kind,
    std::uint32_t sequence,
    std::uint64_t nItems,
    const std::vector< std::vector< Real > > & arrays )
{
    std::ofstream output( path, std::ios::binary | std::ios::app );
    if ( ! output )
    {
        throw std::runtime_error( "cannot open UNs stage trace file" );
    }

    const char magic[ 8 ] = "OFSTG01";
    const std::int32_t outerStep = Iteration::outerSteps;
    const std::int32_t gridLevel = GridState::gridLevel;
    const std::uint32_t nEquations = 1;
    const std::uint32_t nArrays =
        static_cast< std::uint32_t >( arrays.size() );
    output.write( magic, sizeof( magic ) );
    output.write(
        reinterpret_cast< const char * >( & kind ), sizeof( kind ) );
    output.write(
        reinterpret_cast< const char * >( & sequence ), sizeof( sequence ) );
    output.write(
        reinterpret_cast< const char * >( & outerStep ), sizeof( outerStep ) );
    output.write(
        reinterpret_cast< const char * >( & gridLevel ), sizeof( gridLevel ) );
    output.write(
        reinterpret_cast< const char * >( & nEquations ),
        sizeof( nEquations ) );
    output.write(
        reinterpret_cast< const char * >( & nArrays ), sizeof( nArrays ) );
    output.write(
        reinterpret_cast< const char * >( & nItems ), sizeof( nItems ) );

    for ( const auto & values : arrays )
    {
        if ( values.size() < nItems )
        {
            throw std::runtime_error(
                "UNs stage trace values have an invalid item extent" );
        }
        output.write(
            reinterpret_cast< const char * >( values.data() ),
            static_cast< std::streamsize >( nItems * sizeof( Real ) ) );
    }
    if ( ! output )
    {
        throw std::runtime_error( "failed while writing UNs stage trace" );
    }
}

}

UNsInvFlux::UNsInvFlux()
{
    limiter = std::make_unique<NsLimiter>();
    limf = limiter->limf.get();
}

UNsInvFlux::~UNsInvFlux()
{
}

void UNsInvFlux::CalcLimiter()
{
    limiter->CalcLimiter();
}

namespace {

// Diagnostic for the limiter work (B2): write the per-cell factor and the
// reconstructed face values once per run so a CPU run and a HIP run of the same
// case can be compared face by face.  Never enabled on a validated path.
void DumpLimiterState( const char * side, LimField * limf, UnsGrid * grid,
                       int nEquations, bool deviceResident )
{
    const char * dir = std::getenv( "ONEFLOW_UNS_LIMITER_DUMP" );
    if ( dir == nullptr || dir[ 0 ] == '\0' ) return;
    static bool alreadyDone = false;
    if ( alreadyDone ) return;
    alreadyDone = true;

    const int nCells = grid->nCells;
    const int nFaces = grid->nFaces;
    const std::size_t phiCount =
        static_cast< std::size_t >( nEquations ) * nCells;
    const std::size_t faceCount =
        static_cast< std::size_t >( nEquations ) * nFaces;
    std::vector< Real > phi( phiCount, 0.0 );
    std::vector< Real > qf1( faceCount, 0.0 );
    std::vector< Real > qf2( faceCount, 0.0 );

    if ( deviceResident )
    {
#ifdef ONEFLOW_ENABLE_HIP
        HipFluxBackend::Shared().CopyLimiterStateToHost(
            grid, nCells, nFaces, nEquations,
            phi.data(), qf1.data(), qf2.data() );
#endif
    }
    else
    {
        for ( int equation = 0; equation < nEquations; ++ equation )
        {
            for ( int cell = 0; cell < nCells; ++ cell )
                phi[ static_cast< std::size_t >( equation ) * nCells + cell ] =
                    ( * limf->limiter )[ equation ][ cell ];
            for ( int face = 0; face < nFaces; ++ face )
            {
                qf1[ static_cast< std::size_t >( equation ) * nFaces + face ] =
                    ( * limf->qf1 )[ equation ][ face ];
                qf2[ static_cast< std::size_t >( equation ) * nFaces + face ] =
                    ( * limf->qf2 )[ equation ][ face ];
            }
        }
    }

    auto write = [ dir, side ]( const char * suffix,
                                const std::vector< Real > & values,
                                int nEquations, int count )
    {
        std::string path = std::string( dir ) + "/" + side + "-" + suffix;
        std::ofstream out( path, std::ios::out | std::ios::binary );
        const int header[ 2 ] = { nEquations, count };
        out.write( reinterpret_cast< const char * >( header ), sizeof( header ) );
        out.write( reinterpret_cast< const char * >( values.data() ),
                   static_cast< std::streamsize >(
                       values.size() * sizeof( Real ) ) );
    };
    write( "phi.bin", phi, nEquations, nCells );
    write( "qf1.bin", qf1, nEquations, nFaces );
    write( "qf2.bin", qf2, nEquations, nFaces );
}

}

void UNsInvFlux::CalcInvFace()
{
    // Optional CPU-vs-HIP numerical trace at the RHS entry (off by default).
    if ( std::getenv( "ONEFLOW_TRACE_POINTS" ) != nullptr
         && Iteration::outerSteps <= 2 )
    {
        static int invPrintCount = 0;
        if ( invPrintCount < 8 )
        {
            ++invPrintCount;
            const int probeCell = 254656;
            MRField * qF = GetFieldPointer< MRField >( Zone::GetUnsGrid(), "q" );
            MRField * dqF = GetFieldPointer< MRField >( Zone::GetUnsGrid(), "dq" );
            std::cout << std::setprecision( 17 );
            std::cout << "TRACE_INVF step=" << Iteration::outerSteps
                      << " inner=" << Iteration::innerSteps
                      << " cfl=" << Iteration::cfl
                      << " lhscoef=" << ctrl.lhscoef
                      << " rkStage=" << ctrl.rkStage
                      << " hipGrad=" << ( UseHipGradient() ? 1 : 0 )
                      << " hipRecon=" << ( UseHipDeviceReconstruction() ? 1 : 0 )
                      << " hipBatch=" << ( UseHipBatchAdapter() ? 1 : 0 );
            if ( qF != nullptr )
            {
                for ( int eq = 0; eq < 5 && eq < qF->GetNEqu(); ++eq )
                    std::cout << " q" << eq << "=" << (*qF)[eq][probeCell];
            }
            if ( dqF != nullptr )
            {
                for ( int eq = 0; eq < 5 && eq < dqF->GetNEqu(); ++eq )
                    std::cout << " dq" << eq << "=" << (*dqF)[eq][probeCell];
            }
            std::cout << std::endl;
        }
    }
    {
        ScopedStageTimer timer( "gradient" );
        ScopedHostTimer hostGradientTimer( "s_gradient" );
        uns_grad.Init();
        const bool hipGradient = this->UseHipGradient();
        const bool hipReconstruction = this->UseHipDeviceReconstruction();
        // A cross-zone device reconstruction consumes ghost gradients that the
        // interface exchange owns.  The exchange is host-staged, so the host
        // mirror stays the authoritative source for that ghost slice: download
        // the gradient, let Grad::SwapBcGrad apply the interface exchange, and
        // upload only the ghost slice back (see
        // CellFaceReconstructionView::syncGradientGhostsFromHost).
        const bool crossZoneGhostGradient = hipReconstruction
            && ( ZoneState::nZones > 1 || ZoneState::nLocal > 1 );
        // packedExchange is the contract that fills the resident ghost
        // gradient every stage; mirrorToHost only feeds host readers.
        const bool packedExchange =
            hipGradient && hipReconstruction && crossZoneGhostGradient;
        const char * forceDownload =
            std::getenv( "ONEFLOW_GRAD_FORCE_DOWNLOAD" );
        const bool downloadGradient =
            ( forceDownload != nullptr && forceDownload[ 0 ] == '1' )
            || ( hipGradient && ( ! hipReconstruction || crossZoneGhostGradient )
            && ! ( packedExchange && HipPackedGradientExchangeEnabled() ) )
            // A viscous case reads the host cell gradients when assembling the
            // viscous flux (CalcVisFluxHip view.dqdx = unsf.dqdx, and the host
            // CalcVisFlux path).  With device reconstruction active on a single
            // zone the first clause above is false, leaving those host gradients
            // stale and the wall-layer residual ~12% off.  Download the resident
            // gradient for viscous configs so the reader sees the current value.
            || ( hipGradient && vis_model.vismodel != 0 );
        // With a cross-zone device reconstruction the host mirror is only used
        // by the host interface exchange, so mirror just the boundary-adjacent
        // owned gradient instead of downloading the whole field.  On the
        // device-packed path the interface faces travel through the packed
        // buffer instead, which leaves the host consumers (aero forces,
        // boundary-gradient store, field dumps) as the only readers.  Those
        // readers run after the RK loop, so only the last stage's mirror is
        // visible to them.  InnerOk covers the pressure-coefficient task,
        // which is independent of nForceSave.
        // Only two host readers consume the boundary-adjacent gradient: the
        // aerodynamic force (nForceSave cadence), the pressure-coefficient /
        // visualization dump (nVisualSave cadence) and the final step.  The
        // steady inner loop itself does not, so InnerOk is not a reader.
        const bool visualDue =
            Iteration::nVisualSave > 0
            && Iteration::outerSteps % Iteration::nVisualSave == 0;
        const bool postStepGradientReader =
            Iteration::ForceOk()
            || visualDue
            || Iteration::outerSteps + 1 >= Iteration::maxSteps;
        const bool hostMirrorReaderThisStep = HipGradientMirrorOnThisStage(
            HipGradientHostMirrorOnDemand(),
            ctrl.debugFieldDump != 0,
            postStepGradientReader,
            ctrl.rkStage,
            static_cast< int >( ctrl.rk_coef.size() ) );
        const bool mirrorBoundaryGradient =
            packedExchange
            && ( ! HipPackedGradientExchangeEnabled() || hostMirrorReaderThisStep );
        if ( std::getenv( "ONEFLOW_UNS_HIP_STAGE_DIAGNOSTIC" ) != nullptr )
        {
            static int printed = 0;
            if ( printed < 4 )
            {
                std::cout << "HIP_STAGE_DIAGNOSTIC hip_gradient=" << hipGradient
                          << " hip_reconstruction=" << hipReconstruction
                          << " download_gradient=" << downloadGradient
                          << " mirror_boundary_gradient="
                          << mirrorBoundaryGradient
                          << " zones=" << ZoneState::nZones
                          << " local_zones=" << ZoneState::nLocal
                          << std::endl;
                ++printed;
            }
        }
        if ( hipGradient )
            uns_grad.CalcGradHip(
                downloadGradient && ! packedExchange,
                packedExchange,
                mirrorBoundaryGradient );
        else
            uns_grad.CalcGrad();
        // Optional full-gradient dump (off by default) to compare the device
        // Green-Gauss result against the host oracle cell by cell.
        {
            const char * gradDump =
                std::getenv( "ONEFLOW_GRAD_DUMP_FILE" );
            if ( gradDump != nullptr && gradDump[ 0 ] != '\0' )
            {
                // Overwrite on every gradient pass so the dump reflects the
                // same last-stage moment as the qf trace (both used to pin
                // the reconstruction divergence).
                std::ofstream gout( gradDump, std::ios::binary | std::ios::trunc );
                const int nEqu = uns_grad.nEqu;
                const int nCells = ug.nCells;
                gout.write( reinterpret_cast< const char * >( &nEqu ), sizeof(int) );
                gout.write( reinterpret_cast< const char * >( &nCells ), sizeof(int) );
                for ( int eq = 0; eq < nEqu; ++ eq )
                {
                    for ( int c = 0; c < nCells; ++ c )
                    {
                        const double v[3] = { (*uns_grad.dqdx)[eq][c],
                                              (*uns_grad.dqdy)[eq][c],
                                              (*uns_grad.dqdz)[eq][c] };
                        gout.write( reinterpret_cast< const char * >( v ), sizeof(v) );
                    }
                }
            }
        }
    }
    // The device reconstruction computes the Barth factor from the resident
    // gradient. Running the host limiter as well repeated the whole face loop on
    // the host every stage for nothing (measured: 270 ms per call, which turned
    // an 88x device speedup into 2.9x), so skip it exactly when the device owns
    // the factor.
    const bool deviceLimiterOwnsFactor =
        ctrl.ilim == ILMT_BARTH && this->UseHipDeviceReconstruction();
    if ( deviceLimiterOwnsFactor )
    {
        // Host readers (restart/visualization dumps, the interface storage) still
        // expect an initialised field; one fill is enough, no face loop.
        limiter->FillNeutralLimiterField();
    }
    else
    {
        ScopedStageTimer timer( "limiter" );
        this->CalcLimiter();
    }
    if ( this->UseHipDeviceReconstruction() )
    {
#ifdef ONEFLOW_ENABLE_HIP
        ScopedStageTimer timer( "reconstruction" );
        UnsGrid * grid = Zone::GetUnsGrid();
        const int nFaces = ug.nFaces;
        const int nEquations = limf->nEqu;
        const bool topologyChanged =
            hipBoundaryGrid != grid
            || hipBoundaryRecord != ug.bcRecord
            || hipBoundaryFaces != nFaces
            || hipBoundaryBoundaryFaces != ug.nBFaces
            || hipBoundaryEquations != nEquations
            || hipBoundaryOperation.size()
                   != static_cast< std::size_t >( nFaces );
        if ( topologyChanged )
        {
            hipBoundaryOperation.assign(
                nFaces, ReconstructionBoundaryOperation::Preserve );
            hipBoundaryGhostPolicy.assign(
                nFaces, BoundaryGhostUpdatePolicy::ExternalHalo );
            hipBoundaryHasSolid = false;
            for ( int face = 0; face < ug.nBFaces; ++ face )
            {
                const int bcType = ug.bcRecord->bcType[ face ];
                if ( bcType == BC::INTERFACE || bcType == BC::PERIODIC )
                {
                    hipBoundaryOperation[ face ] =
                        ReconstructionBoundaryOperation::Preserve;
                }
                else if ( bcType == BC::SOLID_SURFACE )
                {
                    hipBoundaryOperation[ face ] =
                        ReconstructionBoundaryOperation::SolidOverride;
                    hipBoundaryHasSolid = true;
                }
                else
                {
                    hipBoundaryOperation[ face ] =
                        ReconstructionBoundaryOperation::Average;
                }
                if ( bcType != BC::INTERFACE && bcType != BC::PERIODIC )
                {
                    hipBoundaryGhostPolicy[ face ] =
                        BoundaryGhostUpdatePolicy::LocalBoundaryModel;
                }
            }
            hipBoundaryGrid = grid;
            hipBoundaryRecord = ug.bcRecord;
            hipBoundaryFaces = nFaces;
            hipBoundaryBoundaryFaces = ug.nBFaces;
            hipBoundaryEquations = nEquations;
        }
        if ( hipBoundaryHasSolid )
        {
            if ( unsf.bc_q == nullptr )
            {
                throw std::runtime_error(
                    "HIP reconstruction solid boundary is missing bc_q" );
            }
            // Solid-override boundary faces are a prefix of the face array, so
            // the wall state only needs nEquations*nBFaces values.  Packing the
            // full nFaces stride moved ~13x more data than the kernel reads.
            const int bcQFaces = static_cast< int >( ug.nBFaces );
            hipBoundaryQ.resize(
                static_cast< std::size_t >( nEquations ) * bcQFaces );
            for ( int equation = 0; equation < nEquations; ++ equation )
            {
                for ( int face = 0; face < bcQFaces; ++ face )
                {
                    hipBoundaryQ[ equation * bcQFaces + face ] =
                        ( * unsf.bc_q )[ equation ][ face ];
                }
            }
        }
        CellFaceReconstructionView view;
        view.nCells = ug.nCells;
        view.nGhostCells = ug.nTCell - ug.nCells;
        view.nFaces = nFaces;
        view.nBoundaryFaces = ug.nBFaces;
        view.nEquations = nEquations;
        view.xFace = ug.xfc->data();
        view.yFace = ug.yfc->data();
        view.zFace = ug.zfc->data();
        view.xCell = ug.xcc->data();
        view.yCell = ug.ycc->data();
        view.zCell = ug.zcc->data();
        view.leftCell = ug.lcf->data();
        view.rightCell = ug.rcf->data();
        view.boundaryGhostPolicy = hipBoundaryGhostPolicy.data();
        view.boundaryOperation = hipBoundaryOperation.data();
        // Barth/Venkat factors come from Limiter::CalcLimiter, which the host
        // runs just above for this stage; the device reconstruction consumes
        // that field instead of recomputing it.
        const bool deviceLimiterActive = ctrl.ilim == ILMT_BARTH;
        view.limiterMode = deviceLimiterActive
            ? ReconstructionLimiterMode::DeviceCell
            : ReconstructionLimiterMode::Disabled;
        view.physicality = ReconstructionPhysicalityPolicy::PositiveDensityPressure;
        view.densityComponent = 0;
        view.pressureComponent = nEquations == 5 ? 4 : 2;
        view.hasSolidBoundary = hipBoundaryHasSolid;
        view.syncGradientGhostsFromHost =
            ( ZoneState::nZones > 1 || ZoneState::nLocal > 1 )
            && ! HipPackedGradientExchangeEnabled();
        view.ownerToken = grid;
        view.cacheKey = grid;
        view.topologyGeneration = 1;
        view.fieldGeneration = 1;
        for ( int equation = 0; equation < nEquations; ++ equation )
        {
            view.q[ equation ] = ( * limf->q )[ equation ].data();
            view.dqdx[ equation ] = ( * limf->dqdx )[ equation ].data();
            view.dqdy[ equation ] = ( * limf->dqdy )[ equation ].data();
            view.dqdz[ equation ] = ( * limf->dqdz )[ equation ].data();
            view.qLeft[ equation ] = ( * limf->qf1 )[ equation ].data();
            view.qRight[ equation ] = ( * limf->qf2 )[ equation ].data();
            if ( hipBoundaryHasSolid )
                view.bcQ[ equation ] = hipBoundaryQ.data()
                    + equation * static_cast< int >( ug.nBFaces );
        }
        const bool needFaceTrace =
            ( std::getenv( "ONEFLOW_UNS_TRACE_FILE" ) != nullptr
              && std::getenv( "ONEFLOW_UNS_TRACE_FILE" )[ 0 ] != '\0' )
            || ( std::getenv( "ONEFLOW_UNS_STAGE_TRACE_FILE" ) != nullptr
                 && std::getenv( "ONEFLOW_UNS_STAGE_TRACE_FILE" )[ 0 ] != '\0' );
        HipFluxBackend::Shared().ReconstructFaceValues(
            view, needFaceTrace );
        DumpLimiterState( "dev", limf, grid, nEquations, true );
#else
        throw std::runtime_error(
            "HIP device reconstruction was requested without HIP support" );
#endif
        return;
    }
    { ScopedStageTimer timer( "reconstruction" );
        ScopedHostTimer hostReconTimer( "s_reconstruction" );
        this->GetQlQrField();
        this->ReconstructFaceValueField();
    }
    // Optional host face-reconstruction trace at face 452808 (off by default).
    // Prints qf1/qf2 and the owner/neighbour q + gradients that feed the
    // reconstruction, to test whether identical inputs still make different qf.
    if ( std::getenv( "ONEFLOW_TRACE_POINTS" ) != nullptr
         && Iteration::outerSteps <= 2 && limf != nullptr )
    {
        static int faceTraceCount = 0;
        if ( faceTraceCount < 12 )
        {
            ++faceTraceCount;
            const int face = 452808;
            const int nFaces = ug.nFaces;
            if ( face < nFaces && ug.lcf != nullptr && ug.rcf != nullptr )
            {
                const int lc = ( *ug.lcf )[ face ];
                const int rc = ( *ug.rcf )[ face ];
                std::cout << std::setprecision( 17 );
                std::cout << "TRACE_FACE step=" << Iteration::outerSteps
                          << " face=" << face << " lc=" << lc << " rc=" << rc;
                for ( int eq = 0; eq < 5 && eq < limf->nEqu; ++eq )
                {
                    std::cout << " qf1_" << eq << "=" << ( *limf->qf1 )[ eq ][ face ]
                              << " qf2_" << eq << "=" << ( *limf->qf2 )[ eq ][ face ];
                    std::cout << " ql_" << eq << "=" << ( *limf->q )[ eq ][ lc ]
                              << " qr_" << eq << "=" << ( *limf->q )[ eq ][ rc ];
                    std::cout << " gxl_" << eq << "=" << ( *limf->dqdx )[ eq ][ lc ]
                              << " gyl_" << eq << "=" << ( *limf->dqdy )[ eq ][ lc ]
                              << " gzl_" << eq << "=" << ( *limf->dqdz )[ eq ][ lc ];
                    if ( eq == 0 )
                        std::cout << " phi_l=" << ( *limf->limiter )[ eq ][ lc ]
                                  << " phi_r=" << ( *limf->limiter )[ eq ][ rc ];
                }
                std::cout << std::endl;
            }
        }
    }
    if ( limf != nullptr )
    {
        DumpLimiterState( "host", limf, Zone::GetUnsGrid(),
                          limf->nEqu, false );
    }
    { ScopedStageTimer timer( "boundary_reconstruction" );
        this->BoundaryQlQrFixField();
    }
}

void UNsInvFlux::GetQlQrField()
{
    this->limiter->GetQlQr();
}

void UNsInvFlux::ReconstructFaceValueField()
{
    this->limiter->CalcFaceValue();
}

void UNsInvFlux::BoundaryQlQrFixField()
{
    this->limiter->BcQlQrFix();
}

void UNsInvFlux::CalcFlux()
{
    if ( nscom.icmpInv == 0 ) return;
    // Host accounting for the RK-stage interior.  The parent
    // macro_update_residuals scope (NsResidentMacroStep) covers 5.3 ms/step on
    // the production configuration but contained no ScopedHostTimer at all, so
    // its breakdown was unreadable: the A4 budget needs to know whether that
    // time sits in host-side loops, in blocking copies, or in queue drains.
    // These scopes are children of macro_update_residuals and add no
    // synchronization of their own (ScopedHostTimer only reads steady_clock).
    inv.Init();
    ug.Init();
    unsf.Init();

    {
        ScopedHostTimer hostTimer( "s_init" );
        this->Alloc();
        this->SetPointer( nscom.ischeme );
    }

    this->CalcInvFace();
    if ( this->UseHipBatchAdapter() )
    {
        ScopedStageTimer fluxTimer( "flux_assembly" );
        ScopedHostTimer hostFluxTimer( "s_flux_assembly" );
        this->CalcAndAddInvFluxHipBatch();
        this->DumpInvFluxTrace();
    }
    else
    {
        ScopedHostTimer hostFluxTimer( "s_flux_host" );
        this->CalcInvFlux();
        this->DumpInvFluxTrace();
        this->AddInvFlux();
    }
    // Optional per-term residual trace after inviscid accumulation (off by default).
    if ( std::getenv( "ONEFLOW_TRACE_POINTS" ) != nullptr
         && Iteration::outerSteps <= 2 )
    {
        static int invResCount = 0;
        if ( invResCount < 30 )
        {
            ++invResCount;
            MRField * resF = FieldHome::GetResidualField( Zone::GetUnsGrid(), NS_SOLVER );
            if ( resF != nullptr )
            {
                const int probeCell = 254656;
                std::cout << std::setprecision( 17 );
                std::cout << "TRACE_RES_INV step=" << Iteration::outerSteps
                          << " cId=" << probeCell;
                for ( int eq = 0; eq < 5 && eq < resF->GetNEqu(); ++eq )
                    std::cout << " r" << eq << "=" << (*resF)[eq][probeCell];
                std::cout << std::endl;
            }
        }
    }
    this->DumpInvFluxStageTrace();

}

void UNsInvFlux::CalcInvFlux()
{
    if ( this->UseHipBatchAdapter() )
    {
        this->CalcInvFluxHipBatch();
        return;
    }
    if ( this->UseCpuBatchAdapter() )
    {
        this->CalcInvFluxCpuBatch();
        return;
    }

    for ( int fId = 0; fId < ug.nFaces; ++ fId )
    {
        ug.fId = fId;

        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        this->PrepareFaceValue();

        ( this->*invFluxPointer )();

        this->UpdateFaceInvFlux();
    }
}

bool UNsInvFlux::UseCpuBatchAdapter() const
{
    const char * enabled = std::getenv( "ONEFLOW_ENABLE_UNS_CPU_BATCH" );
    if ( enabled == nullptr || enabled[ 0 ] != '1' ) return false;
    if ( AccelRuntime::Instance().IsAccelerator() ) return false;
    return nscom.ischeme == ISCHEME_LAX_FRIEDRICHS
        && nscom.nEqu == 5 && !limiter->limfIsNullPtr() && limiter->GetNEquations() == 5;
}

bool UNsInvFlux::UseHipGradient() const
{
    const char * feature = "ONEFLOW_ENABLE_UNS_HIP_GRADIENT";
    if ( ! HipProductionFeatureEnabled( feature ) ) return false;
    if ( ! this->UseHipBatchAdapter() )
    {
        if ( ! HipProductionFeatureExplicitlyEnabled( feature ) )
            return false;
        throw std::runtime_error(
            "ONEFLOW_ENABLE_UNS_HIP_GRADIENT=1 requires the HIP batch path" );
    }
    // Barth is computed on the device from the resident gradient.  It must not
    // be computed on the host: with a device-authoritative q the host gradient
    // mirror is deliberately stale, and a host-side factor biased the solution
    // by 8.5e-5 over 200 steps.  Venkat is not implemented on the device yet.
    if ( ctrl.ilim != ILMT_NO && ctrl.ilim != ILMT_BARTH )
    {
        if ( ! HipProductionFeatureExplicitlyEnabled( feature ) )
            return false;
        throw std::runtime_error(
            "ONEFLOW_ENABLE_UNS_HIP_GRADIENT=1 supports limiter off or Barth "
            "on the device" );
    }
    if ( vis_model.vismodel != 0 )
    {
        // The NS gradient kernel is the same five-equation Green-Gauss form
        // for inviscid and viscous cases; SA working-variable gradients are
        // separate fields and are not consumed here.  Allow viscous configs
        // (e.g. the Roe/SA flagship) to opt in explicitly.
        const char * allowViscous =
            std::getenv( "ONEFLOW_ALLOW_VISCOUS_HIP_GRADIENT" );
        if ( allowViscous != nullptr && allowViscous[ 0 ] == '1' )
            return true;
        if ( ! HipProductionFeatureExplicitlyEnabled( feature ) )
            return false;
        throw std::runtime_error(
            "ONEFLOW_ENABLE_UNS_HIP_GRADIENT=1 requires inviscid mode "
            "(set ONEFLOW_ALLOW_VISCOUS_HIP_GRADIENT=1 for Roe/SA)" );
    }
    // Cross-zone gradient needs the resident halo contract: the interface
    // exchange must publish neighbor owned-cell values into the state-owned
    // device q ghost region, and the resident state generation contract must
    // keep host/device ghosts identical.  Without that interface path the
    // device gradient would read stale halo values.
    if ( ( ZoneState::nZones > 1 || ZoneState::nLocal > 1 )
         && ! HipProductionFeatureEnabled(
                "ONEFLOW_ENABLE_UNS_HIP_INTERFACE" ) )
    {
        if ( ! HipProductionFeatureExplicitlyEnabled( feature ) )
            return false;
        throw std::runtime_error(
            "ONEFLOW_ENABLE_UNS_HIP_GRADIENT=1 in a multi-zone domain "
            "requires ONEFLOW_ENABLE_UNS_HIP_INTERFACE=1" );
    }
    return true;
}

bool UNsInvFlux::UseHipDeviceReconstruction() const
{
    if ( ! HipProductionFeatureEnabled(
             "ONEFLOW_ENABLE_UNS_HIP_RECONSTRUCTION" ) )
        return false;
    if ( ! this->UseHipGradient() ) return false;
    // A partitioned multi-zone grid also has cross-rank interface ghosts. The
    // resident q ghost is supplied by the interface exchange contract, but the
    // ghost *gradient* at Preserve faces comes from the neighbor zone and is
    // exchanged by Grad::SwapBcGrad on the host.  The reconstruction view
    // therefore refreshes the device ghost gradient from the host mirror in
    // that case, and the interface contract must be active.
    if ( ( ZoneState::nZones > 1 || ZoneState::nLocal > 1 )
         && ! HipProductionFeatureEnabled(
                "ONEFLOW_ENABLE_UNS_HIP_INTERFACE" ) )
    {
        return false;
    }
    return true;
}

bool UNsInvFlux::UseHipDeviceStateUpdate() const
{
    const char * feature = "ONEFLOW_ENABLE_UNS_HIP_STATE_UPDATE";
    const bool explicitRequest =
        HipProductionFeatureExplicitlyEnabled( feature );
    if ( ! HipProductionFeatureEnabled( feature ) ) return false;
    if ( ! this->UseHipDeviceReconstruction() )
    {
        if ( ! explicitRequest ) return false;
        throw std::runtime_error(
            "ONEFLOW_ENABLE_UNS_HIP_STATE_UPDATE=1 requires HIP "
            "device reconstruction" );
    }
    // Keep this predicate in sync with UNsUpdate::UseHipDeviceStateUpdate.
    // The flagship (vismodel!=0) cannot run the device state update, so the
    // residual must be downloaded to the host every step: otherwise the host
    // res field never receives the device-computed inviscid flux (only the
    // viscous AddF2CField lands there), and LU-SGS / residual output / aero
    // read a partial host residual whose continuity row stays zero.
    if ( vis_model.vismodel != 0 )
    {
        if ( ! explicitRequest ) return false;
        throw std::runtime_error(
            "HIP state update does not support a viscous or turbulence model" );
    }
    return true;
}

bool UNsInvFlux::UseHipBatchAdapter() const
{
    const char * feature = "ONEFLOW_ENABLE_UNS_HIP_BATCH";
    if ( ! HipProductionFeatureEnabled( feature ) ) return false;

    AccelRuntime & runtime = AccelRuntime::Instance();
    EulerInvFluxCapabilityRequest request;
#ifdef ONEFLOW_ENABLE_HIP
    request.hipBuildEnabled = true;
#endif
    request.runtimeInitialized = runtime.IsInitialized();
    request.hipBackendSelected = request.runtimeInitialized
        && runtime.Backend().Kind() == AccelBackendKind::HIP;
    request.supportedSolver = SolverState::solverType == NS_SOLVER;
    request.localZoneCount = ZoneState::nLocal;
    request.gridLevel = GridState::gridLevel;
    request.gridCount = GridState::nGrids;
    request.nEquations = nscom.nEqu;
    request.limiterEquations = limf == nullptr ? 0 : limf->nEqu;
    request.laxFriedrichsScheme =
        nscom.ischeme == ISCHEME_LAX_FRIEDRICHS;
    request.roeScheme = nscom.ischeme == ISCHEME_ROE;
    // The fused flux kernel computes only the inviscid Riemann flux; the
    // viscous flux is assembled separately by UNsVisFlux.  On a viscous case
    // the no-flux device path is disabled by default.  Explicit opt-in lets
    // the flagship (Roe/SA) run its inviscid flux on the device while the
    // viscous term keeps its own (host or HIP-hook) path.
    const char * allowViscous =
        std::getenv( "ONEFLOW_ALLOW_VISCOUS_HIP_BATCH" );
    request.viscousModelActive =
        vis_model.vismodel != 0
        && ! ( allowViscous != nullptr && allowViscous[ 0 ] == '1' );
    request.hasInterfaceExchange =
        ZoneState::nZones > 1 || ZoneState::nLocal > 1;
    request.backendSupportsInterfaceExchange =
        HipProductionFeatureEnabled( "ONEFLOW_ENABLE_UNS_HIP_INTERFACE" );

    const EulerInvFluxCapabilityDecision decision =
        EvaluateEulerInvFluxCapability( request );
    if ( ! decision.enabled )
    {
        if ( ! HipProductionFeatureExplicitlyEnabled( feature ) )
            return false;
        throw std::runtime_error(
            std::string( "ONEFLOW_ENABLE_UNS_HIP_BATCH=1 rejected: " )
            + EulerInvFluxCapabilityReasonName( decision.reason ) );
    }
    return true;
}

void UNsInvFlux::CalcInvFluxCpuBatch()
{
    CpuFluxBackend backend;
    this->CalcInvFluxBatch( backend );
}

void UNsInvFlux::CalcInvFluxHipBatch()
{
#ifdef ONEFLOW_ENABLE_HIP
    HipFluxBackend & backend = HipFluxBackend::Shared();
    this->CalcInvFluxBatch( backend );
#else
    throw std::runtime_error(
        "UNs HIP batch path was called without HIP support." );
#endif
}

void UNsInvFlux::CalcAndAddInvFluxHipBatch()
{
#ifdef ONEFLOW_ENABLE_HIP
    UnsGrid * grid = Zone::GetUnsGrid();
    MRField * res = GetFieldPointer< MRField >( grid, "res" );
    if ( res == nullptr
         || res->GetNEqu() < static_cast< HXSize_t >( nscom.nEqu ) )
    {
        throw std::runtime_error(
            "UNs HIP batch residual field is unavailable." );
    }

    const int nFaces = ug.nFaces;
    const int nCells = ug.nCells;
    const int nEquations = nscom.nEqu;
    for ( int equation = 0; equation < nEquations; ++ equation )
    {
        if ( ( * res )[ equation ].size()
             < static_cast< std::size_t >( nCells ) )
        {
            throw std::runtime_error(
                "UNs HIP batch residual extent is invalid." );
        }
    }

    const char * traceFile = std::getenv( "ONEFLOW_UNS_TRACE_FILE" );
    const char * stageTraceFile =
        std::getenv( "ONEFLOW_UNS_STAGE_TRACE_FILE" );
    const bool needFlux =
        ( traceFile != nullptr && traceFile[ 0 ] != '\0' )
        || ( stageTraceFile != nullptr && stageTraceFile[ 0 ] != '\0' );
    if ( this->UseHipDeviceReconstruction() )
    {
        FaceConnectivityView connectivity;
        connectivity.nFaces = nFaces;
        connectivity.nBoundaryFaces = ug.nBFaces;
        connectivity.leftCell = ug.lcf->data();
        connectivity.rightCell = ug.rcf->data();

        ResidualView residual;
        residual.nCells = nCells;
        residual.nEquations = nEquations;
        residual.values = nullptr;
        for ( int equation = 0; equation < nEquations; ++ equation )
            residual.components[ equation ] = ( * res )[ equation ].data();

        std::vector< Real > faceFlux;
        FaceFluxView hostFlux;
        FaceFluxView * hostFluxPointer = nullptr;
        if ( needFlux )
        {
            faceFlux.resize( nEquations * nFaces );
            hostFlux.nFaces = nFaces;
            hostFlux.nEquations = nEquations;
            hostFlux.values = faceFlux.data();
            hostFluxPointer = & hostFlux;
        }
        HipFluxBackend & backend = HipFluxBackend::Shared();
        // The residual field itself only has to reach the host when a
        // diagnostic trace, a field/visualization dump or the final restart
        // dump needs it. Everything else consumes the device reduction.
        const bool dumpDue =
            ( Iteration::nFieldSave > 0
              && Iteration::outerSteps % Iteration::nFieldSave == 0 )
            || ( Iteration::nVisualSave > 0
                 && Iteration::outerSteps % Iteration::nVisualSave == 0 )
            || Iteration::outerSteps >= Iteration::maxSteps;
        // Optional host-residual trace around the device flux accumulation:
        // prints cell 254656's residual before and after so an equation that
        // becomes exactly zero can be attributed to the base upload vs the
        // scatter/accumulate kernel (off by default).
        const bool resTrace = std::getenv( "ONEFLOW_RES_DUMP" ) != nullptr;
        auto dumpRes = [&]( const char * tag )
        {
            if ( ! resTrace ) return;
            static int resTraceCount = 0;
            if ( resTraceCount >= 12 ) return;
            ++resTraceCount;
            const int probe = 254656;
            std::cout << std::setprecision( 17 );
            std::cout << "RES_DUMP " << tag << " step=" << Iteration::outerSteps;
            for ( int equation = 0; equation < nEquations && equation < 5; ++ equation )
                std::cout << " r" << equation << "=" << ( *res )[ equation ][ probe ];
            std::cout << std::endl;
        };
        dumpRes( "BEFORE" );
        backend.CalcAndAddReconstructedFaceFlux(
            connectivity, residual, DeviceFluxScheme( nscom.ischeme ),
            nscom.gama_ref,
            hostFluxPointer, Zone::GetUnsGrid(),
            ! this->UseHipDeviceStateUpdate() || needFlux || dumpDue,
            DeviceOwnsExplicitRkBookkeeping()
                ? ResidualBase::DeviceZero
                : ResidualBase::HostUpload,
            ctrl.centropy1, ctrl.centropy2, ctrl.ieigenfix == 1 );
        dumpRes( "AFTER" );
        if ( needFlux )
        {
            for ( int equation = 0; equation < nEquations; ++ equation )
                for ( int face = 0; face < nFaces; ++ face )
                    ( * invflux )[ equation ][ face ] =
                        faceFlux[ equation * nFaces + face ];
        }
        return;
    }

    PrimitiveFaceStateView primitiveState;
    primitiveState.nFaces = nFaces;
    primitiveState.nEquations = nEquations;
    primitiveState.primitiveLeft = nullptr;
    primitiveState.primitiveRight = nullptr;
    for ( int equation = 0; equation < nEquations; ++ equation )
    {
        primitiveState.primitiveLeftComponents[ equation ] =
            ( * limf->qf1 )[ equation ].data();
        primitiveState.primitiveRightComponents[ equation ] =
            ( * limf->qf2 )[ equation ].data();
    }
    primitiveState.xNormal = ( * ug.xfn ).data();
    primitiveState.yNormal = ( * ug.yfn ).data();
    primitiveState.zNormal = ( * ug.zfn ).data();
    primitiveState.meshVelocityNormal = ( * ug.vfn ).data();
    primitiveState.faceArea = ( * ug.farea ).data();
    primitiveState.gamma = nscom.gama_ref;

    FaceConnectivityView connectivity;
    connectivity.nFaces = nFaces;
    connectivity.nBoundaryFaces = ug.nBFaces;
    connectivity.leftCell = ug.lcf->data();
    connectivity.rightCell = ug.rcf->data();

    ResidualView residual;
    residual.nCells = nCells;
    residual.nEquations = nEquations;
    residual.values = nullptr;
    for ( int equation = 0; equation < nEquations; ++ equation )
    {
        residual.components[ equation ] = ( * res )[ equation ].data();
    }

    std::vector< Real > faceFlux;
    FaceFluxView hostFlux;
    FaceFluxView * hostFluxPointer = nullptr;
    if ( needFlux )
    {
        faceFlux.resize( nEquations * nFaces );
        hostFlux.nFaces = nFaces;
        hostFlux.nEquations = nEquations;
        hostFlux.values = faceFlux.data();
        hostFluxPointer = & hostFlux;
    }

    HipFluxBackend & backend = HipFluxBackend::Shared();
    backend.CalcAndAddPrimitiveFaceFlux(
        primitiveState, connectivity, residual,
        DeviceFluxScheme( nscom.ischeme ), hostFluxPointer, grid,
        ctrl.centropy1, ctrl.centropy2, ctrl.ieigenfix == 1 );

    if ( needFlux )
    {
        for ( int equation = 0; equation < nEquations; ++ equation )
        {
            for ( int face = 0; face < nFaces; ++ face )
            {
                ( * invflux )[ equation ][ face ] =
                    faceFlux[ equation * nFaces + face ];
            }
        }
    }

#else
    throw std::runtime_error(
        "UNs HIP batch path was called without HIP support." );
#endif
}

void UNsInvFlux::CalcInvFluxBatch( FluxBackend & backend )
{
    const int nFaces = ug.nFaces;
    const int nEquations = limiter->GetNEquations();
    std::vector< Real > primitiveLeft( nEquations * nFaces );
    std::vector< Real > primitiveRight( nEquations * nFaces );
    std::vector< Real > xNormal( nFaces );
    std::vector< Real > yNormal( nFaces );
    std::vector< Real > zNormal( nFaces );
    std::vector< Real > meshVelocityNormal( nFaces );
    std::vector< Real > faceArea( nFaces );
    std::vector< Real > faceFlux( nEquations * nFaces );

    { ScopedStageTimer timer( "host_pack" );
    MRField * qf1 = limiter->GetLeftField();
    MRField * qf2 = limiter->GetRightField();

    for ( int face = 0; face < nFaces; ++ face )
    {
        xNormal[ face ] = ( * ug.xfn )[ face ];
        yNormal[ face ] = ( * ug.yfn )[ face ];
        zNormal[ face ] = ( * ug.zfn )[ face ];
        meshVelocityNormal[ face ] = ( * ug.vfn )[ face ];
        faceArea[ face ] = ( * ug.farea )[ face ];
        for ( int equation = 0; equation < nEquations; ++ equation )
        {
            primitiveLeft[ equation * nFaces + face ] =
                ( * qf1 )[ equation ][ face ];
            primitiveRight[ equation * nFaces + face ] =
                ( * qf2 )[ equation ][ face ];
        }
    }
    }

    PrimitiveFaceStateView primitiveState;
    primitiveState.nFaces = nFaces;
    primitiveState.nEquations = nEquations;
    primitiveState.primitiveLeft = primitiveLeft.data();
    primitiveState.primitiveRight = primitiveRight.data();
    primitiveState.xNormal = xNormal.data();
    primitiveState.yNormal = yNormal.data();
    primitiveState.zNormal = zNormal.data();
    primitiveState.meshVelocityNormal = meshVelocityNormal.data();
    primitiveState.faceArea = faceArea.data();
    primitiveState.gamma = nscom.gama_ref;

    FaceFluxView flux;
    flux.nFaces = nFaces;
    flux.nEquations = nEquations;
    flux.values = faceFlux.data();

    { ScopedStageTimer timer( "flux_backend" );
        EulerCpuAdapter adapter;
        adapter.CalcInvFlux( primitiveState, flux, backend, 1 );
    }
    for ( int equation = 0; equation < nEquations; ++ equation )
    {
        for ( int face = 0; face < nFaces; ++ face )
        {
            ( * invflux )[ equation ][ face ] =
                faceFlux[ equation * nFaces + face ];
        }
    }
}

void UNsInvFlux::PrepareFaceValue()
{
    gcom.xfn   = ( * ug.xfn   )[ ug.fId ];
    gcom.yfn   = ( * ug.yfn   )[ ug.fId ];
    gcom.zfn   = ( * ug.zfn   )[ ug.fId ];
    gcom.vfn   = ( * ug.vfn   )[ ug.fId ];
    gcom.farea = ( * ug.farea )[ ug.fId ];

    nscom.gama1 = ( * unsf.gama )[ 0 ][ ug.lc ];
    nscom.gama2 = ( * unsf.gama )[ 0 ][ ug.rc ];
    nscom.gama  = half * ( nscom.gama1 + nscom.gama2 );

    inv.gama1 = nscom.gama1;
    inv.gama2 = nscom.gama2;
    inv.gama  = half * ( inv.gama1 + inv.gama2 );

    MRField * qf1 = limiter->GetLeftField();
    MRField * qf2 = limiter->GetRightField();

    int nEquations = limiter->GetNEquations();

    for ( int iEqu = 0; iEqu < nEquations; ++ iEqu )
    {
        inv.prim1[ iEqu ] = ( * qf1 )[ iEqu ][ ug.fId ];
        inv.prim2[ iEqu ] = ( * qf2 )[ iEqu ][ ug.fId ];
    }
}

void UNsInvFlux::UpdateFaceInvFlux()
{
    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        ( * invflux )[ iEqu ][ ug.fId ] = gcom.farea * inv.flux[ iEqu ];
    }
}

void UNsInvFlux::DumpInvFluxTrace()
{
    const char * configured = std::getenv( "ONEFLOW_UNS_TRACE_FILE" );
    if ( configured == nullptr || configured[ 0 ] == '\0' ) return;
    const std::string tracePath = RankAwareTracePath( configured );

    MRField * qf1 = limiter->GetLeftField();
    MRField * qf2 = limiter->GetRightField();

    if ( limiter->limfIsNullPtr() || qf1 == nullptr || qf2 == nullptr
         || invflux == nullptr )
    {
        throw std::runtime_error(
            "UNsInvFlux trace requested before face fields are available" );
    }

    std::ofstream output( tracePath, std::ios::binary | std::ios::trunc );
    if ( ! output )
    {
        throw std::runtime_error( "cannot open UNsInvFlux trace file" );
    }

    const char magic[ 8 ] = { 'O', 'F', 'T', 'R', 'C', '0', '1', '\0' };
    const std::uint64_t nFaces = static_cast< std::uint64_t >( ug.nFaces );
    const std::uint32_t nEquations =
        static_cast< std::uint32_t >( limiter->GetNEquations() );
    const std::uint32_t nArrays = 3;
    output.write( magic, sizeof( magic ) );
    output.write(
        reinterpret_cast< const char * >( & nFaces ), sizeof( nFaces ) );
    output.write(
        reinterpret_cast< const char * >( & nEquations ),
        sizeof( nEquations ) );
    output.write(
        reinterpret_cast< const char * >( & nArrays ), sizeof( nArrays ) );

    auto writeField = [&]( const MRField & field )
    {
        for ( std::uint32_t equation = 0; equation < nEquations; ++ equation )
        {
            const auto & values = field[ equation ];
            if ( values.size() < nFaces )
            {
                throw std::runtime_error(
                    "UNsInvFlux trace field has an invalid face extent" );
            }
            output.write(
                reinterpret_cast< const char * >( values.data() ),
                static_cast< std::streamsize >(
                    nFaces * sizeof( Real ) ) );
        }
    };

    writeField( *limiter->GetLeftField());
    writeField( *limiter->GetRightField());
    writeField( *invflux );
    if ( ! output )
    {
        throw std::runtime_error( "failed while writing UNsInvFlux trace" );
    }
}

void UNsInvFlux::DumpInvFluxStageTrace()
{
    const char * configured = std::getenv( "ONEFLOW_UNS_STAGE_TRACE_FILE" );
    if ( configured == nullptr || configured[ 0 ] == 0 ) return;
    const std::string tracePath = RankAwareTracePath( configured );
    if ( limf == nullptr || limf->qf1 == nullptr || limf->qf2 == nullptr
         || invflux == nullptr )
    {
        throw std::runtime_error(
            "UNs stage trace requested before face fields are available" );
    }

    UnsGrid * grid = Zone::GetUnsGrid();
    MRField * residual = GetFieldPointer< MRField >( grid, "res" );
    if ( residual == nullptr )
    {
        throw std::runtime_error(
            "UNs stage trace requested before residual is available" );
    }

    static std::uint32_t sequence = 0;
    const std::uint32_t nEquations =
        static_cast< std::uint32_t >( limf->nEqu );
    AppendStageTraceRecord(
        tracePath.c_str(), 1, sequence,
        static_cast< std::uint64_t >( ug.nFaces ), nEquations,
        { limf->qf1.get(), limf->qf2.get(), invflux.get() } );
    std::vector< std::vector< Real > > faceMetadata(
        8, std::vector< Real >( ug.nFaces ) );
    for ( int face = 0; face < ug.nFaces; ++ face )
    {
        faceMetadata[ 0 ][ face ] = ( * ug.lcf )[ face ];
        faceMetadata[ 1 ][ face ] = ( * ug.rcf )[ face ];
        faceMetadata[ 2 ][ face ] = face < ug.nBFaces ? 1.0 : 0.0;
        faceMetadata[ 3 ][ face ] = ( * ug.xfn )[ face ];
        faceMetadata[ 4 ][ face ] = ( * ug.yfn )[ face ];
        faceMetadata[ 5 ][ face ] = ( * ug.zfn )[ face ];
        faceMetadata[ 6 ][ face ] = ( * ug.vfn )[ face ];
        faceMetadata[ 7 ][ face ] = ( * ug.farea )[ face ];
    }
    AppendStageTraceValueRecord(
        tracePath.c_str(), 4, sequence,
        static_cast< std::uint64_t >( ug.nFaces ), faceMetadata );
    AppendStageTraceRecord(
        tracePath.c_str(), 2, sequence,
        static_cast< std::uint64_t >( ug.nCells ), nEquations,
        { residual } );
    ++ sequence;
}

void UNsInvFlux::AddInvFlux()
{
    ScopedStageTimer timer( "residual_state_update" );
    UnsGrid * grid = Zone::GetUnsGrid();
    MRField * res = FieldHome::GetResidualField( grid, NS_SOLVER );

    ONEFLOW::AddF2CField( res, invflux.get() );
}

void UNsInvFlux::Alloc()
{
    const bool needsAllocation =
        invflux == nullptr
        || invflux->GetNEqu() != nscom.nEqu
        || nscom.nEqu == 0
        || ( nscom.nEqu > 0
             && (*invflux)[ 0 ].size()
                    != static_cast< HXSize_t >( ug.nFaces ) );
    if ( needsAllocation )
    {
        invflux = std::make_unique<MRField>( nscom.nEqu, ug.nFaces );
    }
}

void UNsInvFlux::DeAlloc()
{
    // Keep the host face-flux buffer for reuse by the next RK stage.
}

EndNameSpace
