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

#include "UNsUpdate.h"
#include "UCom.h"
#include "NsCom.h"
#include "UNsCom.h"
#include "NsIdx.h"
#include "Zone.h"
#include "ZoneState.h"
#include "HXMath.h"
#include "Iteration.h"
#include "GridState.h"
#include "SolverDef.h"
#include "AccelViews.h"
#include "AccelRuntime.h"
#include "HipProductionPolicy.h"
#include "Ctrl.h"
#include "NsRhs.h"
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

void AppendStateStageTrace(
    const char * path,
    std::uint32_t sequence,
    std::uint64_t nItems,
    std::uint32_t nEquations,
    const MRField & field )
{
    std::ofstream output( path, std::ios::binary | std::ios::app );
    if ( ! output )
    {
        throw std::runtime_error( "cannot open UNs stage trace file" );
    }

    const char magic[ 8 ] = "OFSTG01";
    const std::uint32_t kind = 3;
    const std::int32_t outerStep = Iteration::outerSteps;
    const std::int32_t gridLevel = GridState::gridLevel;
    const std::uint32_t nArrays = 1;
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

    if ( field.GetNEqu() < nEquations )
    {
        throw std::runtime_error(
            "UNs state trace field has an invalid equation extent" );
    }
    for ( std::uint32_t equation = 0;
          equation < nEquations; ++ equation )
    {
        const auto & values = field[ equation ];
        if ( values.size() < nItems )
        {
            throw std::runtime_error(
                "UNs state trace field has an invalid cell extent" );
        }
        output.write(
            reinterpret_cast< const char * >( values.data() ),
            static_cast< std::streamsize >( nItems * sizeof( Real ) ) );
    }
    if ( ! output )
    {
        throw std::runtime_error( "failed while writing UNs state trace" );
    }
}

}

UNsUpdate::UNsUpdate()
{
}

UNsUpdate::~UNsUpdate()
{
}

namespace
{

bool UseHipDeviceStateUpdate( int solverType )
{
    const char * feature = "ONEFLOW_ENABLE_UNS_HIP_STATE_UPDATE";
    const bool explicitRequest =
        HipProductionFeatureExplicitlyEnabled( feature );
    if ( ! HipProductionFeatureEnabled( feature ) ) return false;
    if ( solverType != NS_SOLVER )
    {
        if ( ! explicitRequest ) return false;
        throw std::runtime_error(
            "HIP state update is only supported for the NS solver" );
    }
    if ( ! HipProductionFeatureEnabled(
             "ONEFLOW_ENABLE_UNS_HIP_RECONSTRUCTION" ) )
    {
        if ( ! explicitRequest ) return false;
        throw std::runtime_error(
            "HIP state update requires HIP device reconstruction" );
    }
    if ( nscom.chemModel != 0 || nscom.nTModel != 1 )
    {
        if ( ! explicitRequest ) return false;
        throw std::runtime_error(
            "HIP state update requires a single ideal-gas temperature field" );
    }
#ifndef ONEFLOW_ENABLE_HIP
    throw std::runtime_error(
        "HIP state update was requested without HIP support" );
#else
    return true;
#endif
}

bool RetainHipStateForBoundary()
{
    return HipProductionFeatureEnabled( "ONEFLOW_ENABLE_UNS_HIP_BOUNDARY" )
        || HipProductionFeatureEnabled( "ONEFLOW_ENABLE_UNS_HIP_INTERFACE" );
}

// The device state may own the newest q for a whole RK step without a host
// mirror.  This is only valid for the steady explicit device chain: the
// dual-time path shifts the host history levels every step and the CPU oracle
// must stay authoritative, so it is excluded.
bool UseDeviceAuthoritativeHipState()
{
    return HipDeviceAuthoritativeStateEnabled()
        && HipProductionBackendSelected()
        && ctrl.startStrategy != 2;
}

}

void UNsUpdate::UpdateFlowField( int solverType )
{
    if ( UseHipDeviceStateUpdate( solverType ) )
    {
#ifdef ONEFLOW_ENABLE_HIP
        UnsGrid * grid = Zone::GetUnsGrid();
        ug.Init();
        unsf.Init();

        CellStateUpdateView view;
        view.nCells = ug.nCells;
        view.nGhostCells = ug.nTCell - ug.nCells;
        view.nEquations = nscom.nEqu;
        view.timeStep = ( * unsf.timestep )[ 0 ].data();
        view.cellVolume = ug.cvol->data();
        view.gamma = nscom.gama_ref;
        view.rkCoefficient = 1.0;
        view.cacheKey = grid;
        view.retainResidentState = RetainHipStateForBoundary();
        const bool deviceAuthoritative =
            view.retainResidentState && UseDeviceAuthoritativeHipState();
        view.deviceAuthoritative = deviceAuthoritative;
        for ( int equation = 0; equation < view.nEquations; ++ equation )
        {
            view.primitive[ equation ] = ( * unsf.q )[ equation ].data();
        }
        // With a device-authoritative state the primitive field stays on the
        // device and the host mirror is refreshed on demand; the 5-equation
        // download per stage was the single largest fixed cost.
        HipFluxBackend::Shared().UpdatePrimitiveState(
            view, ! deviceAuthoritative );
        if ( ! deviceAuthoritative )
        {
            for ( int cId = 0; cId < ug.nCells; ++ cId )
            {
                const Real density = ( * unsf.q )[ IDX::IR ][ cId ];
                const Real pressure = ( * unsf.q )[ IDX::IP ][ cId ];
                ( * unsf.tempr )[ 0 ][ cId ] =
                    pressure / ( nscom.statecoef * density );
            }
        }
        this->DumpUpdatedStateTrace();
        return;
#endif
    }

    GetUpdateField( solverType, this->q, this->dq );

    ug.Init();
    unsf.Init();

    for ( int cId = 0; cId < ug.nCells; ++ cId )
    {
        ug.cId = cId;

        this->PrepareData();

        this->CalcFlowField();

        this->UpdateFlowFieldValue();
    }

    // The host RK update rewrote owned cells of the host q mirror.  Publish
    // that mutation through the resident-state generation contract so the next
    // device consumer cannot reuse a stale resident primitive state.
#ifdef ONEFLOW_ENABLE_HIP
    if ( HipProductionBackendSelected() )
    {
        // A host write is only valid when the mirror is current; refresh it
        // first if the device owned the field.
        NsRefreshHostStateMirror();
        HipFluxBackend::Shared().NotifyHostFieldStateModified(
            Zone::GetUnsGrid() );
    }
#endif

    this->DumpUpdatedStateTrace();
}

void UNsUpdate::DumpUpdatedStateTrace()
{
    const char * configured = std::getenv( "ONEFLOW_UNS_STAGE_TRACE_FILE" );
    if ( configured == nullptr || configured[ 0 ] == 0 ) return;
#ifdef ONEFLOW_ENABLE_HIP
    // The trace reads the host mirror; pull the complete field back first when
    // the device owns it.
    NsRefreshHostStateMirror();
#endif
    const std::string tracePath = RankAwareTracePath( configured );
    if ( unsf.q == nullptr )
    {
        throw std::runtime_error(
            "UNs state trace requested before flow state is available" );
    }

    static std::uint32_t sequence = 0;
    AppendStateStageTrace(
        tracePath.c_str(), sequence,
        static_cast< std::uint64_t >( ug.nCells ),
        static_cast< std::uint32_t >( nscom.nTEqu ), * unsf.q );
    ++ sequence;
}

void UNsUpdate::PrepareData()
{
    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.prim [ iEqu ] = ( * unsf.q )[ iEqu ][ ug.cId ];
        nscom.prim0[ iEqu ] = ( * unsf.q )[ iEqu ][ ug.cId ];
    }

    for ( int iEqu = 0; iEqu < nscom.nTModel; ++ iEqu )
    {
        nscom.t [ iEqu ] = ( * unsf.tempr )[ iEqu ][ ug.cId ];
        nscom.t0[ iEqu ] = ( * unsf.tempr )[ iEqu ][ ug.cId ];
    }

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.dq[ iEqu ] = ( * unsf.dq )[ iEqu ][ ug.cId ];
    }

    nscom.gama = ( * unsf.gama )[ 0 ][ ug.cId ];
}

void UNsUpdate::DumpProbeInfo()
{
    std::cout << std::setprecision( 17 );
    std::cout << "Warning : non-physical state"
              << " step = " << Iteration::outerSteps
              << ", zid = " << ZoneState::zid
              << ", cid = " << ug.cId
              << ", density = " << nscom.prim[ IDX::IR ]
              << ", pressure = " << nscom.prim[ IDX::IP ]
              << ", density0 = " << nscom.prim0[ IDX::IR ]
              << ", pressure0 = " << nscom.prim0[ IDX::IP ]
              << ", dq_density = " << nscom.dq[ IDX::IR ]
              << ", dq_pressure = " << nscom.dq[ IDX::IP ]
              << ", timestep = " << nscom.timestep
              << ", dt = " << nscom.dt
              << std::endl;
}

void UNsUpdate::SolutionFix()
{
    nscom.prim = 0;

    Real sumV = 0.0;

    int fn = ( * ug.c2f )[ ug.cId ].size();

    for ( int iFace = 0; iFace < fn; ++ iFace )
    {
        int fId = ( * ug.c2f )[ ug.cId ][ iFace ];

        ug.lc = ( * ug.lcf )[ fId ];
        ug.rc = ( * ug.rcf )[ fId ];

        int iNei = ug.lc;
        if ( ug.cId == ug.lc  ) iNei = ug.rc;

        Real volN = one;

        sumV += volN;

        for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
        {
            Real f = ( * unsf.q )[ iEqu ][ iNei ];
            if ( iEqu == IDX::IR || iEqu == IDX::IP )
            {
                f = ABS( f );
            }
            nscom.prim[ iEqu ] += f * volN;
        }
    }

    Real rVol = 1.0 / sumV;

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.prim[ iEqu ] *= rVol;
    }
}

void UNsUpdate::UpdateFlowFieldValue()
{
    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        ( * unsf.q )[ iEqu ][ ug.cId ] = nscom.prim[ iEqu ];
    }

    for ( int iEqu = 0; iEqu < nscom.nTModel; ++ iEqu )
    {
        ( * unsf.tempr )[ iEqu ][ ug.cId ] = nscom.t[ iEqu ];
    }
}

EndNameSpace
