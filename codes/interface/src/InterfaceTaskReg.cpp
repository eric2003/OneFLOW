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

#include "InterfaceTaskReg.h"
#include <memory>
#include "InterField.h"
#include "Fatal.h"
#include "ActionState.h"
#include "HXMath.h"
#include "SolverDef.h"
#include "SolverInfo.h"
#include "Restart.h"
#include "Unsteady.h"
#include "UnsteadyImp.h"
#include "Update.h"
#include "FieldWrap.h"
#include "CmxTask.h"
#include "DataBase.h"
#include "DataBook.h"
#include "Lusgs.h"
#include "Lhs.h"
#include "FieldManager.h"
#include "StageProfiler.h"
#include "HipProductionPolicy.h"
#include "AccelRuntime.h"
#include "HipFluxBackend.h"
#include "FaceTopo.h"
#include "Ctrl.h"
#include "SolverState.h"
#include "Zone.h"
#include "Grid.h"
#include "GridState.h"
#include "UnsGrid.h"
#include "InterFace.h"
#include "RegisterUtils.h"
#include "FieldRecord.h"
#include "UResidual.h"
#include "UNsCom.h"
#include "TaskRegister.h"

BeginNameSpace( ONEFLOW )

REGISTER_TASK( RegisterInterfaceTask )

void RegisterInterfaceTask()
{
    REGISTER_DATA_CLASS( CalcInterfaceGrad );
    REGISTER_DATA_CLASS( UploadInterfaceData );
    REGISTER_DATA_CLASS( DownloadInterfaceData );
    REGISTER_DATA_CLASS( PrepareInterfaceField );
    REGISTER_DATA_CLASS( PrepareOversetInterfaceField );
    REGISTER_DATA_CLASS( ExchangeInterfaceDq );
}


void CalcInterfaceGrad( StringField & data )
{
    ;
}

//void UploadInterfaceData( StringField & data )
//{
//    int solverType = SolverState::solverType;
//
//    FieldManager * fieldManager = FieldManagerRegistry::GetFieldManager( solverType );
//
//    fieldManager->GetInterfaceFieldProperty().UploadInterfaceValue();
//
//}

//void DownloadInterfaceData( StringField & data )
//{
//    int solverType = SolverState::solverType;
//
//    FieldManager * fieldManager = FieldManagerRegistry::GetFieldManager( solverType );
//    fieldManager->GetInterfaceFieldProperty().DownloadInterfaceValue();
//}

namespace {

// Only the gradient is device-owned: UNsGrad already packed it into the same
// send storage from the resident field, and the receive side unpacks it
// directly into the device ghosts.  Re-packing it from the host mirror would
// overwrite the device payload, which is exactly what made the host mirror
// mandatory on every step.  The other interface fields, dataBook layout and
// message sizes stay untouched, so peers remain in step.
bool SkipHostFieldPackFor( const std::string & fieldName )
{
#ifdef ONEFLOW_ENABLE_HIP
    if ( vis_model.vismodel != 0 ) return false;
    if ( fieldName == "dqdx" || fieldName == "dqdy" || fieldName == "dqdz" )
    {
        return HipDeviceGradientInterfaceEnabled( false );
    }
    // When the device owns q, the host q/tempr/gama ghosts are stale and the
    // exchange itself runs through the resident device pack.  Packing and
    // downloading stale host copies is pure cost (and the host q pack used to
    // cost several milliseconds per step).
    if ( HipSkipHostInterfaceExchange()
         && HipDeviceAuthoritativeStateEnabled()
         && ( fieldName == "q" || fieldName == "tempr"
              || fieldName == "gama" ) )
    {
        return true;
    }
    if ( ! HipUnusedInterfaceFieldsSkipped() ) return false;
    return fieldName == "limiter" || fieldName == "dq"
        || fieldName == "dtdx" || fieldName == "dtdy" || fieldName == "dtdz"
        || fieldName == "heatflux"
        || fieldName == "visl" || fieldName == "vist";
#else
    (void) fieldName;
    return false;
#endif
}

// Fill the interface send storage for the device-owned gradient.  This runs
// once per outer iteration, from UploadInterfaceData, because that is the round
// trip that actually exchanges the storage.  Packing inside every RK stage only
// repeated the same work (nStages-1 of the packs were overwritten unused) and
// put a device->host copy on the RK critical path.  The receive side still
// unpacks every stage, so the CPU one-round lag is unchanged.
void PackResidentGradientToInterface( UnsGrid * grid )
{
#ifdef ONEFLOW_ENABLE_HIP
    if ( grid == nullptr || grid->interFace == nullptr ) return;
    MRField * dqdx = ONEFLOW::GetFieldPointer< MRField >( grid, "dqdx" );
    MRField * dqdy = ONEFLOW::GetFieldPointer< MRField >( grid, "dqdy" );
    MRField * dqdz = ONEFLOW::GetFieldPointer< MRField >( grid, "dqdz" );
    if ( dqdx == nullptr || dqdy == nullptr || dqdz == nullptr ) return;
    const int nEqu = static_cast< int >( dqdx->GetNEqu() );
    if ( nEqu <= 0 ) return;

    InterFace * interFace = grid->interFace.get();
    MRField * fields[ 3 ] = { dqdx, dqdy, dqdz };
    const std::string names[ 3 ] = { "dqdx", "dqdy", "dqdz" };
    std::vector< int > endpointCells;
    std::vector< Real > packed;
    for ( int ghostId = MAX_GHOST_LEVELS - 1; ghostId >= 0; -- ghostId )
    {
        const int nIFaces = interFace->nIFaces;
        if ( nIFaces <= 0 ) continue;
        endpointCells.resize( static_cast< std::size_t >( nIFaces ) );
        packed.resize( static_cast< std::size_t >( 3 ) * nEqu * nIFaces );
        for ( int iFace = 0; iFace < nIFaces; ++ iFace )
        {
            int iCell = -1;
            grid->GetFaceTopo().GetSId( iFace, ghostId + 1, iCell );
            endpointCells[ iFace ] = iCell;
        }
        HipFluxBackend::Shared().PackGradientFacesToHost(
            grid, endpointCells.data(), nIFaces, nEqu,
            packed.data() );
        UploadInterfaceValuePacked(
            grid, fields, names, 3, nEqu, ghostId,
            packed.data() );
    }
#else
    ( void ) grid;
#endif
}

}

void UploadInterfaceData( StringField & data )
{
    Grid & baseGrid = Zone::GetGridReference();
    if ( ! ONEFLOW::IsUnsGrid( baseGrid.type ) )
    {
        return;
    }

    UnsGrid * grid = ONEFLOW::UnsGridCast( &baseGrid );

    int solverType = SolverState::solverType;
    FieldManager * fieldManager =
        FieldManagerRegistry::GetFieldManager( solverType );

    const auto & fieldData =
        fieldManager->GetInterfaceFieldProperty().GetData();

    ScopedStageTimer uploadTimer( "post_upload" );
#ifdef ONEFLOW_ENABLE_HIP
    if ( SkipHostFieldPackFor( "dqdx" ) )
    {
        ScopedStageTimer packTimer( "swap_upload" );
        PackResidentGradientToInterface( grid );
    }
#endif
    for ( const auto & [ fieldName, nEqu ] : fieldData )
    {
        if ( SkipHostFieldPackFor( fieldName ) ) continue;
        MRField * field =
            ONEFLOW::GetFieldPointer< MRField >(
                grid,
                fieldName );

        if ( field == nullptr )
        {
            Fatal(
                "Grid field is not allocated for interface upload: "
                + fieldName );
        }

        ONEFLOW::UploadInterfaceValue(
            grid,
            field,
            fieldName,
            nEqu );
    }
}

void DownloadInterfaceData( StringField & data )
{
    Grid & baseGrid = Zone::GetGridReference();
    if ( ! ONEFLOW::IsUnsGrid( baseGrid.type ) )
    {
        return;
    }

    UnsGrid * grid = ONEFLOW::UnsGridCast( &baseGrid );

    int solverType = SolverState::solverType;
    FieldManager * fieldManager =
        FieldManagerRegistry::GetFieldManager( solverType );

    const auto & fieldData =
        fieldManager->GetInterfaceFieldProperty().GetData();

    ScopedStageTimer downloadTimer( "post_download" );
#ifdef ONEFLOW_ENABLE_HIP
    if ( SkipHostFieldPackFor( "dqdx" ) )
    {
        // The receive payload only changes once per exchange round.  Upload it
        // here so every RK stage can unpack it device-to-device instead of
        // repeating the same host->device copy three times.
        MRField * dx = ONEFLOW::GetFieldPointer< MRField >( grid, "dqdx" );
        MRField * dy = ONEFLOW::GetFieldPointer< MRField >( grid, "dqdy" );
        MRField * dz = ONEFLOW::GetFieldPointer< MRField >( grid, "dqdz" );
        if ( dx != nullptr && dy != nullptr && dz != nullptr )
        {
            MRField * fields[ 3 ] = { dx, dy, dz };
            const std::string names[ 3 ] = { "dqdx", "dqdy", "dqdz" };
            const int nEqu = static_cast< int >( dx->GetNEqu() );
            std::vector< Real > packed;
            InterFace * interFace = grid->interFace.get();
            for ( int ghostId = MAX_GHOST_LEVELS - 1; ghostId >= 0; -- ghostId )
            {
                const int nIFaces = interFace->nIFaces;
                if ( nIFaces <= 0 ) continue;
                packed.resize( static_cast< std::size_t >( 3 ) * nEqu * nIFaces );
                DownloadInterfaceValuePacked(
                    grid, fields, names, 3, nEqu, ghostId, packed.data() );
                HipFluxBackend::Shared().UploadGradientFacesToDevice(
                    grid, nIFaces, nEqu, packed.data() );
            }
        }
    }
#endif
    for ( const auto & [ fieldName, nEqu ] : fieldData )
    {
        if ( SkipHostFieldPackFor( fieldName ) ) continue;
        MRField * field =
            ONEFLOW::GetFieldPointer< MRField >(
                grid,
                fieldName );

        if ( field == nullptr )
        {
            Fatal(
                "Grid field is not allocated for interface download: "
                + fieldName );
        }

        ONEFLOW::DownloadInterfaceValue(
            grid,
            field,
            fieldName,
            nEqu );
    }
}

void PrepareInterfaceField( StringField & data )
{
    Grid & grid = Zone::GetGridReference();
    InterFace * interFace = grid.interFace.get();
    if ( ! ONEFLOW::IsValid( interFace ) ) return;

    int solverType = SolverState::solverType;
    int iFk  = ( * interfaceMap )[ data[ 0 ] ];
    int iSr  = ( * sendRecvMap )[ data[ 1 ] ];

    auto fieldRecord = std::make_unique<FieldRecord>();

    PrepareInterfaceFieldRecord( solverType, iFk, iSr, fieldRecord.get() );

    SetInterfaceFieldData( iSr, fieldRecord.get() );
}

void PrepareOversetInterfaceField( StringField & data )
{
    ;
}

void ExchangeInterfaceDq( StringField & data )
{
    ;
}

EndNameSpace
