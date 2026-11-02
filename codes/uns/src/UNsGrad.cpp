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

#include "UNsGrad.h"
#include "UCom.h"
#include "NsCom.h"
#include "DataBase.h"
#include "FieldWrap.h"
#include "SolverDef.h"
#include "FieldManager.h"
#include "FaceTopo.h"
#include "BcRecord.h"
#include "UnsGrid.h"
#include "Zone.h"
#include "FaceMesh.h"
#include "CellMesh.h"
#include "CellTopo.h"
#include "HipProductionPolicy.h"
#include "Ctrl.h"
#include "StageProfiler.h"
#include "InterFace.h"
#include "InterField.h"
#ifdef ONEFLOW_ENABLE_HIP
#include "HipFluxBackend.h"
#include <cstdlib>
#endif
#include <stdexcept>
#include <algorithm>
#include <vector>

BeginNameSpace( ONEFLOW )

UNsGrad uns_grad;
UTGrad ut_grad;

UNsGrad::UNsGrad()
{
    ;
}

UNsGrad::~UNsGrad()
{
    ;
}

void UNsGrad::Init()
{
    UnsGrid * grid = Zone::GetUnsGrid();

    // Pin solver identity (do not rely on ambient SolverState).
    const int sid = NS_SOLVER;
    name  = FieldHome::GetGradName( sid, 0 );
    namex = FieldHome::GetGradName( sid, 1 );
    namey = FieldHome::GetGradName( sid, 2 );
    namez = FieldHome::GetGradName( sid, 3 );

    q    = FieldHome::GetSolutionField( grid, sid );
    dqdx = FieldHome::GetGradField( grid, sid, 1 );
    dqdy = FieldHome::GetGradField( grid, sid, 2 );
    dqdz = FieldHome::GetGradField( grid, sid, 3 );
    bdqdx = FieldHome::GetBcGradField( grid, sid, 0 );
    bdqdy = FieldHome::GetBcGradField( grid, sid, 1 );
    bdqdz = FieldHome::GetBcGradField( grid, sid, 2 );

    this->nEqu = nscom.nTEqu;

    this->istore = 1;
}

void UNsGrad::CalcGradHip( bool downloadToHost, bool packedExchange,
                          bool mirrorToHost )
{
#ifdef ONEFLOW_ENABLE_HIP
    if ( q == nullptr || dqdx == nullptr || dqdy == nullptr || dqdz == nullptr )
    {
        throw std::runtime_error( "UNs HIP gradient fields are unavailable" );
    }
    UnsGrid * grid = Zone::GetUnsGrid();
    CellGradientView view;
    view.nCells = ug.nCells;
    view.nGhostCells = ug.nTCell - ug.nCells;
    view.nFaces = ug.nFaces;
    view.nBoundaryFaces = ug.nBFaces;
    view.nEquations = nEqu;
    view.xFace = ug.xfc->data();
    view.yFace = ug.yfc->data();
    view.zFace = ug.zfc->data();
    view.xNormal = ug.xfn->data();
    view.yNormal = ug.yfn->data();
    view.zNormal = ug.zfn->data();
    view.faceArea = ug.farea->data();
    view.xCell = ug.xcc->data();
    view.yCell = ug.ycc->data();
    view.zCell = ug.zcc->data();
    view.cellVolume = ug.cvol->data();
    view.leftCell = ug.lcf->data();
    view.rightCell = ug.rcf->data();
    view.cacheKey = grid;
    view.downloadToHost = downloadToHost;
    for ( int equation = 0; equation < nEqu; ++ equation )
    {
        view.q[ equation ] = ( * q )[ equation ].data();
        view.dqdx[ equation ] = ( * dqdx )[ equation ].data();
        view.dqdy[ equation ] = ( * dqdy )[ equation ].data();
        view.dqdz[ equation ] = ( * dqdz )[ equation ].data();
    }
    HipFluxBackend::Shared().CalcGradient( view );
    if ( packedExchange && ! view.downloadToHost )
    {
        // The host mirror is consumed by Grad::SwapBcGrad only at the
        // boundary-adjacent owned cells: UploadInterfaceValue packs
        // dqdx[leftCell] for interface faces and StoreBcGrad stores the same
        // value for every boundary face.  Mirror exactly that set, then run the
        // unchanged host swap so the receive buffer (and therefore the
        // one-round lag) keeps the CPU oracle semantics.
        if ( mirrorToHost )
        {
            if ( ug.bcRecord == nullptr )
            {
                throw std::runtime_error(
                    "UNs HIP gradient mirror requires a boundary record" );
            }
            std::vector< int > cells( static_cast< std::size_t >( ug.nBFaces ) );
            for ( int face = 0; face < ug.nBFaces; ++ face )
                cells[ face ] = ( * ug.lcf )[ face ];
            std::sort( cells.begin(), cells.end() );
            cells.erase( std::unique( cells.begin(), cells.end() ),
                         cells.end() );
            Real * dx[ 5 ] = {};
            Real * dy[ 5 ] = {};
            Real * dz[ 5 ] = {};
            for ( int equation = 0; equation < nEqu; ++ equation )
            {
                dx[ equation ] = ( * dqdx )[ equation ].data();
                dy[ equation ] = ( * dqdy )[ equation ].data();
                dz[ equation ] = ( * dqdz )[ equation ].data();
            }
            HipFluxBackend::Shared().MirrorGradientCellsToHost(
                grid, nEqu, cells, dx, dy, dz );
        }
        if ( HipPackedGradientExchangeEnabled()
             && grid->interFace != nullptr )
        {
            const bool deviceUnpack =
                HipDeviceGradientInterfaceEnabled( vis_model.vismodel != 0 );
            // Only the receive side runs here.  The send storage is packed
            // once per outer iteration by UploadInterfaceData, right before
            // the round trip that exchanges it.
            InterFace * interFace = grid->interFace.get();
            MRField * fields[ 3 ] = { dqdx, dqdy, dqdz };
            const std::string names[ 3 ] = { namex, namey, namez };
            std::vector< int > endpointCells;
            std::vector< Real > packed;
            for ( int ghostId = MAX_GHOST_LEVELS - 1;
                  ghostId >= 0; -- ghostId )
            {
                const int nIFaces = interFace->nIFaces;
                if ( nIFaces <= 0 ) continue;
                endpointCells.resize( static_cast< std::size_t >( nIFaces ) );
                packed.resize(
                    static_cast< std::size_t >( 3 ) * nEqu * nIFaces );
                {
                    ScopedStageTimer timer( "swap_download" );
                    if ( ! deviceUnpack )
                    {
                        DownloadInterfaceValuePacked(
                            grid, fields, names, 3, nEqu, ghostId,
                            packed.data() );
                    }
                    for ( int iFace = 0; iFace < nIFaces; ++ iFace )
                    {
                        int iCell = -1;
                        grid->GetFaceTopo().GetTId(
                            iFace, ghostId + 1, iCell );
                        endpointCells[ iFace ] = iCell;
                    }
                    if ( deviceUnpack )
                    {
                        HipFluxBackend::Shared().UnpackGradientFacesFromDevice(
                            grid, endpointCells.data(), nIFaces, nEqu );
                    }
                    else
                    {
                        HipFluxBackend::Shared().UnpackGradientFacesFromHost(
                            grid, endpointCells.data(), nIFaces, nEqu,
                            packed.data() );
                    }
                }
            }
        }
        else
        {
            this->SwapBcGrad();
        }
    }
    else if ( view.downloadToHost )
    {
        this->SwapBcGrad();
    }
#else
    throw std::runtime_error(
        "UNs HIP gradient was called without HIP support" );
#endif
}

UTGrad::UTGrad()
{
    ;
}

UTGrad::~UTGrad()
{
    ;
}

void UTGrad::Init()
{
    UnsGrid * grid = Zone::GetUnsGrid();

    // Shared NS temperature field + gradient (not solver-prefixed; not in SolverInfo).
    name  = "tempr";
    namex = "dtdx";
    namey = "dtdy";
    namez = "dtdz";

    q    = GetFieldPointer< MRField > ( grid, name  );
    dqdx = GetFieldPointer< MRField > ( grid, namex );
    dqdy = GetFieldPointer< MRField > ( grid, namey );
    dqdz = GetFieldPointer< MRField > ( grid, namez );

    this->nEqu = nscom.nTModel;

    this->istore = 0;
}

EndNameSpace
