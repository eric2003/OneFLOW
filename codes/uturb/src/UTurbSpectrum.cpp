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

#include "UTurbSpectrum.h"
#include "NsCom.h"
#include "TurbCom.h"
#include "TurbLusgs.h"
#include "UTurbCom.h"
#include "UNsCom.h"
#include "UCom.h"
#include "Zone.h"
#include "Grid.h"
#include "Fatal.h"
#include "UnsGrid.h"
#include "FaceTopo.h"
#include "FaceMesh.h"
#include "CellMesh.h"
#include "HXMath.h"
#include "DataBase.h"
#include "FieldBase.h"
#include "Unsteady.h"
#include "TimeIntegration.h"
#include "Ctrl.h"
#include "NsIdx.h"
#include "HipProductionPolicy.h"
#ifdef ONEFLOW_ENABLE_HIP
#include "HipTurbulenceBackend.h"
#endif
#include <iostream>


BeginNameSpace( ONEFLOW )


UTurbSpectrum::UTurbSpectrum()
{
}

UTurbSpectrum::~UTurbSpectrum()
{
}

void UTurbSpectrum::CalcSpectrum()
{
    ug.Init();
    turblu.Init();
    turbsp.Init();
    uturbf.Init();
    this->ZeroSpectrum();
    if ( turbcom.nEqu == 1 )
    {
        this->CalcSpectrum1Equ();
    }
    else if ( turbcom.nEqu >= 2 )
    {
        this->CalcSpectrum2Equ();
    }
}

void UTurbSpectrum::ZeroSpectrum()
{
    ( * uturbf.impsr ) = 0;
    ( * uturbf.matrix_l ) = 0;
    ( * uturbf.matrix_r ) = 0;
}

void UTurbSpectrum::CalcSpectrum1Equ()
{
    if ( this->CalcSpectrum1EquHip() ) return;

    this->CalcUnsteadySpectrum();

    for ( int fId = 0; fId < ug.nFaces; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        this->PrepareFaceValue1Equ();

        this->CalcFaceSpectrum1Equ();

        this->UpdateSpectrumRadius();
    }
}

bool UTurbSpectrum::CalcSpectrum1EquHip()
{
#ifdef ONEFLOW_ENABLE_HIP
    if ( turbcom.nEqu != 1 ) return false;
    if ( ! HipProductionFeatureEnabled( "ONEFLOW_ENABLE_UTURB_HIP_SPECTRUM" ) )
    {
        return false;
    }
    if ( ctrl.idualtime != 0 )
    {
        // The device diagonal reproduces the idualtime == 0 branch only; the
        // dual-time variant adds a sp2/pdt1 term the view does not carry.
        return false;
    }
    UnsGrid * grid = Zone::GetUnsGrid();
    if ( grid == nullptr ) return false;
    if ( uturbf.impsr == nullptr || uturbf.matrix_l == nullptr
         || uturbf.matrix_r == nullptr || uturbf.q == nullptr
         || uturbf.visl == nullptr || uturbf.timestep == nullptr
         || uturbf.meanFlow.q == nullptr )
    {
        return false;
    }

    TurbSaSpectrumView view;
    view.nFaces = ug.nFaces;
    view.nBoundaryFaces = ug.nBFaces;
    view.nCells = ug.nCells;
    view.nGhostCells = ug.nTCell - ug.nCells;
    view.xNormal = ug.xfn->data();
    view.yNormal = ug.yfn->data();
    view.zNormal = ug.zfn->data();
    view.meshVelocityNormal = ug.vfn->data();
    view.faceArea = ug.farea->data();
    view.leftCell = ug.lcf->data();
    view.rightCell = ug.rcf->data();
    view.density = ( * uturbf.meanFlow.q )[ IDX::IR ].data();
    view.velocity[ 0 ] = ( * uturbf.meanFlow.q )[ IDX::IU ].data();
    view.velocity[ 1 ] = ( * uturbf.meanFlow.q )[ IDX::IV ].data();
    view.velocity[ 2 ] = ( * uturbf.meanFlow.q )[ IDX::IW ].data();
    view.workingVariable = ( * uturbf.q )[ ISA ].data();
    view.laminarViscosity = ( * uturbf.visl )[ 0 ].data();
    view.cellVolume = ug.cvol->data();
    view.unsteadyDiagonal = true;
    view.timestep = ( * uturbf.timestep )[ 0 ].data();
    view.turbCflRatio = turbcom.turb_cfl_ratio;
    view.unsteadySp2 = timeIntegration.sp2;
    view.invReynolds = turbcom.oreynolds;
    view.invSigma = turbcom.osigma;
    view.implicitSpectrum = ( * uturbf.impsr )[ 0 ].data();
    view.matrixLeft = ( * uturbf.matrix_l )[ 0 ].data();
    view.matrixRight = ( * uturbf.matrix_r )[ 0 ].data();
    view.cacheKey = grid;

    // The host path starts from a zeroed spectrum; the backend does the same
    // internally before accumulating.
    this->ZeroSpectrum();
    HipTurbulenceBackend::Shared().CalcSpectrum( view );
    return true;
#else
    return false;
#endif
}

void UTurbSpectrum::CalcSpectrum2Equ()
{
    this->CalcUnsteadySpectrum();

    for ( int fId = 0; fId < ug.nFaces; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        this->PrepareFaceValue2Equ();

        this->CalcFaceSpectrum2Equ();

        this->UpdateSpectrumRadius();
    }
}

void UTurbSpectrum::CalcUnsteadySpectrum()
{
    if ( ctrl.idualtime == 0 )//Single time step, note: Yes timeIntegration.sp2 !
    {
        for ( int cId = 0; cId < ug.nCells; ++ cId )
        {
            Real vol = ( * ug.cvol )[ cId ];
            Real ts  = ( * uturbf.timestep )[ 0 ][ cId ] * turbcom.turb_cfl_ratio;
            Real unsteadyTerm = (  timeIntegration.sp2 / ts ) * vol;

            for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
            {
                ( * uturbf.impsr )[ iEqu ][ cId ] += unsteadyTerm;
            }
        }
    }
    else
    {
        for ( int cId = 0; cId < ug.nCells; ++ cId )
        {
            Real vol = ( * ug.cvol )[ cId ];
            Real ts  = ( * unsf.timestep )[ 0 ][ cId ] * turbcom.turb_cfl_ratio;
            Real unsteadyTerm = (  timeIntegration.sp1 / ts + timeIntegration.sp2 / ctrl.pdt1 ) * vol;

            for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
            {
                ( * uturbf.impsr )[ iEqu ][ cId ] += unsteadyTerm;
            }

        }
    }
}


void UTurbSpectrum::UpdateSpectrumRadius()
{
    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        ( * uturbf.impsr )[ iEqu ][ ug.lc ] += turbsp.radius1[ iEqu ];
        ( * uturbf.impsr )[ iEqu ][ ug.rc ] += turbsp.radius2[ iEqu ];
    }

    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        ( * uturbf.matrix_l )[ iEqu ][ ug.fId ] += turbsp.matrix1[ iEqu ];
        ( * uturbf.matrix_r )[ iEqu ][ ug.fId ] += turbsp.matrix2[ iEqu ];
    }
}

void UTurbSpectrum::PrepareFaceValue1Equ()
{
    gcom.xfn   = ( * ug.xfn   )[ ug.fId ];
    gcom.yfn   = ( * ug.yfn   )[ ug.fId ];
    gcom.zfn   = ( * ug.zfn   )[ ug.fId ];
    gcom.vfn   = ( * ug.vfn   )[ ug.fId ];
    gcom.farea = ( * ug.farea )[ ug.fId ];

    gcom.cvol1 = ( * ug.cvol )[ ug.lc ];
    gcom.cvol2 = ( * ug.cvol )[ ug.rc ];

    for ( int iEqu = 0; iEqu < nscom.nEqu; ++ iEqu )
    {
        turbsp.q1_ns[ iEqu ] = ( * uturbf.meanFlow.q )[ iEqu ][ ug.lc ];
        turbsp.q2_ns[ iEqu ] = ( * uturbf.meanFlow.q )[ iEqu ][ ug.rc ];
    }

    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turbsp.q1[ iEqu ] = ( * uturbf.q )[ iEqu ][ ug.lc ];
        turbsp.q2[ iEqu ] = ( * uturbf.q )[ iEqu ][ ug.rc ];
    }

    turbcom.visl1 = ( * uturbf.visl )[ 0 ][ ug.lc ];
    turbcom.visl2 = ( * uturbf.visl )[ 0 ][ ug.rc ];
}

void UTurbSpectrum::PrepareFaceValue2Equ()
{
    gcom.xfn   = ( * ug.xfn   )[ ug.fId ];
    gcom.yfn   = ( * ug.yfn   )[ ug.fId ];
    gcom.zfn   = ( * ug.zfn   )[ ug.fId ];
    gcom.vfn   = ( * ug.vfn   )[ ug.fId ];
    gcom.farea = ( * ug.farea )[ ug.fId ];

    gcom.cvol1 = ( * ug.cvol )[ ug.lc ];
    gcom.cvol2 = ( * ug.cvol )[ ug.rc ];

    for ( int iEqu = 0; iEqu < nscom.nEqu; ++ iEqu )
    {
        turbsp.q1_ns[ iEqu ] = ( * uturbf.meanFlow.q )[ iEqu ][ ug.lc ];
        turbsp.q2_ns[ iEqu ] = ( * uturbf.meanFlow.q )[ iEqu ][ ug.rc ];
    }

    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turbsp.q1[ iEqu ] = ( * uturbf.q )[ iEqu ][ ug.lc ];
        turbsp.q2[ iEqu ] = ( * uturbf.q )[ iEqu ][ ug.rc ];
    }

    turbcom.visl1 = ( * uturbf.visl )[ 0 ][ ug.lc ];
    turbcom.visl2 = ( * uturbf.visl )[ 0 ][ ug.rc ];

    turbcom.vist1 = ( * uturbf.vist )[ 0 ][ ug.lc ];
    turbcom.vist2 = ( * uturbf.vist )[ 0 ][ ug.rc ];

    turbcom.bld1 = ( * uturbf.bld )[ 0 ][ ug.lc ];
    turbcom.bld2 = ( * uturbf.bld )[ 0 ][ ug.rc ];
}

EndNameSpace
