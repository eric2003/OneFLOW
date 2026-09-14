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
#include "UCom.h"
#include "UNsCom.h"
#include "NsCom.h"
#include "NsIdx.h"
#include "HXMath.h"
#include "Fatal.h"
#include "Boundary.h"
#include "BcRecord.h"
#include "UNsLimiter.h"
#include "FieldImp.h"
#include "Iteration.h"
#include "TurbCom.h"
#include "UTurbCom.h"
#include "AccelRuntime.h"
#include "EulerCpuAdapter.h"
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>


BeginNameSpace( ONEFLOW )

UNsInvFlux::UNsInvFlux()
{
    limiter = new NsLimiter();
    limf = limiter->limf;
}

UNsInvFlux::~UNsInvFlux()
{
    delete limiter;
}

void UNsInvFlux::CalcLimiter()
{
    limiter->CalcLimiter();
}

void UNsInvFlux::CalcInvFace()
{
    uns_grad.Init();
    uns_grad.CalcGrad();

    this->CalcLimiter();

    this->GetQlQrField();

    this->ReconstructFaceValueField();

    this->BoundaryQlQrFixField();
}

void UNsInvFlux::GetQlQrField()
{
    limf->GetQlQr();
}

void UNsInvFlux::ReconstructFaceValueField()
{
    limf->CalcFaceValue();
    //limf->CalcFaceValueWeighted();
    if ( Iteration::outerSteps == -31 )
    {
        Real mindiff = 1.0e-10;
        int idumpface = 1;
        int idumpcell = 0;

        HXDebug::DumpField( "limf.dqdx.debug", limf->dqdx );
        HXDebug::CompareFile( mindiff, idumpcell );
        HXDebug::DumpField( "limf.dqdy.debug", limf->dqdy );
        HXDebug::CompareFile( mindiff, idumpcell );
        HXDebug::DumpField( "limf.dqdz.debug", limf->dqdz );
        HXDebug::CompareFile( mindiff, idumpcell );

        HXDebug::DumpField( "limf.qf1_recon.debug", limf->qf1 );
        HXDebug::CompareFile( mindiff, idumpface );
        HXDebug::DumpField( "limf.qf2_recon.debug", limf->qf2 );
        HXDebug::CompareFile( mindiff, idumpface );
    }
}

void UNsInvFlux::BoundaryQlQrFixField()
{
    limf->BcQlQrFix();

    if ( Iteration::outerSteps == -31 )
    {
        Real mindiff = 1.0e-10;
        int idumpface = 1;
        int idumpcell = 0;

        HXDebug::DumpField( "limf.qf1.debug", limf->qf1 );
        HXDebug::CompareFile( mindiff, idumpface );
        HXDebug::DumpField( "limf.qf2.debug", limf->qf2 );
        HXDebug::CompareFile( mindiff, idumpface );
    }
}

void UNsInvFlux::CalcFlux()
{
    if ( nscom.icmpInv == 0 ) return;
    inv.Init();
    ug.Init();
    unsf.Init();
    Alloc();

    this->SetPointer( nscom.ischeme );

    //ReadTmp();
    this->CalcInvFace();
    this->CalcInvFlux();
    this->AddInvFlux();

    DeAlloc();
}

void UNsInvFlux::CalcInvFlux()
{
    if ( this->UseCpuBatchAdapter() )
    {
        this->CalcInvFluxCpuBatch();
        return;
    }

    for ( int fId = 0; fId < ug.nFaces; ++ fId )
    {
        ug.fId = fId;

        if ( fId == 24 )
        {
            int kkk = 1;
        }

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
        && nscom.nEqu == 5 && limf != nullptr && limf->nEqu == 5;
}

void UNsInvFlux::CalcInvFluxCpuBatch()
{
    const int nFaces = ug.nFaces;
    const int nEquations = limf->nEqu;
    std::vector< Real > primitiveLeft( nEquations * nFaces );
    std::vector< Real > primitiveRight( nEquations * nFaces );
    std::vector< Real > xNormal( nFaces );
    std::vector< Real > yNormal( nFaces );
    std::vector< Real > zNormal( nFaces );
    std::vector< Real > meshVelocityNormal( nFaces );
    std::vector< Real > faceArea( nFaces );
    std::vector< Real > faceFlux( nEquations * nFaces );

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
                ( * limf->qf1 )[ equation ][ face ];
            primitiveRight[ equation * nFaces + face ] =
                ( * limf->qf2 )[ equation ][ face ];
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

    EulerCpuAdapter adapter;
    adapter.CalcInvFlux( primitiveState, flux, 1 );
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

    for ( int iEqu = 0; iEqu < limf->nEqu; ++ iEqu )
    {
        inv.prim1[ iEqu ] = ( * limf->qf1 )[ iEqu ][ ug.fId ];
        inv.prim2[ iEqu ] = ( * limf->qf2 )[ iEqu ][ ug.fId ];
    }
}

void UNsInvFlux::UpdateFaceInvFlux()
{
    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        ( * invflux )[ iEqu ][ ug.fId ] = gcom.farea * inv.flux[ iEqu ];
    }
}

void UNsInvFlux::AddInvFlux()
{
    UnsGrid * grid = Zone::GetUnsGrid();
    MRField * res = GetFieldPointer< MRField >( grid, "res" );

    ONEFLOW::AddF2CField( res, invflux );
    if ( Iteration::outerSteps == -31 )
    {
        HXDebug::CheckNANField( res );
        Real mindiff = 1.0e-10;
        int idumpface = 1;
        int idumpcell = 0;
        MRField * q = GetFieldPointer< MRField >( grid, "q" );
        HXDebug::DumpField( "flow.debug", q );
        HXDebug::CompareFile( 1.0e-12, idumpcell );

        HXDebug::DumpField( "InvFaceFlux.debug", invflux );
        HXDebug::CompareFile( mindiff, idumpface );
        HXDebug::DumpResField( "InvResFlux.debug" );
        HXDebug::CompareFile( mindiff, idumpcell );
    }
}

void UNsInvFlux::Alloc()
{
    invflux = new MRField( nscom.nEqu, ug.nFaces );
}

void UNsInvFlux::DeAlloc()
{
    delete invflux;
}

void UNsInvFlux::ReadTmp()
{
    static int iii = 0;
    if ( iii ) return;
    iii = 1;
    std::fstream file;
    file.open( "nsflow.dat", std::ios_base::in | std::ios_base::binary );
    if ( ! file )
    {
        Fatal( "Failed to open file: nsflow.dat" );
    }

    unsf.Init();

    for ( int cId = 0; cId < ug.nTCell; ++ cId )
    {
        for ( int iEqu = 0; iEqu < 5; ++ iEqu )
        {
            file.read( reinterpret_cast< char * >( & ( * unsf.q )[ iEqu ][ cId ] ), sizeof( double ) );
        }
    }

    for ( int cId = 0; cId < ug.nTCell; ++ cId )
    {
        file.read( reinterpret_cast< char * >( & ( * unsf.visl )[ 0 ][ cId ] ), sizeof( double ) );
    }

    for ( int cId = 0; cId < ug.nTCell; ++ cId )
    {
        file.read( reinterpret_cast< char * >( & ( * unsf.vist )[ 0 ][ cId ] ), sizeof( double ) );
    }

    std::vector< Real > tmp1( ug.nTCell ), tmp2( ug.nTCell );

    for ( int cId = 0; cId < ug.nTCell; ++ cId )
    {
        tmp1[ cId ] = ( * unsf.timestep )[ 0 ][ cId ];
    }

    for ( int cId = 0; cId < ug.nTCell; ++ cId )
    {
        file.read( reinterpret_cast< char * >( & ( * unsf.timestep )[ 0 ][ cId ] ), sizeof( double ) );
    }

    for ( int cId = 0; cId < ug.nTCell; ++ cId )
    {
        tmp2[ cId ] = ( * unsf.timestep )[ 0 ][ cId ];
    }

    turbcom.Init();
    uturbf.Init();
    for ( int iCell = 0; iCell < ug.nTCell; ++ iCell )
    {
        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            file.read( reinterpret_cast< char * >( & ( * uturbf.q )[ iEqu ][ iCell ] ), sizeof( double ) );
        }
    }
    file.close();
    file.clear();
}



EndNameSpace
