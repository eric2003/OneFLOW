/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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
#include "Boundary.h"
#include "BcRecord.h"
#include "UNsLimiter.h"
#include "FieldImp.h"
#include "Iteration.h"
#include "TurbCom.h"
#include "UTurbCom.h"
#include <iostream>
#include <iomanip>
using namespace std;

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
    if ( ! file ) exit( 0 );

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
