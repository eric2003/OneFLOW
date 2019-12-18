/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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

#include "UTurbInvFlux.h"
#include "UTurbGrad.h"
#include "TurbCom.h"
#include "UNsGrad.h"
#include "UTurbLimiter.h"
#include "UNsLimiter.h"
#include "NsIdx.h"
#include "Zone.h"
#include "Com.h"
#include "UnsGrid.h"
#include "DataBase.h"
#include "UCom.h"

BeginNameSpace( ONEFLOW )

UTurbInvFlux::UTurbInvFlux()
{
    limiter = new TurbLimiter();
    nslimiter = new NsLimiter();
    nslimiter->limflag = turbcom.tns_ilim;
    limf = limiter->limf;
    limiter->limflag = turbcom.turb_ilim;
}

UTurbInvFlux::~UTurbInvFlux()
{
    delete limiter;
    delete nslimiter;
}

void UTurbInvFlux::CmpLimiter()
{
     limiter->CmpLimiter();
     nslimiter->CmpLimiter();
}

void UTurbInvFlux::CmpInvFace()
{
    this->CmpLimiter();
    this->GetQlQrField();

    this->ReconstructFaceValueField();

    this->BoundaryQlQrFixField();
}

void UTurbInvFlux::GetQlQrField()
{
    limf->GetQlQr();
    nslimiter->limf->GetQlQr();
}

void UTurbInvFlux::ReconstructFaceValueField()
{
    limf->CmpFaceValue();
    nslimiter->limf->CmpFaceValue();
}

void UTurbInvFlux::BoundaryQlQrFixField()
{
    limf->BcQlQrFix();
    nslimiter->limf->BcQlQrFix();
}

void UTurbInvFlux::AddInvFlux()
{
    UnsGrid * grid = Zone::GetUnsGrid();
    MRField * res = GetFieldPointer< MRField >( grid, "turbres" );

    ONEFLOW::AddF2CField( res, invflux );

    vector< vector< Real > > tmp( turbcom.nEqu );
    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        tmp[ iEqu ].resize( ug.nCell );
        for ( int cId = 0; cId < ug.nCell; ++ cId )
        {
            tmp[ iEqu ][ cId ] = ( * invflux  )[ iEqu ][ cId ];
        }
    }
    int kkk = 1;
    {
    vector< vector< Real > > tmp( turbcom.nEqu );
    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        tmp[ iEqu ].resize( ug.nCell );
        for ( int cId = 0; cId < ug.nCell; ++ cId )
        {
            tmp[ iEqu ][ cId ] = ( * res   )[ iEqu ][ cId ];
        }
    }
    int kkk = 1;
    }
}

void UTurbInvFlux::Alloc()
{
    invflux = new MRField( limf->nEqu, ug.nFace );
}

void UTurbInvFlux::DeAlloc()
{
    delete invflux;
}

void UTurbInvFlux::CmpFlux()
{
    TurbInv & inv = turbInv;
    inv.Init();
    ug.Init();
    Alloc();

    this->CmpInvFace();
    this->CmpInvFlux();
    this->AddInvFlux();

    DeAlloc();
}

void UTurbInvFlux::CmpInvFlux()
{
    for ( int fId = 0; fId < ug.nFace; ++ fId )
    {
        ug.fId = fId;

        if ( fId == 384 )
        {
            int kkk = 1;
        }

        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        this->PrepareFaceValue();
        this->RoeFlux();
        this->UpdateFaceInvFlux();
    }
}

void UTurbInvFlux::PrepareFaceValue()
{
    TurbInv & inv = turbInv;

    gcom.xfn   = ( * ug.xfn   )[ ug.fId ];
    gcom.yfn   = ( * ug.yfn   )[ ug.fId ];
    gcom.zfn   = ( * ug.zfn   )[ ug.fId ];
    gcom.vfn   = ( * ug.vfn   )[ ug.fId ];
    gcom.farea = ( * ug.farea )[ ug.fId ];

    for ( int iEqu = 0; iEqu < limf->nEqu; ++ iEqu )
    {
        inv.prim1[ iEqu ] = ( * limf->qf1 )[ iEqu ][ ug.fId ];
        inv.prim2[ iEqu ] = ( * limf->qf2 )[ iEqu ][ ug.fId ];
    }

    inv.rl = ( * nslimiter->limf->qf1 )[ IDX::IR ][ ug.fId ];
    inv.ul = ( * nslimiter->limf->qf1 )[ IDX::IU ][ ug.fId ];
    inv.vl = ( * nslimiter->limf->qf1 )[ IDX::IV ][ ug.fId ];
    inv.wl = ( * nslimiter->limf->qf1 )[ IDX::IW ][ ug.fId ];

    inv.rr = ( * nslimiter->limf->qf2 )[ IDX::IR ][ ug.fId ];
    inv.ur = ( * nslimiter->limf->qf2 )[ IDX::IU ][ ug.fId ];
    inv.vr = ( * nslimiter->limf->qf2 )[ IDX::IV ][ ug.fId ];
    inv.wr = ( * nslimiter->limf->qf2 )[ IDX::IW ][ ug.fId ];
}

void UTurbInvFlux::UpdateFaceInvFlux()
{
    TurbInv & inv = turbInv;

    for ( int iEqu = 0; iEqu < limf->nEqu; ++ iEqu )
    {
        ( * invflux )[ iEqu ][ ug.fId ] = gcom.farea * inv.flux[ iEqu ];
    }
}

EndNameSpace