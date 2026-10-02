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
#include "UTurbCom.h"
#include "UNsCom.h"

BeginNameSpace( ONEFLOW )

UTurbInvFlux::UTurbInvFlux()
{
    limiter = std::make_unique<TurbLimiter>();
    nslimiter = std::make_unique<NsLimiter>();
    nslimiter->limflag = turbcom.tns_ilim;
    limiter->limflag = turbcom.turb_ilim;
}

UTurbInvFlux::~UTurbInvFlux()
{
}

void UTurbInvFlux::CalcLimiter()
{
     limiter->CalcLimiter();
     nslimiter->CalcLimiter();
}

void UTurbInvFlux::CalcInvFace()
{
    this->CalcLimiter();
    this->GetQlQrField();

    this->ReconstructFaceValueField();

    this->BoundaryQlQrFixField();
}

void UTurbInvFlux::GetQlQrField()
{
    limiter->GetQlQr();
    nslimiter->GetQlQr();
}

void UTurbInvFlux::ReconstructFaceValueField()
{
    limiter->CalcFaceValue();
    nslimiter->CalcFaceValue();
}

void UTurbInvFlux::BoundaryQlQrFixField()
{
    limiter->BcQlQrFix();
    nslimiter->BcQlQrFix();
}

void UTurbInvFlux::AddInvFlux()
{
    UnsGrid * grid = Zone::GetUnsGrid();
    MRField * res = GetFieldPointer< MRField >( grid, "turbres" );

    ONEFLOW::AddF2CField( res, invflux );
}

void UTurbInvFlux::CalcFlux()
{
    TurbInv & inv = turbInv;
    inv.Init();
    ug.Init();
    unsf.Init();
    uturbf.Init();

    invflux = new MRField( limiter->GetNEquations(), ug.nFaces);

    this->CalcInvFace();
    this->CalcInvFlux();
    this->AddInvFlux();

    delete invflux;
}

void UTurbInvFlux::CalcInvFlux()
{
    for ( int fId = 0; fId < ug.nFaces; ++ fId )
    {
        ug.fId = fId;

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

    MRField * qf1 = limiter->GetLeftField();
    MRField * qf2 = limiter->GetRightField();

    int nEquations = limiter->GetNEquations();

    for ( int iEqu = 0; iEqu < nEquations; ++ iEqu )
    {
        inv.prim1[ iEqu ] = ( * qf1 )[ iEqu ][ ug.fId ];
        inv.prim2[ iEqu ] = ( * qf2 )[ iEqu ][ ug.fId ];
    }

    MRField * ns_qf1 = nslimiter->GetLeftField();
    MRField * ns_qf2 = nslimiter->GetRightField();

    inv.rl = ( * ns_qf1 )[ IDX::IR ][ ug.fId ];
    inv.ul = ( * ns_qf1 )[ IDX::IU ][ ug.fId ];
    inv.vl = ( * ns_qf1 )[ IDX::IV ][ ug.fId ];
    inv.wl = ( * ns_qf1 )[ IDX::IW ][ ug.fId ];

    inv.rr = ( * ns_qf2 )[ IDX::IR ][ ug.fId ];
    inv.ur = ( * ns_qf2 )[ IDX::IU ][ ug.fId ];
    inv.vr = ( * ns_qf2 )[ IDX::IV ][ ug.fId ];
    inv.wr = ( * ns_qf2 )[ IDX::IW ][ ug.fId ];
}

void UTurbInvFlux::UpdateFaceInvFlux()
{
    TurbInv & inv = turbInv;

    int nEquations = limiter->GetNEquations();

    for ( int iEqu = 0; iEqu < nEquations; ++ iEqu )
    {
        ( * invflux )[ iEqu ][ ug.fId ] = gcom.farea * inv.flux[ iEqu ];
    }
}

EndNameSpace
