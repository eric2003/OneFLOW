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
#include "FieldWrap.h"
#include "SolverDef.h"
#include "UCom.h"
#include "UTurbCom.h"
#include "UNsCom.h"
#include "HipProductionPolicy.h"
#ifdef ONEFLOW_ENABLE_HIP
#include "HipTurbulenceBackend.h"
#include <vector>
#endif

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
    MRField * res = FieldHome::GetResidualField( grid, TURB_SOLVER );

    ONEFLOW::AddF2CField( res, invflux.get() );
}

void UTurbInvFlux::CalcFlux()
{
    TurbInv & inv = turbInv;
    inv.Init();
    ug.Init();
    // UNsField::Init pins NS_SOLVER (safe while SolverState is TURB).
    unsf.Init();
    uturbf.Init();

    invflux = std::make_unique<MRField>( limiter->GetNEquations(), ug.nFaces );

    this->CalcInvFace();
    this->CalcInvFlux();
    this->AddInvFlux();

    invflux.reset();
}

void UTurbInvFlux::CalcInvFlux()
{
    // Device counterpart of PrepareFaceValue -> RoeFlux -> UpdateFaceInvFlux,
    // one thread per face.  invflux is written only if the whole view resolves,
    // so a declining device path leaves it for the host loop below.
#ifdef ONEFLOW_ENABLE_HIP
    if ( HipProductionFeatureEnabled( "ONEFLOW_ENABLE_UTURB_HIP_INVFLUX" )
         && turbcom.nEqu == 1 && invflux != nullptr
         // The backend reports whether it can take these views at all; the face
         // entries currently cannot, because their pointers are host addresses
         // and there is no upload path yet.  Ask rather than assume.
         && HipTurbulenceBackend::Shared().SupportsFaceFlux()
         // The kernel implements only the VIS_TEST normal-gradient variant,
         // which the shipped cases set (iturb_visflux = 2).  Declining keeps any
         // other scheme on the working host path instead of letting the view
         // validator throw.
         && turbcom.iturb_visflux == 2 )
    {
        UnsGrid * grid = Zone::GetUnsGrid();
        MRField * qf1 = limiter->GetLeftField();
        MRField * qf2 = limiter->GetRightField();
        MRField * ns_qf1 = nslimiter->GetLeftField();
        MRField * ns_qf2 = nslimiter->GetRightField();
        if ( qf1 != nullptr && qf2 != nullptr
             && ns_qf1 != nullptr && ns_qf2 != nullptr )
        {
            TurbSaFaceView view;
            view.nFaces = ug.nFaces;
            view.nBoundaryFaces = ug.nBFaces;
            // The per-cell planes below are indexed through leftCell/rightCell,
            // which reach into the ghost range, so the backend needs the full
            // extent to copy them.
            view.nCells = ug.nCells;
            view.nGhostCells = ug.nTCell - ug.nCells;
            view.xNormal = ug.xfn->data();
            view.yNormal = ug.yfn->data();
            view.zNormal = ug.zfn->data();
            view.meshVelocityNormal = ug.vfn->data();
            view.faceArea = ug.farea->data();
            view.leftCell = ug.lcf->data();
            view.rightCell = ug.rcf->data();
            view.xFace = ug.xfc->data();
            view.yFace = ug.yfc->data();
            view.zFace = ug.zfc->data();
            view.xCell = ug.xcc->data();
            view.yCell = ug.ycc->data();
            view.zCell = ug.zcc->data();
            // The turbulence limiter carries exactly one plane (its nEqu is
            // turbcom.nEqu == 1), so ISA is the only index that is in range.
            view.faceWorkingVariableLeft = ( * qf1 )[ ISA ].data();
            view.faceWorkingVariableRight = ( * qf2 )[ ISA ].data();
            // Density and velocity come from the NS limiter, matching
            // PrepareFaceValue:224-232.  Reading IR/IU/IV/IW out of qf1 would
            // run four planes past the end of a one-equation field.
            view.faceDensityLeft = ( * ns_qf1 )[ IDX::IR ].data();
            view.faceDensityRight = ( * ns_qf2 )[ IDX::IR ].data();
            view.faceVelocityLeft[ 0 ] = ( * ns_qf1 )[ IDX::IU ].data();
            view.faceVelocityLeft[ 1 ] = ( * ns_qf1 )[ IDX::IV ].data();
            view.faceVelocityLeft[ 2 ] = ( * ns_qf1 )[ IDX::IW ].data();
            view.faceVelocityRight[ 0 ] = ( * ns_qf2 )[ IDX::IU ].data();
            view.faceVelocityRight[ 1 ] = ( * ns_qf2 )[ IDX::IV ].data();
            view.faceVelocityRight[ 2 ] = ( * ns_qf2 )[ IDX::IW ].data();
            view.cellWorkingVariable = ( * uturbf.q )[ ISA ].data();
            view.cellLaminarViscosity = ( * uturbf.visl )[ 0 ].data();
            view.cellDensity = ( * uturbf.meanFlow.q )[ IDX::IR ].data();
            view.workingVariableGradientX = ( * uturbf.dqdx )[ ISA ].data();
            view.workingVariableGradientY = ( * uturbf.dqdy )[ ISA ].data();
            view.workingVariableGradientZ = ( * uturbf.dqdz )[ ISA ].data();
            view.densityWeight = 1.0;
            view.invSigma = turbcom.osigma;
            view.invReynolds = turbcom.oreynolds;
            view.viscousFluxVariant = turbcom.iturb_visflux;
            view.diffusiveFlux = nullptr;
            view.cacheKey = grid;

            invfluxScratch_.assign( static_cast< std::size_t >( ug.nFaces ), 0.0 );
            view.inviscidFlux = invfluxScratch_.data();
            HipTurbulenceBackend::Shared().CalcFaceFlux( view );
            // UpdateFaceInvFlux:172 scales by face area; the kernel leaves that
            // to the caller because the diffusive side applies its own factor.
            for ( int face = 0; face < ug.nFaces; ++ face )
            {
                ( * invflux )[ ISA ][ face ] =
                    invfluxScratch_[ face ] * ( * ug.farea )[ face ];
            }
            return;
        }
    }
#endif
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
