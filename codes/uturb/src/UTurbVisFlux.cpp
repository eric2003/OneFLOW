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

#include "UTurbVisFlux.h"
#include "TurbVisFluxPolicy.h"
#include "UTurbGrad.h"
#include "TurbCom.h"
#include "UTurbCom.h"
#include "NsIdx.h"
#include "UTurbLimiter.h"
#include "Zone.h"
#include "Com.h"
#include "UnsGrid.h"
#include "DataBase.h"
#include "FieldWrap.h"
#include "SolverDef.h"
#include "UCom.h"
#include "Ctrl.h"
#include "VisGrad.h"
#include "HXMath.h"
#include "HipProductionPolicy.h"
#ifdef ONEFLOW_ENABLE_HIP
#include "HipTurbulenceBackend.h"
#include <vector>
#endif
#include <iostream>


BeginNameSpace( ONEFLOW )

VisGrad visTurb;

UTurbVisFlux::UTurbVisFlux()
{

}

UTurbVisFlux::~UTurbVisFlux()
{

}

void UTurbVisFlux::CalcVisFlux()
{
    ug.Init();
    uturbf.Init();
    visTurb.Init( turbcom.nEqu );
    this->SetVisPointer();
    visflux = std::make_unique<MRField>( turbcom.nEqu, ug.nFaces );
    if ( turbcom.nEqu == 1 )
    {
        this->CalcVisFlux1Equ();
    }
    else if ( turbcom.nEqu >= 2 )
    {
        this->CalcVisFlux2Equ();
    }
    visflux.reset();
}

void UTurbVisFlux::SetVisPointer()
{
    if ( turbcom.iturb_visflux == VIS_STD )
    {
        this->visPointer = & UTurbVisFlux::StdMethod;
    }
    else if ( turbcom.iturb_visflux == VIS_TEST )
    {
        this->visPointer = & UTurbVisFlux::TestMethod;
    }
    else if ( turbcom.iturb_visflux == VIS_NEW1 )
    {
        this->visPointer = & UTurbVisFlux::New1Method;
    }
    else if ( turbcom.iturb_visflux == VIS_NEW2 )
    {
        this->visPointer = & UTurbVisFlux::New2Method;
    }
    else if ( turbcom.iturb_visflux == VIS_AVER )
    {
        this->visPointer = & UTurbVisFlux::AverMethod;
    }
    else
    {
        this->visPointer = & UTurbVisFlux::TestMethod;
    }
}

void UTurbVisFlux::AverGrad()
{
    visTurb.AverGrad();
}

void UTurbVisFlux::ZeroNormalGrad()
{
    visTurb.ZeroNormalGrad();
}

void UTurbVisFlux::AverFaceValue()
{
    visTurb.AverFaceValue();
    this->AverOtherFaceValue();
}

void UTurbVisFlux::AverOtherFaceValue()
{
}

void UTurbVisFlux::CalcNormalGrad()
{
    visTurb.CalcNormalGrad();
}

void UTurbVisFlux::CalcTestMethod()
{
    visTurb.CalcTestMethod();
}

void UTurbVisFlux::CalcNew1Method()
{
    visTurb.CalcNew1Method();
}

void UTurbVisFlux::CalcNew2Method()
{
    visTurb.CalcNew2Method();
}

void UTurbVisFlux::CorrectFaceGrad()
{
    visTurb.CorrectFaceGrad();
}

void UTurbVisFlux::ModifyFaceGrad()
{
    visTurb.ModifyFaceGrad();
}

void UTurbVisFlux::AccurateFaceValue()
{
    this->AccurateOtherFaceValue();
    visTurb.AccurateFaceValue();
}

void UTurbVisFlux::AccurateOtherFaceValue()
{

}

void UTurbVisFlux::AverMethod()
{
    this->ZeroNormalGrad();

    this->AverFaceValue();

    this->AverGrad();

    this->CalcNormalGrad();
}

void UTurbVisFlux::StdMethod()
{
    this->CalcGradCoef();

    this->ZeroNormalGrad();

    this->AverFaceValue();

    this->AverGrad();

    this->CorrectFaceGrad();

    this->CalcNormalGrad();
}

void UTurbVisFlux::TestMethod()
{
    this->ZeroNormalGrad();

    this->AverFaceValue();

    this->PrepareCellGeom();

    this->CalcTestMethod();

    this->ModifyFaceGrad();
}

void UTurbVisFlux::New1Method()
{
    this->ZeroNormalGrad();

    this->AccurateFaceValue();

    this->PrepareCellGeom();

    this->CalcNew1Method();

    this->ModifyFaceGrad();
}

void UTurbVisFlux::New2Method()
{
    this->ZeroNormalGrad();

    this->AccurateFaceValue();

    this->PrepareCellGeom();

    this->CalcNew2Method();

    this->ModifyFaceGrad();
}

void UTurbVisFlux::CalcGradCoef()
{
    vgg.CalcGradCoef();
}


void UTurbVisFlux::PrepareCellGeom()
{
    vgg.PrepareCellGeom();
}

void UTurbVisFlux::CalcVisFlux1Equ()
{
    // Device counterpart of PrepareFaceValue -> CalcFaceVisFlux1Equ ->
    // UpdateFaceVisFlux.  Only the diffusive side is requested here; the
    // inviscid face flux is UTurbInvFlux's.  visflux is left untouched when the
    // path declines, so the host loop below still owns it.
#ifdef ONEFLOW_ENABLE_HIP
    if ( HipProductionFeatureEnabled( "ONEFLOW_ENABLE_UTURB_HIP_VISFLUX" )
         && turbcom.nEqu == 1 && visflux != nullptr
         // The kernel implements only the VIS_TEST normal-gradient variant,
         // which is what the shipped cases set (iturb_visflux = 2).  Any other
         // scheme must run on the host rather than be silently substituted, so
         // decline here instead of letting the view validator throw: a throw
         // would abort a case that works fine on the host path.
         && turbcom.iturb_visflux == 2 )
    {
        UnsGrid * grid = Zone::GetUnsGrid();
        {
            TurbSaFaceView view;
            view.nFaces = ug.nFaces;
            view.nBoundaryFaces = ug.nBFaces;
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
            // The diffusive flux reads the cell-centred working variable at each
            // face's left/right cell (PrepareFaceValue:303-311 puts
            // (*uturbf.q)[IKE][lc] into visTurb.q1), not the reconstructed face
            // value.  For the one-equation model IKE == ISA.
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
            view.inviscidFlux = nullptr;
            view.cacheKey = grid;

            visfluxScratch_.assign( static_cast< std::size_t >( ug.nFaces ), 0.0 );
            view.diffusiveFlux = visfluxScratch_.data();
            HipTurbulenceBackend::Shared().CalcFaceFlux( view );
            for ( int face = 0; face < ug.nFaces; ++ face )
            {
                ( * visflux )[ ISA ][ face ] = visfluxScratch_[ face ];
            }
            this->AddVisFlux();
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

        this->CalcFaceVisFlux1Equ();

        this->UpdateFaceVisFlux();
    }
    this->AddVisFlux();
}

void UTurbVisFlux::CalcVisFlux2Equ()
{
    for ( int fId = 0; fId < ug.nFaces; ++ fId )
    {
        ug.fId = fId;

        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        this->PrepareFaceValue();

        this->CalcFaceVisFlux2Equ();

        this->UpdateFaceVisFlux();
    }
    this->AddVisFlux();
}

void UTurbVisFlux::AddVisFlux()
{
    UnsGrid * grid = Zone::GetUnsGrid();
    MRField * res = FieldHome::GetResidualField( grid, TURB_SOLVER );

    ONEFLOW::AddF2CField( res, visflux.get() );
}

void UTurbVisFlux::PrepareFaceValue()
{
    gcom.xfn   = ( * ug.xfn   )[ ug.fId ];
    gcom.yfn   = ( * ug.yfn   )[ ug.fId ];
    gcom.zfn   = ( * ug.zfn   )[ ug.fId ];
    gcom.vfn   = ( * ug.vfn   )[ ug.fId ];
    gcom.farea = ( * ug.farea )[ ug.fId ];

    gcom.CalcTangent();

    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        visTurb.dqdx1[ iEqu ] = ( * uturbf.dqdx )[ iEqu ][ ug.lc ];
        visTurb.dqdy1[ iEqu ] = ( * uturbf.dqdy )[ iEqu ][ ug.lc ];
        visTurb.dqdz1[ iEqu ] = ( * uturbf.dqdz )[ iEqu ][ ug.lc ];

        visTurb.dqdx2[ iEqu ] = ( * uturbf.dqdx )[ iEqu ][ ug.rc ];
        visTurb.dqdy2[ iEqu ] = ( * uturbf.dqdy )[ iEqu ][ ug.rc ];
        visTurb.dqdz2[ iEqu ] = ( * uturbf.dqdz )[ iEqu ][ ug.rc ];
    }

    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        visTurb.q1[ iEqu ] = ( * uturbf.q )[ iEqu ][ ug.lc ];
        visTurb.q2[ iEqu ] = ( * uturbf.q )[ iEqu ][ ug.rc ];
    }

    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        visTurb.q11[ iEqu ] = visTurb.q1[ iEqu ];
        visTurb.q22[ iEqu ] = visTurb.q2[ iEqu ];
    }

    turbcom.visl1 = ( * uturbf.visl )[ 0 ][ ug.lc ];
    turbcom.visl2 = ( * uturbf.visl )[ 0 ][ ug.rc ];

    turbcom.vist1 = ( * uturbf.vist )[ 0 ][ ug.lc ];
    turbcom.vist2 = ( * uturbf.vist )[ 0 ][ ug.rc ];

    if ( turbcom.sst_type )
    {
        turbcom.bld1 = ( * uturbf.bld )[ 0 ][ ug.lc ];
        turbcom.bld2 = ( * uturbf.bld )[ 0 ][ ug.rc ];
    }

    turbcom.rho1 = ( * uturbf.meanFlow.q )[ IDX::IR ][ ug.lc ];
    turbcom.rho2 = ( * uturbf.meanFlow.q )[ IDX::IR ][ ug.rc ];

    this->AverGrad();
    this->CalcFaceWeight();

    ( this->* visPointer )();
}

void UTurbVisFlux::CalcFaceVisFlux1Equ()
{
    const Real orl = 1.0 / ( turbcom.rho1 + SMALL );
    const Real orr = 1.0 / ( turbcom.rho2 + SMALL );

    const Real visl_nv = ArithmeticMean( turbcom.visl1 * orl, turbcom.visl2 * orr );
    const Real vist_nv = ArithmeticMean( visTurb.q1[ IKE ], visTurb.q2[ IKE ] );

    turbcom.comVis[ IKE ] = SaFaceEffectiveDiffusivity(
        visl_nv, vist_nv, turbcom.osigma );

    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turbcom.flux[ iEqu ] = TurbScalarDiffusiveFlux(
            turbcom.oreynolds, turbcom.comVis[ iEqu ],
            visTurb.dqdn[ iEqu ], gcom.farea );
    }
}

void UTurbVisFlux::CalcFaceVisFlux2Equ()
{
    turbcom.CalcSigkw();

    const Real visl = ArithmeticMean( turbcom.visl1, turbcom.visl2 );
    const Real vist = ArithmeticMean( turbcom.vist1, turbcom.vist2 );

    turbcom.comVis[ IKE ] = SstFaceEffectiveDiffusivity( visl, vist, turbcom.sigk );
    turbcom.comVis[ IKW ] = SstFaceEffectiveDiffusivity( visl, vist, turbcom.sigw );

    //Re-gama transition model
    if ( turbcom.transition_model == ITReGama )
    {
        turbcom.comVis[ ITGama ] = SstFaceEffectiveDiffusivity(
            visl, vist, one / turbcom.trans_df );
        turbcom.comVis[ ITRect ] = ( visl + vist ) * turbcom.trans_dct;
    }

    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turbcom.flux[ iEqu ] = TurbScalarDiffusiveFlux(
            turbcom.oreynolds, turbcom.comVis[ iEqu ],
            visTurb.dqdn[ iEqu ], gcom.farea );
    }
}

void UTurbVisFlux::UpdateFaceVisFlux()
{
    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        ( * visflux )[ iEqu ][ ug.fId ] = turbcom.flux[ iEqu ];
    }
}

void UTurbVisFlux::CalcFaceWeight()
{
    vgg.CalcFaceWeight();
}

EndNameSpace
