/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2020 He Xin and the OneFLOW contributors.
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
#include "UTurbGrad.h"
#include "TurbCom.h"
#include "UTurbCom.h"
#include "NsIdx.h"
#include "UTurbLimiter.h"
#include "Zone.h"
#include "Com.h"
#include "UnsGrid.h"
#include "DataBase.h"
#include "UCom.h"
#include "Ctrl.h"
#include "VisGrad.h"
#include "HXMath.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

VisGrad visTurb;

UTurbVisFlux::UTurbVisFlux()
{

}

UTurbVisFlux::~UTurbVisFlux()
{

}

void UTurbVisFlux::Alloc()
{
    visflux = new MRField( turbcom.nEqu, ug.nFace );
}

void UTurbVisFlux::DeAlloc()
{
    delete visflux;
}

void UTurbVisFlux::CalcVisFlux()
{
    ug.Init();
    uturbf.Init();
    visTurb.Init( turbcom.nEqu );
    this->SetVisPointer();
    Alloc();
    if ( turbcom.nEqu == 1 )
    {
        this->CalcVisFlux1Equ();
    }
    else if ( turbcom.nEqu >= 2 )
    {
        this->CalcVisFlux2Equ();
    }
    DeAlloc();
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
    for ( int fId = 0; fId < ug.nFace; ++ fId )
    {
        ug.fId = fId;

        if ( fId == 866 )
        {
            int kkk = 1;
        }

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
    for ( int fId = 0; fId < ug.nFace; ++ fId )
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
    MRField * res = GetFieldPointer< MRField >( grid, "turbres" );

    ONEFLOW::AddF2CField( res, visflux );
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

    turbcom.rho1 = ( * uturbf.q_ns )[ IDX::IR ][ ug.lc ];
    turbcom.rho2 = ( * uturbf.q_ns )[ IDX::IR ][ ug.rc ];

    this->AverGrad();
    this->CalcFaceWeight();

    ( this->* visPointer )();
}

void UTurbVisFlux::CalcFaceVisFlux1Equ()
{
    Real orl = 1.0 / ( turbcom.rho1 + SMALL );
    Real orr = 1.0 / ( turbcom.rho2 + SMALL );

    Real visl_nv = half * ( turbcom.visl1 * orl  + turbcom.visl2 * orr  );
    Real vist_nv = half * ( visTurb.q1[ IKE ] + visTurb.q2[ IKE ] );

    // flux
    turbcom.comVis[ IKE ] = ( visl_nv  + vist_nv ) * turbcom.osigma;

    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        if ( ABS( visTurb.dqdn[ iEqu ] ) > 1.0e-10 )
        {
            int kkk = 1;
        }
        turbcom.flux[ iEqu ] = - turbcom.oreynolds * turbcom.comVis[ iEqu ] * visTurb.dqdn[ iEqu ] * gcom.farea;
    }
}

void UTurbVisFlux::CalcFaceVisFlux2Equ()
{
    turbcom.CalcSigkw();

    Real visl = half * ( turbcom.visl1 + turbcom.visl2 );
    Real vist = half * ( turbcom.vist1 + turbcom.vist2 );

    turbcom.comVis[ IKE ] = visl + vist * turbcom.sigk;
    turbcom.comVis[ IKW ] = visl + vist * turbcom.sigw;

    //Re-gama transition model
    if ( turbcom.transition_model == ITReGama )
    {
        turbcom.comVis[ ITGama ] = visl + vist / turbcom.trans_df;
        turbcom.comVis[ ITRect ] = ( visl + vist ) * turbcom.trans_dct;
    }

    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turbcom.flux[ iEqu ] = - turbcom.oreynolds * turbcom.comVis[ iEqu ] * visTurb.dqdn[ iEqu ] * gcom.farea;
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