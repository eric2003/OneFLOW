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

#include "UNsVisFlux.h"
#include "HeatFlux.h"
#include "Zone.h"
#include "ZoneState.h"
#include "UnsGrid.h"
#include "DataBase.h"
#include "NsCtrl.h"
#include "UCom.h"
#include "UNsCom.h"
#include "NsCom.h"
#include "VisGrad.h"
#include "UNsGrad.h"
#include "NsIdx.h"
#include "HXMath.h"
#include "Boundary.h"
#include "BcRecord.h"
#include "ULimiter.h"
#include "FieldImp.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )


UNsVisFlux::UNsVisFlux()
{
    ;
}

UNsVisFlux::~UNsVisFlux()
{
    ;
}

void UNsVisFlux::SetVisPointer()
{
    if ( nscom.ivischeme == VIS_STD )
    {
        this->visPointer = & UNsVisFlux::StdMethod;
    }
    else if ( nscom.ivischeme == VIS_TEST )
    {
        this->visPointer = & UNsVisFlux::TestMethod;
    }
    else if ( nscom.ivischeme == VIS_NEW1 )
    {
        this->visPointer = & UNsVisFlux::New1Method;
    }
    else if ( nscom.ivischeme == VIS_NEW2 )
    {
        this->visPointer = & UNsVisFlux::New2Method;
    }
    else if ( nscom.ivischeme == VIS_AVER )
    {
        this->visPointer = & UNsVisFlux::AverMethod;
    }
    else
    {
        this->visPointer = & UNsVisFlux::TestMethod;
    }
}

void UNsVisFlux::CalcFlux()
{
    if ( vis_model.vismodel == 0 ) return;
    ug.Init();
    unsf.Init();
    visQ.Init( nscom.nEqu );
    visT.Init( nscom.nTModel );
    vis.Init();
    heat_flux.Init();

    Alloc();

    this->SetVisPointer();

    this->PrepareField();
    this->CalcVisFlux();
    this->AddVisFlux();

    DeAlloc();
}

void UNsVisFlux::Alloc()
{
    visflux = new MRField( nscom.nEqu, ug.nFace );
}

void UNsVisFlux::DeAlloc()
{
    delete visflux;
}

void UNsVisFlux::PrepareField()
{
    ut_grad.Init();
    ut_grad.CalcGrad();
    //ut_grad.CalcGradDebug();
}

void UNsVisFlux::CalcVisFlux()
{
    for ( int fId = 0; fId < ug.nFace; ++ fId )
    {
        ug.fId = fId;

        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        if ( fId == 147489 )
        {
            int kkk = 1;
        }

        if ( ug.lc == 11 || ug.rc == 11 )
        {
        }

        this->PrepareFaceValue();

        this->CalcFaceVisFlux();

        this->UpdateFaceVisFlux();
    }
}

void UNsVisFlux::CalcFaceVisFlux()
{
    this->CalcHeatFlux();

    this->CalcStress();

    this->CalcNsVisFlux();
}

void UNsVisFlux::CalcHeatFlux()
{
    this->ZeroHeatFlux();

    this->AddChemHeatFlux();

    this->AddHeatFlux();

    SaveHeatFlux();
}

void UNsVisFlux::SaveHeatFlux()
{
    if ( ug.fId >= ug.nBFace ) return;
    if ( ug.bcRecord->bcType[ ug.fId ] != BC::SOLID_SURFACE ) return;
    SurfaceValue * heat_sur = heat_flux.heatflux[ ZoneState::zid ];
    Real non_dim_heatflux = - nscom.oreynolds * vis.qNormal;
    heat_sur->var->push_back( non_dim_heatflux );
}

void UNsVisFlux::CalcStress()
{
    Real divv2p3 = two3rd * ( vis.dudx + vis.dvdy + vis.dwdz );

    vis.txx = nscom.vis * ( two * vis.dudx - divv2p3 );
    vis.tyy = nscom.vis * ( two * vis.dvdy - divv2p3 );
    vis.tzz = nscom.vis * ( two * vis.dwdz - divv2p3 );
    vis.txy = nscom.vis * ( vis.dudy + vis.dvdx );
    vis.txz = nscom.vis * ( vis.dudz + vis.dwdx );
    vis.tyz = nscom.vis * ( vis.dvdz + vis.dwdy );

    this->CalcAniStress();
}

void UNsVisFlux::CalcAniStress()
{
    if ( ctrl.nrokplus <= 0 ) return;
    Real two3rdRhok = two3rd * vis.rhok;
    vis.txx += vis.b11 - two3rdRhok;
    vis.tyy += vis.b22 - two3rdRhok;
    vis.tzz += vis.b33 - two3rdRhok;
    vis.txy += vis.b12;
    vis.txz += vis.b13;
    vis.tyz += vis.b23;
}

void UNsVisFlux::CalcNsVisFlux()
{
    vis.fvis[ IDX::IR  ] = 0.0;
    vis.fvis[ IDX::IRU ] = gcom.xfn * vis.txx + gcom.yfn * vis.txy + gcom.zfn * vis.txz;
    vis.fvis[ IDX::IRV ] = gcom.xfn * vis.txy + gcom.yfn * vis.tyy + gcom.zfn * vis.tyz;
    vis.fvis[ IDX::IRW ] = gcom.xfn * vis.txz + gcom.yfn * vis.tyz + gcom.zfn * vis.tzz;
    vis.fvis[ IDX::IRE ] = vis.um * vis.fvis[ IDX::IRU ] + 
                           vis.vm * vis.fvis[ IDX::IRV ] + 
                           vis.wm * vis.fvis[ IDX::IRW ] + vis.qNormal;
}

void UNsVisFlux::ZeroHeatFlux()
{
    vis.qNormal = 0.0;
    vis.qx      = 0.0;
    vis.qy      = 0.0;
    vis.qz      = 0.0;
}

void UNsVisFlux::AddChemHeatFlux()
{
    if ( nscom.chemModel == 1 )
    {
    }
}

void UNsVisFlux::AddHeatFlux()
{
    vis.qNormal = 0.0;
    vis.qx      = 0.0;
    vis.qy      = 0.0;
    vis.qz      = 0.0;

    nscom.kcp = ( nscom.visl * nscom.oprl + nscom.vist * nscom.oprt ) * nscom.const_cp;
    vis.qNormal += gcom.xfn * vis.qx + gcom.yfn * vis.qy + gcom.zfn * vis.qz;
    vis.qNormal += nscom.kcp * vis.dtdn;
}

void UNsVisFlux::AddVisFlux()
{
    UnsGrid * grid = Zone::GetUnsGrid();
    MRField * res = GetFieldPointer< MRField >( grid, "res" );

    ONEFLOW::AddF2CField( res, visflux );
}

void UNsVisFlux::PrepareFaceValue()
{
    gcom.xfn   = ( * ug.xfn   )[ ug.fId ];
    gcom.yfn   = ( * ug.yfn   )[ ug.fId ];
    gcom.zfn   = ( * ug.zfn   )[ ug.fId ];
    gcom.vfn   = ( * ug.vfn   )[ ug.fId ];
    gcom.farea = ( * ug.farea )[ ug.fId ];

    gcom.CalcTangent();

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        visQ.dqdx1[ iEqu ] = ( * unsf.dqdx )[ iEqu ][ ug.lc ];
        visQ.dqdy1[ iEqu ] = ( * unsf.dqdy )[ iEqu ][ ug.lc ];
        visQ.dqdz1[ iEqu ] = ( * unsf.dqdz )[ iEqu ][ ug.lc ];

        visQ.dqdx2[ iEqu ] = ( * unsf.dqdx )[ iEqu ][ ug.rc ];
        visQ.dqdy2[ iEqu ] = ( * unsf.dqdy )[ iEqu ][ ug.rc ];
        visQ.dqdz2[ iEqu ] = ( * unsf.dqdz )[ iEqu ][ ug.rc ];
    }

    for ( int iEqu = 0; iEqu < nscom.nTModel; ++ iEqu )
    {
        visT.dqdx1[ iEqu ] = ( * unsf.dtdx )[ iEqu ][ ug.lc ];
        visT.dqdy1[ iEqu ] = ( * unsf.dtdy )[ iEqu ][ ug.lc ];
        visT.dqdz1[ iEqu ] = ( * unsf.dtdz )[ iEqu ][ ug.lc ];

        visT.dqdx2[ iEqu ] = ( * unsf.dtdx )[ iEqu ][ ug.rc ];
        visT.dqdy2[ iEqu ] = ( * unsf.dtdy )[ iEqu ][ ug.rc ];
        visT.dqdz2[ iEqu ] = ( * unsf.dtdz )[ iEqu ][ ug.rc ];
    }

    nscom.visl1 = ( * unsf.visl )[ 0 ][ ug.lc ];
    nscom.visl2 = ( * unsf.visl )[ 0 ][ ug.rc ];

    nscom.vist1 = ( * unsf.vist )[ 0 ][ ug.lc ];
    nscom.vist2 = ( * unsf.vist )[ 0 ][ ug.rc ];

    nscom.visl = half * ( nscom.visl1 + nscom.visl2 );
    nscom.vist = half * ( nscom.vist1 + nscom.vist2 );
    nscom.vis  = nscom.visl + nscom.vist;

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        visQ.q1[ iEqu ] = ( * unsf.q )[ iEqu ][ ug.lc ];
        visQ.q2[ iEqu ] = ( * unsf.q )[ iEqu ][ ug.rc ];
    }

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        visQ.q11[ iEqu ] = visQ.q1[ iEqu ];
        visQ.q22[ iEqu ] = visQ.q2[ iEqu ];
    }

    for ( int iEqu = 0; iEqu < nscom.nTModel; ++ iEqu )
    {
        visT.q1[ iEqu ] = ( * unsf.tempr )[ iEqu ][ ug.lc ];
        visT.q2[ iEqu ] = ( * unsf.tempr )[ iEqu ][ ug.rc ];
    }

    for ( int iEqu = 0; iEqu < nscom.nTModel; ++ iEqu )
    {
        visT.q11[ iEqu ] = visT.q1[ iEqu ];
        visT.q22[ iEqu ] = visT.q2[ iEqu ];
    }

    this->AverGrad();
    this->CalcFaceWeight();

    ( this->* visPointer )();

    this->SaveFacePara();
}

void UNsVisFlux::SaveFacePara()
{
    vis.dudx  = visQ.dqdx[ IDX::IU ];
    vis.dudy  = visQ.dqdy[ IDX::IU ];
    vis.dudz  = visQ.dqdz[ IDX::IU ];

    vis.dvdx  = visQ.dqdx[ IDX::IV ];
    vis.dvdy  = visQ.dqdy[ IDX::IV ];
    vis.dvdz  = visQ.dqdz[ IDX::IV ];

    vis.dwdx  = visQ.dqdx[ IDX::IW ];
    vis.dwdy  = visQ.dqdy[ IDX::IW ];
    vis.dwdz  = visQ.dqdz[ IDX::IW ];

    vis.um  = visQ.q[ IDX::IU ];
    vis.vm  = visQ.q[ IDX::IV ];
    vis.wm  = visQ.q[ IDX::IW ];

    vis.dtdn = visT.dqdn[ IDX::ITT ];
    vis.tmid = visT.q[ IDX::ITT ];
}

void UNsVisFlux::CalcFaceWeight()
{
    vgg.CalcFaceWeight();
}

void UNsVisFlux::AverMethod()
{
    this->ZeroNormalGrad();

    this->AverFaceValue();

    this->AverGrad();

    this->CalcNormalGrad();
}

void UNsVisFlux::StdMethod()
{
    this->CalcGradCoef();

    this->ZeroNormalGrad();

    this->AverFaceValue();

    this->AverGrad();

    this->CorrectFaceGrad();

    this->CalcNormalGrad();
}

void UNsVisFlux::TestMethod()
{
    this->ZeroNormalGrad();

    this->AverFaceValue();

    this->PrepareCellGeom();

    this->CalcTestMethod();

    this->ModifyFaceGrad();
}

void UNsVisFlux::New1Method()
{
    this->ZeroNormalGrad();

    this->AccurateFaceValue();

    this->PrepareCellGeom();

    this->CalcNew1Method();

    this->ModifyFaceGrad();
}

void UNsVisFlux::New2Method()
{
    this->ZeroNormalGrad();

    this->AccurateFaceValue();

    this->PrepareCellGeom();

    this->CalcNew2Method();

    this->ModifyFaceGrad();
}

void UNsVisFlux::CalcGradCoef()
{
    vgg.CalcGradCoef();
}


void UNsVisFlux::PrepareCellGeom()
{
    vgg.PrepareCellGeom();
}

void UNsVisFlux::UpdateFaceVisFlux()
{
    Real coeff = - nscom.oreynolds * gcom.farea;
    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        ( * visflux )[ iEqu ][ ug.fId ] = coeff * vis.fvis[ iEqu ];
    }
}

void CalcLaminarViscosity( int flag )
{
    ug.Init();
    unsf.Init();
    ug.SetStEd( flag );

    Real minLimit = 0.0;

    for ( int cId = ug.ist; cId < ug.ied; ++ cId )
    {
        Real temperature = ( * unsf.tempr )[ IDX::ITT ][ cId ];
        Real visl = Sutherland::CalcViscosity( temperature );
        ( * unsf.visl )[ 0 ][ cId ] = MAX( minLimit, visl );
    }
}

EndNameSpace