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
#include "TurbSrcFlux.h"
#include "NsCtrl.h"
#include "DataBase.h"
#include "HXMath.h"
#include "TurbCom.h"
#include "Com.h"

BeginNameSpace( ONEFLOW )

TurbSrcFlux::TurbSrcFlux()
{
    ;
}

TurbSrcFlux::~TurbSrcFlux()
{
    ;
}

void TurbSrcFlux::SetSrcFluxPointer()
{
    if ( vis_model.visname.substr( 0, 6 ) == "2eq-kw" )
    {
        // MENTER'S K-OMEGA SST MODEL
        this->cmpBeta = & TurbSrcFlux::CmpFbetaDefault;
        if ( vis_model.visname.substr( 0, 13 ) == "2eq-kw-menter"  )
        {
            this->srcFlux = & TurbSrcFlux::CmpSrc2EquKwMenter;
            this->cmpBeta = & TurbSrcFlux::CmpFbetaDefault;
            this->cmpProd = & TurbSrcFlux::CmpProdwKwMenter;
        }
        else if ( vis_model.visname.substr( 0, 18 ) == "2eq-kw-wilcox-1998" )
        {
            this->srcFlux = & TurbSrcFlux::CmpSrc2EquKwWilcox1998;
            this->cmpBeta = & TurbSrcFlux::CmpFbetaOfKwWilcox1998;
            this->cmpProd = & TurbSrcFlux::CmpProdwKwWilcox1998;
        }
        else if ( vis_model.visname.substr( 0, 18 ) == "2eq-kw-wilcox-2006" )
        {
            this->srcFlux = & TurbSrcFlux::CmpSrc2EquKwWilcox2006;
            this->cmpBeta = & TurbSrcFlux::CmpFbetaOfKwWilcox2006;
            this->cmpProd = & TurbSrcFlux::CmpProdwKwWilcox2006;
        }
        else
        {
            this->srcFlux = & TurbSrcFlux::CmpSrc2EquKwDefault;
            this->cmpBeta = & TurbSrcFlux::CmpFbetaDefault;
            this->cmpProd = & TurbSrcFlux::CmpProdwKwDefault;
        }
    }
    else if ( vis_model.visname.substr( 0, 4 ) == "easm" )
    {
        if ( vis_model.visname.substr( 0, 12 ) == "easm-kw-2001" ||
             vis_model.visname.substr( 0, 12 ) == "easm-kw-2003" )
        {
            this->srcFlux = & TurbSrcFlux::CmpSrc2EquEasmKw2003;
            this->cmpBeta = & TurbSrcFlux::CmpFbetaOfEasmKw2003;
            this->cmpProd = & TurbSrcFlux::CmpProdwEasmKw2003;
        }
        else if ( vis_model.visname.substr( 0, 12 ) == "easm-kw-2005" )
        {
            this->srcFlux = & TurbSrcFlux::CmpSrc2EquEasmKw2005;
            this->cmpBeta = & TurbSrcFlux::CmpFbetaDefault;
            this->cmpProd = & TurbSrcFlux::CmpProdwEasmKw2005;
        }
    }
}

void TurbSrcFlux::CmpSrcSa()
{
    turbcom.CmpSrcSa();
}

void TurbSrcFlux::CmpFbetaCoef()
{
    ( this->* cmpBeta )();
}

void TurbSrcFlux::CmpProdW()
{
    ( this->* cmpProd )();
}


void TurbSrcFlux::CmpSrc2Equ()
{
    this->CmpVGrad();
    this->CmpTransition();
    this->CmpProdk();
    this->CmpFbetaCoef();
    this->CmpDissk();
    this->LimitProdk();
    this->CmpProdW();
    this->ModifyPd();
    this->CmpSrc();
}

void TurbSrcFlux::CmpSrc2EquKwMenter()
{
    this->CmpVGrad();
    this->CmpTransition();
    this->CmpProdk();
    this->CmpDissk();
    this->LimitProdk();
    this->CmpProdwKwMenter();
    this->ModifyPd();
    this->CmpSrc();
}

void TurbSrcFlux::CmpSrc2EquKwWilcox1998()
{
    this->CmpVGrad();
    this->CmpTransition();
    this->CmpProdk();
    this->CmpFbetaOfKwWilcox1998();
    this->CmpDissk();
    this->LimitProdk();
    this->CmpProdwKwWilcox1998();
    this->ModifyPd();
    this->CmpSrc();
}

void TurbSrcFlux::CmpSrc2EquKwWilcox2006()
{
    this->CmpVGrad();
    this->CmpTransition();
    this->CmpProdk();
    this->CmpFbetaOfKwWilcox2006();
    this->CmpDissk();
    this->LimitProdk();
    this->CmpProdwKwWilcox2006();
    this->ModifyPd();
    this->CmpSrc();
}

void TurbSrcFlux::CmpSrc2EquKwDefault()
{
    this->CmpVGrad();
    this->CmpTransition();
    this->CmpProdk();
    this->CmpDissk();
    this->LimitProdk();
    this->CmpProdwKwDefault();
    this->ModifyPd();
    this->CmpSrc();
}

void TurbSrcFlux::CmpSrc2EquEasmKw2003()
{
    this->CmpVGrad();
    this->CmpTransition();
    this->CmpProdk();
    this->CmpFbetaOfEasmKw2003();
    this->CmpDissk();
    this->LimitProdk();
    this->CmpProdwEasmKw2003();
    this->ModifyPd();
    this->CmpSrc();
}

void TurbSrcFlux::CmpSrc2EquEasmKw2005()
{
    this->CmpVGrad();
    this->CmpTransition();
    this->CmpProdk();
    this->CmpDissk();
    this->LimitProdk();
    this->CmpProdwEasmKw2005();
    this->ModifyPd();
    this->CmpSrc();
}


void TurbSrcFlux::CmpFbetaDefault()
{
}

void TurbSrcFlux::CmpFbetaOfKwWilcox1998()
{
    turbcom.CmpFbetaOfKwWilcox1998();
}

void TurbSrcFlux::CmpFbetaOfKwWilcox2006()
{
    turbcom.CmpFbetaOfKwWilcox2006();
}

void TurbSrcFlux::CmpFbetaOfEasmKw2003()
{
    turbcom.CmpFbetaOfEasmKw2003();
}

void TurbSrcFlux::CmpVGrad()
{
    turbcom.CmpVGrad();
}

void TurbSrcFlux::CmpTransition()
{
    turbcom.RGamaTransition();
}

void TurbSrcFlux::CmpProdk()
{
    turbcom.CmpProdk();
}

void TurbSrcFlux::CmpDissk()
{
    turbcom.CmpDissk();
}

void TurbSrcFlux::LimitProdk()
{
    turbcom.LimitProdk();
}

void TurbSrcFlux::CmpProdwKwMenter()
{
    turbcom.CmpProdwKwMenter();
}

void TurbSrcFlux::CmpProdwKwWilcox1998()
{
    turbcom.CmpProdwKwWilcox1998();
}

void TurbSrcFlux::CmpProdwKwWilcox2006()
{
    turbcom.CmpProdwKwWilcox2006();
}

void TurbSrcFlux::CmpProdwKwDefault()
{
    turbcom.CmpProdwKwDefault();
}

void TurbSrcFlux::CmpProdwEasmKw2003()
{
    turbcom.CmpProdwEasmKw2003();
}

void TurbSrcFlux::CmpProdwEasmKw2005()
{
    turbcom.CmpProdwEasmKw2005();
}

void TurbSrcFlux::ModifyPd()
{
    turbcom.ModifyPd();
}

void TurbSrcFlux::CmpSrc()
{
    turbcom.CmpSrc();
}

void TurbSrcFlux::CmpCellVist1Equ()
{
    Real olam = turbcom.rho / ( turbcom.visl + SMALL );
    turbcom.xsi  = turbcom.nuet * olam;
    turbcom.xsi3 = POWER3( turbcom.xsi );

    turbcom.fv1  = turbcom.xsi3 / ( turbcom.xsi3 + turbcom.cv13 );

    turbcom.vist = turbcom.fv1 * turbcom.rho * turbcom.nuet;

    turbcom.vist = MIN( turbcom.visl * turbcom.max_vis_ratio, turbcom.vist );
}

void TurbSrcFlux::CmpCellVist2Equ()
{
    Real dist2 = SQR( turbcom.dist );
    Real part1 = 2.0   * sqrt( turbcom.ke ) / ( turbcom.betas * turbcom.kw * turbcom.dist );
    Real part2 = 500.0 * turbcom.visl  / ( turbcom.rho  * turbcom.kw * dist2 * turbcom.reynolds );
    Real arg2  = MAX( part1 , part2 );
    Real f2    = tanh( arg2 * arg2 );

    Real vort = 0.0;

    if ( turbcom.iprod_sst == 1 )
    {
        Real s11 = turbcom.dudx;
        Real s22 = turbcom.dvdy;
        Real s33 = turbcom.dwdz;
        Real s12 = half * ( turbcom.dudy + turbcom.dvdx );
        Real s13 = half * ( turbcom.dudz + turbcom.dwdx );
        Real s23 = half * ( turbcom.dvdz + turbcom.dwdy );
        Real sij2         = two * ( SQR( s11, s22, s33 ) + two * SQR( s12, s13, s23 ) );
        Real magnitudeSij = sqrt( sij2 );
        vort = magnitudeSij;
    }
    else
    {
        vort = DIST(  turbcom.dwdy - turbcom.dvdz,  turbcom.dudz - turbcom.dwdx,  turbcom.dvdx - turbcom.dudy );
    }

    turbcom.vist = turbcom.rho * turbcom.ke / MAX( turbcom.kw, vort * f2 / turbcom.a1 ) * turbcom.reynolds;
    turbcom.vist = MIN( turbcom.visl * turbcom.max_vis_ratio, turbcom.vist );
}



EndNameSpace