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
        this->cmpBeta = & TurbSrcFlux::CalcFbetaDefault;
        if ( vis_model.visname.substr( 0, 13 ) == "2eq-kw-menter"  )
        {
            this->srcFlux = & TurbSrcFlux::CalcSrc2EquKwMenter;
            this->cmpBeta = & TurbSrcFlux::CalcFbetaDefault;
            this->cmpProd = & TurbSrcFlux::CalcProdwKwMenter;
        }
        else if ( vis_model.visname.substr( 0, 18 ) == "2eq-kw-wilcox-1998" )
        {
            this->srcFlux = & TurbSrcFlux::CalcSrc2EquKwWilcox1998;
            this->cmpBeta = & TurbSrcFlux::CalcFbetaOfKwWilcox1998;
            this->cmpProd = & TurbSrcFlux::CalcProdwKwWilcox1998;
        }
        else if ( vis_model.visname.substr( 0, 18 ) == "2eq-kw-wilcox-2006" )
        {
            this->srcFlux = & TurbSrcFlux::CalcSrc2EquKwWilcox2006;
            this->cmpBeta = & TurbSrcFlux::CalcFbetaOfKwWilcox2006;
            this->cmpProd = & TurbSrcFlux::CalcProdwKwWilcox2006;
        }
        else
        {
            this->srcFlux = & TurbSrcFlux::CalcSrc2EquKwDefault;
            this->cmpBeta = & TurbSrcFlux::CalcFbetaDefault;
            this->cmpProd = & TurbSrcFlux::CalcProdwKwDefault;
        }
    }
    else if ( vis_model.visname.substr( 0, 4 ) == "easm" )
    {
        if ( vis_model.visname.substr( 0, 12 ) == "easm-kw-2001" ||
             vis_model.visname.substr( 0, 12 ) == "easm-kw-2003" )
        {
            this->srcFlux = & TurbSrcFlux::CalcSrc2EquEasmKw2003;
            this->cmpBeta = & TurbSrcFlux::CalcFbetaOfEasmKw2003;
            this->cmpProd = & TurbSrcFlux::CalcProdwEasmKw2003;
        }
        else if ( vis_model.visname.substr( 0, 12 ) == "easm-kw-2005" )
        {
            this->srcFlux = & TurbSrcFlux::CalcSrc2EquEasmKw2005;
            this->cmpBeta = & TurbSrcFlux::CalcFbetaDefault;
            this->cmpProd = & TurbSrcFlux::CalcProdwEasmKw2005;
        }
    }
}

void TurbSrcFlux::CalcSrcSa()
{
    turbcom.CalcSrcSa();
}

void TurbSrcFlux::CalcFbetaCoef()
{
    ( this->* cmpBeta )();
}

void TurbSrcFlux::CalcProdW()
{
    ( this->* cmpProd )();
}


void TurbSrcFlux::CalcSrc2Equ()
{
    this->CalcVGrad();
    this->CalcTransition();
    this->CalcProdk();
    this->CalcFbetaCoef();
    this->CalcDissk();
    this->LimitProdk();
    this->CalcProdW();
    this->ModifyPd();
    this->CalcSrc();
}

void TurbSrcFlux::CalcSrc2EquKwMenter()
{
    this->CalcVGrad();
    this->CalcTransition();
    this->CalcProdk();
    this->CalcDissk();
    this->LimitProdk();
    this->CalcProdwKwMenter();
    this->ModifyPd();
    this->CalcSrc();
}

void TurbSrcFlux::CalcSrc2EquKwWilcox1998()
{
    this->CalcVGrad();
    this->CalcTransition();
    this->CalcProdk();
    this->CalcFbetaOfKwWilcox1998();
    this->CalcDissk();
    this->LimitProdk();
    this->CalcProdwKwWilcox1998();
    this->ModifyPd();
    this->CalcSrc();
}

void TurbSrcFlux::CalcSrc2EquKwWilcox2006()
{
    this->CalcVGrad();
    this->CalcTransition();
    this->CalcProdk();
    this->CalcFbetaOfKwWilcox2006();
    this->CalcDissk();
    this->LimitProdk();
    this->CalcProdwKwWilcox2006();
    this->ModifyPd();
    this->CalcSrc();
}

void TurbSrcFlux::CalcSrc2EquKwDefault()
{
    this->CalcVGrad();
    this->CalcTransition();
    this->CalcProdk();
    this->CalcDissk();
    this->LimitProdk();
    this->CalcProdwKwDefault();
    this->ModifyPd();
    this->CalcSrc();
}

void TurbSrcFlux::CalcSrc2EquEasmKw2003()
{
    this->CalcVGrad();
    this->CalcTransition();
    this->CalcProdk();
    this->CalcFbetaOfEasmKw2003();
    this->CalcDissk();
    this->LimitProdk();
    this->CalcProdwEasmKw2003();
    this->ModifyPd();
    this->CalcSrc();
}

void TurbSrcFlux::CalcSrc2EquEasmKw2005()
{
    this->CalcVGrad();
    this->CalcTransition();
    this->CalcProdk();
    this->CalcDissk();
    this->LimitProdk();
    this->CalcProdwEasmKw2005();
    this->ModifyPd();
    this->CalcSrc();
}


void TurbSrcFlux::CalcFbetaDefault()
{
}

void TurbSrcFlux::CalcFbetaOfKwWilcox1998()
{
    turbcom.CalcFbetaOfKwWilcox1998();
}

void TurbSrcFlux::CalcFbetaOfKwWilcox2006()
{
    turbcom.CalcFbetaOfKwWilcox2006();
}

void TurbSrcFlux::CalcFbetaOfEasmKw2003()
{
    turbcom.CalcFbetaOfEasmKw2003();
}

void TurbSrcFlux::CalcVGrad()
{
    turbcom.CalcVGrad();
}

void TurbSrcFlux::CalcTransition()
{
    turbcom.RGamaTransition();
}

void TurbSrcFlux::CalcProdk()
{
    turbcom.CalcProdk();
}

void TurbSrcFlux::CalcDissk()
{
    turbcom.CalcDissk();
}

void TurbSrcFlux::LimitProdk()
{
    turbcom.LimitProdk();
}

void TurbSrcFlux::CalcProdwKwMenter()
{
    turbcom.CalcProdwKwMenter();
}

void TurbSrcFlux::CalcProdwKwWilcox1998()
{
    turbcom.CalcProdwKwWilcox1998();
}

void TurbSrcFlux::CalcProdwKwWilcox2006()
{
    turbcom.CalcProdwKwWilcox2006();
}

void TurbSrcFlux::CalcProdwKwDefault()
{
    turbcom.CalcProdwKwDefault();
}

void TurbSrcFlux::CalcProdwEasmKw2003()
{
    turbcom.CalcProdwEasmKw2003();
}

void TurbSrcFlux::CalcProdwEasmKw2005()
{
    turbcom.CalcProdwEasmKw2005();
}

void TurbSrcFlux::ModifyPd()
{
    turbcom.ModifyPd();
}

void TurbSrcFlux::CalcSrc()
{
    turbcom.CalcSrc();
}

void TurbSrcFlux::CalcCellVist1Equ()
{
    Real olam = turbcom.rho / ( turbcom.visl + SMALL );
    turbcom.xsi  = turbcom.nuet * olam;
    turbcom.xsi3 = POWER3( turbcom.xsi );

    turbcom.fv1  = turbcom.xsi3 / ( turbcom.xsi3 + turbcom.cv13 );

    turbcom.vist = turbcom.fv1 * turbcom.rho * turbcom.nuet;

    turbcom.vist = MIN( turbcom.visl * turbcom.max_vis_ratio, turbcom.vist );
}

void TurbSrcFlux::CalcCellVist2Equ()
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