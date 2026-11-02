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
#include "TurbSrcFlux.h"
#include "TurbSstPolicy.h"
#include "TurbSaPolicy.h"
#include "TurbViscosityPolicy.h"
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
        this->calcBeta = & TurbSrcFlux::CalcFbetaDefault;
        if ( vis_model.visname.substr( 0, 13 ) == "2eq-kw-menter"  )
        {
            this->srcFlux = & TurbSrcFlux::CalcSrc2EquKwMenter;
            this->calcBeta = & TurbSrcFlux::CalcFbetaDefault;
            this->calcProd = & TurbSrcFlux::CalcProdwKwMenter;
        }
        else if ( vis_model.visname.substr( 0, 18 ) == "2eq-kw-wilcox-1998" )
        {
            this->srcFlux = & TurbSrcFlux::CalcSrc2EquKwWilcox1998;
            this->calcBeta = & TurbSrcFlux::CalcFbetaOfKwWilcox1998;
            this->calcProd = & TurbSrcFlux::CalcProdwKwWilcox1998;
        }
        else if ( vis_model.visname.substr( 0, 18 ) == "2eq-kw-wilcox-2006" )
        {
            this->srcFlux = & TurbSrcFlux::CalcSrc2EquKwWilcox2006;
            this->calcBeta = & TurbSrcFlux::CalcFbetaOfKwWilcox2006;
            this->calcProd = & TurbSrcFlux::CalcProdwKwWilcox2006;
        }
        else
        {
            this->srcFlux = & TurbSrcFlux::CalcSrc2EquKwDefault;
            this->calcBeta = & TurbSrcFlux::CalcFbetaDefault;
            this->calcProd = & TurbSrcFlux::CalcProdwKwDefault;
        }
    }
    else if ( vis_model.visname.substr( 0, 4 ) == "easm" )
    {
        if ( vis_model.visname.substr( 0, 12 ) == "easm-kw-2001" ||
             vis_model.visname.substr( 0, 12 ) == "easm-kw-2003" )
        {
            this->srcFlux = & TurbSrcFlux::CalcSrc2EquEasmKw2003;
            this->calcBeta = & TurbSrcFlux::CalcFbetaOfEasmKw2003;
            this->calcProd = & TurbSrcFlux::CalcProdwEasmKw2003;
        }
        else if ( vis_model.visname.substr( 0, 12 ) == "easm-kw-2005" )
        {
            this->srcFlux = & TurbSrcFlux::CalcSrc2EquEasmKw2005;
            this->calcBeta = & TurbSrcFlux::CalcFbetaDefault;
            this->calcProd = & TurbSrcFlux::CalcProdwEasmKw2005;
        }
    }
}

void TurbSrcFlux::CalcSrcSa()
{
    turbcom.CalcSrcSa();
}

void TurbSrcFlux::CalcFbetaCoef()
{
    ( this->* calcBeta )();
}

void TurbSrcFlux::CalcProdW()
{
    ( this->* calcProd )();
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
    SaEddyViscosityInput input;
    input.density = turbcom.rho;
    input.workingVariable = turbcom.nuet;
    input.laminarViscosity = turbcom.visl;
    input.cv1Cubed = turbcom.cv13;
    input.maxViscosityRatio = turbcom.max_vis_ratio;

    const SaEddyViscosityResult result = SaComputeEddyViscosity( input );
    turbcom.xsi = result.chi;
    turbcom.xsi3 = POWER3( result.chi );
    turbcom.fv1 = result.fv1;
    turbcom.vist = result.eddyViscosity;
}

void TurbSrcFlux::CalcCellVist2Equ()
{
    SstEddyViscosityInput input;
    input.kineticEnergy = turbcom.ke;
    input.omega = turbcom.kw;
    input.wallDistance = turbcom.dist;
    input.betaStar = turbcom.betas;
    input.laminarViscosity = turbcom.visl;
    input.density = turbcom.rho;
    input.reynolds = turbcom.reynolds;
    input.a1 = turbcom.a1;
    input.maxViscosityRatio = turbcom.max_vis_ratio;
    input.productionOption = turbcom.iprod_sst;
    input.dudx = turbcom.dudx;
    input.dudy = turbcom.dudy;
    input.dudz = turbcom.dudz;
    input.dvdx = turbcom.dvdx;
    input.dvdy = turbcom.dvdy;
    input.dvdz = turbcom.dvdz;
    input.dwdx = turbcom.dwdx;
    input.dwdy = turbcom.dwdy;
    input.dwdz = turbcom.dwdz;

    const SstEddyViscosityResult result = SstComputeEddyViscosity( input );
    turbcom.vist = result.eddyViscosity;
}



EndNameSpace
