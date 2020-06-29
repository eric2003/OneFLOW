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

#include "NsVisFlux.h"
#include "NsCom.h"
#include "UCom.h"
#include "HXMath.h"
#include "BcRecord.h"
#include "Boundary.h"
#include "Grad.h"

BeginNameSpace( ONEFLOW )

NsVis vis;

VisGrad visQ;
VisGrad visT;

NsVis::NsVis()
{
    ;
}

NsVis::~NsVis()
{
    ;
}

void NsVis::Init()
{
    fvis.resize( nscom.nEqu );
}

NsVisFlux::NsVisFlux()
{
    ;
}

NsVisFlux::~NsVisFlux()
{
    ;
}

void NsVisFlux::AverGrad()
{
    visQ.AverGrad();
    visT.AverGrad();
}

void NsVisFlux::ZeroNormalGrad()
{
    visQ.ZeroNormalGrad();
    visT.ZeroNormalGrad();
}

void NsVisFlux::AverFaceValue()
{
    visQ.AverFaceValue();
    visT.AverFaceValue();
    this->AverOtherFaceValue();
}

void NsVisFlux::AverOtherFaceValue()
{
    nscom.visl   = half * ( nscom.visl1 + nscom.visl2 );
    nscom.vist   = half * ( nscom.vist1 + nscom.vist2 );
}

void NsVisFlux::AccurateFaceValue()
{
    this->AccurateOtherFaceValue();
    visQ.AccurateFaceValue();
    visT.AccurateFaceValue();
}

void NsVisFlux::AccurateOtherFaceValue()
{
    nscom.visl   = half * ( nscom.visl1 + nscom.visl2 );
    nscom.vist   = half * ( nscom.vist1 + nscom.vist2 );
}

void NsVisFlux::CorrectFaceGrad()
{
    visQ.CorrectFaceGrad();
    visT.CorrectFaceGrad();
}

void NsVisFlux::CalcNormalGrad()
{
    visQ.CalcNormalGrad();
    visT.CalcNormalGrad();
}

void NsVisFlux::CalcTestMethod()
{
    visQ.CalcTestMethod();
    visT.CalcTestMethod();
}

void NsVisFlux::CalcNew1Method()
{
    visQ.CalcNew1Method();
    visT.CalcNew1Method();
}

void NsVisFlux::CalcNew2Method()
{
    visQ.CalcNew2Method();
    visT.CalcNew2Method();
}

void NsVisFlux::ModifyFaceGrad()
{
    visQ.ModifyFaceGrad();
    visT.ModifyFaceGrad();
}

Real Sutherland::cdim = 110.4;
Real Sutherland::c;

Sutherland::Sutherland()
{
}

Sutherland::~Sutherland()
{

}

void Sutherland::CalcConst()
{
    Sutherland::c = Sutherland::cdim / nscom.tref_dim;
}

Real Sutherland::CalcViscosity( Real t )
{
    Real t3 = t * t * t;
    return sqrt( t3 ) * ( 1 + Sutherland::c ) / ( t + Sutherland::c );
}

EndNameSpace