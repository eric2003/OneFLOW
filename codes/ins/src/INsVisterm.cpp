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

#include "INsVisterm.h"
#include "INsCom.h"
#include "UCom.h"
#include "HXMath.h"
#include "BcRecord.h"
#include "Boundary.h"
#include "Grad.h"

BeginNameSpace( ONEFLOW )

INsVis Ivis;

INsVis::INsVis()
{
    ;
}

INsVis::~INsVis()
{
    ;
}

void INsVis::Init()
{
    fvis.resize( nscom.nEqu );
}

INsVisterm::INsVisterm()
{
    ;
}

INsVisterm::~INsVisterm()
{
    ;
}

void INsVisterm::AverGrad()
{
    visQ.AverGrad();
    visT.AverGrad();
}

void INsVisterm::ZeroNormalGrad()
{
    visQ.ZeroNormalGrad();
    visT.ZeroNormalGrad();
}

void INsVisterm::AverFaceValue()
{
    visQ.AverFaceValue();
    visT.AverFaceValue();
    this->AverOtherFaceValue();
}

void INsVisterm::AverOtherFaceValue()
{
    nscom.visl   = half * ( nscom.visl1 + nscom.visl2 );
    nscom.vist   = half * ( nscom.vist1 + nscom.vist2 );
}

void INsVisterm::AccurateFaceValue()
{
    this->AccurateOtherFaceValue();
    visQ.AccurateFaceValue();
    visT.AccurateFaceValue();
}

void INsVisterm::AccurateOtherFaceValue()
{
    nscom.visl   = half * ( nscom.visl1 + nscom.visl2 );
    nscom.vist   = half * ( nscom.vist1 + nscom.vist2 );
}

void INsVisterm::CorrectFaceGrad()
{
    visQ.CorrectFaceGrad();
    visT.CorrectFaceGrad();
}

void INsVisterm::CalcNormalGrad()
{
    visQ.CalcNormalGrad();
    visT.CalcNormalGrad();
}

void INsVisterm::CalcTestMethod()
{
    visQ.CalcTestMethod();
    visT.CalcTestMethod();
}

void INsVisterm::CalcNew1Method()
{
    visQ.CalcNew1Method();
    visT.CalcNew1Method();
}

void INsVisterm::CalcNew2Method()
{
    visQ.CalcNew2Method();
    visT.CalcNew2Method();
}

void INsVisterm::ModifyFaceGrad()
{
    visQ.ModifyFaceGrad();
    visT.ModifyFaceGrad();
}

Real Iutherland::Icdim = 110.4;
Real Iutherland::Ic;

Iutherland::Iutherland()
{
}

Iutherland::~Iutherland()
{

}

void Iutherland::ICalcConst()
{
	Iutherland::Ic = Iutherland::Icdim / nscom.tref_dim;
}

Real Iutherland::ICalcViscosity( Real t )
{
	Real t3 = t * t * t;
	return sqrt(t3) * (1 + Iutherland::Ic) / ( t + Iutherland::Ic );
}




EndNameSpace
