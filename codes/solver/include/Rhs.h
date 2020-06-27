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


#pragma once
#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

class Rhs
{
public:
    Rhs ();
    ~Rhs();
public:
    void UpdateNsResiduals();
	void UpdateINsResiduals();
};

void NsCompBc();
void NsCompGamaT( int flag );
void NsCompRHS();
void NsCompInvFlux();
void NsCompVisFlux();
void NsCompSrcFlux();
void NsCompChemSrc();
void NsCompTurbEnergy();
void NsCompDualTimeStepSrc();

void INsCompBc();
void INsCompGamaT(int flag);
void INsCompRHS();
//void INsCompInvFlux();
//void INsCompVisFlux();
//void INsCompSrcFlux();
void INsCompChemSrc();
void INsCompTurbEnergy();
//void INsCompDualTimeStepSrc();
void INsCorrectPresscoef();
//void INsCorrectSpeed();
void INsCompInv();
void INsCompVis();
void INsCompSrc();
void INsMomPred();
void INsCompFaceflux();
void INsCompPressCorrectEquandUpdatePress();
void INsUpdateFaceflux();
void INsCompSpeedCorrectandUpdateSpeed();


EndNameSpace