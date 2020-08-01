/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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

void NsCalcBc();
void NsCalcGamaT(int flag);
void NsCalcRHS();
void NsCalcInvFlux();
void NsCalcVisFlux();
void NsCalcSrcFlux();
void NsCalcChemSrc();
void NsCalcTurbEnergy();
void NsCalcDualTimeStepSrc();

void INsCmpBc();
void INSCmpGamaT(int flag);
void INsCmpRHS();
void INsCmpChemSrc();
void INsCmpTurbEnergy();
void INsCorrectPresscoef();
void INsCmpTimestep();
void INsPreflux();
void INsCmpInv();
void INsCmpVis();
void INsCmpUnstead();
void INsCmpSrc();
void INsMomPre();
void INsCmpFaceflux();
void INsCmpPressCorrectEquandUpdatePress();
void INsUpdateFaceflux();
void INsUpdateRes();
void INsCmpSpeedCorrectandUpdateSpeed();


EndNameSpace