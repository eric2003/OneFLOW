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
#include "INsInvterm.h"
#include "systemSolver.h"
#include "poisson.h"

BeginNameSpace(ONEFLOW)

class UINsFField;
class Limiter;
class LimField;
class SolveMRhs;


class UINsInvterm : public INsInvterm
{
public:
    UINsInvterm();
    ~UINsInvterm();
public:
    void Alloc();
    void DeAlloc();
	void CmpINsTimestep();
	void CmpINsPreflux();
	void INsPreflux();
	void Initflux();
    void CmpInvcoff();
    void CmpInvMassFlux();
    void CmpInvFace();
    void CmpLimiter();
	void CmpFaceflux();
	void CmpINsMomRes();
	void CmpINsPreRes();
	void CmpCorrectPresscoef();
	void CmpNewMomCoe();
	void CmpPressCorrectEqu();
	void UpdateFaceflux();
	void CmpUpdateINsFaceflux();
	void CmpUpdateINsBcFaceflux();
	void UpdateSpeed();
	void UpdateINsRes();
    void AddFlux();
    void PrepareFaceValue();
	void PrepareProFaceValue();
	void CmpPreGrad();
	//void CmpINsinvTerm();
    //void UpdateFaceInvFlux();
    void ReadTmp();
public:
    void GetQlQrField();
    void ReconstructFaceValueField();
    void BoundaryQlQrFixField();
    void Init();
    void MomPre();
public:
    Limiter* limiter;
    LimField* limf;
    MRField* iinvflux;
public:
    Real Number;
};
//void PrimToQ(RealField & prim, Real gama, RealField & q);
extern UINsInvterm NonZero;
EndNameSpace