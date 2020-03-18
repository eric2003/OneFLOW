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
#include "VisGrad.h"
BeginNameSpace( ONEFLOW )

class INsVis
{
public:
    INsVis();
    ~INsVis();
public:
    void Init();
public:
    Real dudx, dudy, dudz;
    Real dvdx, dvdy, dvdz;
    Real dwdx, dwdy, dwdz;
	Real dpdx, dpdy, dpdz;
    Real dtdx, dtdy, dtdz;
    Real dtdn;

    Real txx, tyy, tzz;
    Real txy, txz, tyz;

    Real rhok;

    Real b11, b22, b33;
    Real b12, b13, b23;

    Real qx, qy, qz, qNormal;
    Real tmid;
    RealField fvis;
	Real p1, p2;
    Real um, vm, wm;
};

extern INsVis Ivis;

class INsVisterm
{
public:
	INsVisterm();
    ~INsVisterm();
public:
    void AverGrad();
    void ZeroNormalGrad();
    void AverFaceValue();
    void AverOtherFaceValue();
    void AccurateFaceValue();
    void AccurateOtherFaceValue();
    void CorrectFaceGrad();
    void CmpNormalGrad();
    void CmpTestMethod();
    void CmpNew1Method();
    void CmpNew2Method();
    void ModifyFaceGrad();
};

extern VisGrad visQ;
extern VisGrad visT;

class Iutherland
{
public:
	Iutherland();
	~Iutherland();
public:
	static Real Icdim;
	static Real Ic;
	static void ICmpConst();
	static Real ICmpViscosity(Real t);
};

EndNameSpace
