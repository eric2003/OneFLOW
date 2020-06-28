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

const int VIS_STD  = 0;
const int VIS_AVER = 1;
const int VIS_TEST = 2;
const int VIS_NEW1 = 3;
const int VIS_NEW2 = 4;

class VisGradGeom
{
public:
    VisGradGeom();
    ~VisGradGeom();
public:
    void CalcFaceWeight();
    void CalcAngle( Real dx, Real dy, Real dz, Real dist, Real & angle );
    void PrepareCellGeom();
    void CalcGradCoef();
public:
    Real dxl, dyl, dzl;
    Real dxr, dyr, dzr;
    Real delta, delt1, delt2;

    Real dxnl, dynl, dznl;
    Real dxnr, dynr, dznr;

    Real fw1, fw2;
    Real c1, c2;
    Real d1, d2, d, od;
    Real angle1, angle2;
    Real skewAngle;

    Real dx, dy, dz, ods;
};

extern VisGradGeom vgg;

class VisGrad
{
public:
    VisGrad ();
    ~VisGrad();
public:
    int nEqu;
    void Init( int nEqu );
public:
    void AverGrad();
    void ZeroNormalGrad();
    void AverFaceValue ();
    void CorrectFaceGrad();
    void CalcNormalGrad();
    void CalcTestMethod();
    void CalcNew1Method();
    void CalcNew2Method();
    bool FaceAngleIsValid();
    bool TestSatisfied();
    bool New1Satisfied();
    bool New2Satisfied();
    void CalcC1C2();
    void AccurateSideValue();
    void AccurateFaceValue();
    void ModifyFaceGrad();
public:
    RealField q, q1, q2;
    RealField q11, q22;

    RealField dqdx, dqdy, dqdz, dqdn;
    RealField dqdx1, dqdy1, dqdz1;
    RealField dqdx2, dqdy2, dqdz2;
    RealField dqdn1, dqdn2;
    RealField dqdt1, dqdt2;
};

void CorrectGrad( Real fl, Real fr, Real & dfdx, Real & dfdy, Real & dfdz, Real dx, Real dy, Real dz, Real ods );

EndNameSpace