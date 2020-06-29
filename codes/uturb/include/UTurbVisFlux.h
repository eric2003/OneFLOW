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
#include "TurbVisFlux.h"
#include "VisGrad.h"

BeginNameSpace( ONEFLOW )

class UTurbVisFlux : public TurbVisFlux
{
public:
    UTurbVisFlux ();
    ~UTurbVisFlux();
public:
    typedef void ( UTurbVisFlux:: * VisPointer )();
    VisPointer visPointer;
    MRField * visflux;
public:
    void Alloc();
    void DeAlloc();
    void AddVisFlux();
    void SetVisPointer();
    void CalcVisFlux();
    void CalcVisFlux1Equ();
    void CalcVisFlux2Equ();

    void PrepareFaceValue();
    void CalcFaceVisFlux1Equ();
    void CalcFaceVisFlux2Equ();
    void UpdateFaceVisFlux();
public:
    void CalcFaceWeight();
    void CalcGradCoef();
    void PrepareCellGeom();
public:
    void AverMethod();
    void StdMethod();
    void TestMethod();
    void New1Method();
    void New2Method();
public:
    void ZeroNormalGrad();
    void AverGrad();
    void AverFaceValue();
    void AverOtherFaceValue();
    void CalcNormalGrad();
    void CalcTestMethod();
    void CalcNew1Method();
    void CalcNew2Method();
    void CorrectFaceGrad();
    void ModifyFaceGrad();
    void AccurateFaceValue();
    void AccurateOtherFaceValue();
};

extern VisGrad visTurb;


EndNameSpace