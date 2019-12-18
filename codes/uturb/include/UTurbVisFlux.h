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
    void CmpVisFlux();
    void CmpVisFlux1Equ();
    void CmpVisFlux2Equ();

    void PrepareFaceValue();
    void CmpFaceVisFlux1Equ();
    void CmpFaceVisFlux2Equ();
    void UpdateFaceVisFlux();
public:
    void CmpFaceWeight();
    void CmpGradCoef();
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
    void CmpNormalGrad();
    void CmpTestMethod();
    void CmpNew1Method();
    void CmpNew2Method();
    void CorrectFaceGrad();
    void ModifyFaceGrad();
    void AccurateFaceValue();
    void AccurateOtherFaceValue();
};

extern VisGrad visTurb;


EndNameSpace