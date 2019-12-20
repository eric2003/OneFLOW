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
#include "NsVisFlux.h"
#include "HXArray.h"

BeginNameSpace( ONEFLOW )

class UNsVisFlux : public NsVisFlux
{
public:
    UNsVisFlux ();
    ~UNsVisFlux();
public:
    typedef void ( UNsVisFlux:: * VisPointer )();
    VisPointer visPointer;
    MRField * visflux;
public:
    void SetVisPointer();
    void CmpFlux();
    void PrepareField();
    void CmpVisFlux();
    void AddVisFlux();
    void CmpFaceVisFlux();
    void UpdateFaceVisFlux();
    void CmpHeatFlux();
    void CmpStress();
    void CmpAniStress();
    void CmpNsVisFlux();
    void ZeroHeatFlux();
    void AddChemHeatFlux();
    void AddHeatFlux();
    void SaveHeatFlux();

    void Alloc();
    void DeAlloc();
public:
    void PrepareFaceValue();
    void SaveFacePara();
    void CmpFaceWeight();
public:
    void AverMethod();
    void StdMethod();
    void TestMethod();
    void New1Method();
    void New2Method();
    void CmpGradCoef();
    void PrepareCellGeom();
};

void CmpLaminarViscosity( int flag );



EndNameSpace