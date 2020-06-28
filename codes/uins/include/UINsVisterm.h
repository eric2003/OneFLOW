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
#include "INsVisterm.h"
#include "HXArray.h"

BeginNameSpace( ONEFLOW )

class UINsVisterm : public INsVisterm
{
public:
	UINsVisterm();
    ~UINsVisterm();
public:
    typedef void ( UINsVisterm:: * VisPointer )();
    VisPointer visPointer;
    MRField * visflux;
public:
    void CalcViscoff();
    void PrepareField();
    void CalcNsVisterm();
	void CalcVisterm();
	void CalcFaceVisterm();

    void Alloc();
    void DeAlloc();
	void CalcSrc();
public:
    void PrepareFaceValue();
    void SaveFacePara();
    void CalcFaceWeight();
public:
    void CalcGradCoef();
    void PrepareCellGeom();
};

void ICalcLaminarViscosity( int flag );



EndNameSpace