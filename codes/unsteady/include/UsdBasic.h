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

class UsdBasic
{
public:
    UsdBasic();
    ~UsdBasic();
public:
    Real bsc1, bsc2, bsc3;
    Real sc1, sc2, sc3;
    Real sp1, sp2;
    Real resc1, resc2, resc3;
    RealField coeff;
public:
    void InitCoef();
    void CalcResCoef();
    void CalcSpectrumCoeff();
    void CalcSrcCoeffBasic();
    void CalcSrcCoeff();
public:
    void InitBasic();
};

EndNameSpace