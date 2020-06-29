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
#include "UsdBasic.h"

BeginNameSpace( ONEFLOW )

class UsdData : public UsdBasic
{
public:
    UsdData();
    ~UsdData();
public:
    int nEqu;

    Real vol, vol1, vol2;
    RealField res, res0, res1, res2;
    RealField prim, prim1, prim2;
    RealField q, q1, q2;
    RealField dualtimeRes;
    RealField dualtimeSrc;
public:
    RealField normList;
public:
    Real sum1, sum2, norm0, totalNorm;
    Real conv;
    int  iConv;
public:
    void CalcCellDualTimeResidual();
    void CalcCellDualTimeSrc();
public:
    virtual void Init();
    void InitSub( int nEqu );
    void ZeroData();
    void CalcCellUnsteadyCri();
    void CalcCvg();
};

extern UsdBasic usd;

EndNameSpace