/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/


#pragma once
#include "HXDefine.h"
#include "UnsteadyConvergence.h"

BeginNameSpace( ONEFLOW )

class UsdData
{
public:
    UsdData();
    ~UsdData();
public:
    int nEqu;

    RealField res, res1, res2;
    RealField dualtimeRes;
    RealField dualtimeSrc;
public:
    UnsteadyConvergence convergence;
public:
    const RealField & GetQ() const;
    const RealField & GetQ1() const;
    const RealField & GetQ2() const;

    void CalcCellDualTimeResidual();
    void CalcCellDualTimeSrc( Real vol, Real vol1, Real vol2 );
public:
    virtual void Init();
    void InitSub( int nEqu );
};

EndNameSpace
