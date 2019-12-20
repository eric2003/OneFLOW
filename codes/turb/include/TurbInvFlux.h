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
#include "HXArray.h"
BeginNameSpace( ONEFLOW )

class TurbInv
{
public:
    TurbInv();
    ~TurbInv();
public:
    void Init();
public:
    RealField prim1, prim2;
    RealField q1, q2;
    RealField flux;
public:
    Real rl, ul, vl, wl;
    Real rr, ur, vr, wr;
    Real rho_coef;
};

extern TurbInv turbInv;

class TurbInvFlux
{
public:
    TurbInvFlux ();
    ~TurbInvFlux();
public:
    void RoeFlux();
};

EndNameSpace