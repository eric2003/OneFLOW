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

BeginNameSpace( ONEFLOW )
class TurbSpecData
{
public:
    TurbSpecData();
    ~TurbSpecData();
public:
    void Init();
public:
    int  nEqu;

    RealField matrix1, matrix2;
    RealField radius1, radius2;
    RealField work;
    RealField q1, q2;
    RealField q1_ns, q2_ns;
    Real rl, ul, vl, wl;
    Real rr, ur, vr, wr;
};

extern TurbSpecData turbsp;

class TurbSpectrum
{
public:
    TurbSpectrum();
    ~TurbSpectrum();
public:
    void CmpFaceSpectrum1Equ();
    void CmpFaceSpectrum2Equ();
};


EndNameSpace