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
#include "HXArray.h"

BeginNameSpace( ONEFLOW )

const int ILMT_ZERO   = -1;
const int ILMT_NO     = 0;
const int ILMT_BARTH  = 1;
const int ILMT_VENCAT = 2;

Real BarthFunction( Real dq, Real minv, Real maxv, Real dot );
Real VencatFunction( Real dq, Real minv, Real maxv, Real eps );
Real VencatFunctionEric( Real dq, Real minv, Real maxv, Real eps );
Real Vencat( Real x, Real y, Real eps );
Real VencatEric( Real x, Real y, Real eps );

class Lim
{
public:
    Lim();
    ~Lim();
public:
    RealField *q, *dqdx, *dqdy, *dqdz;
    RealField *limiter;
    RealField * minvf, * maxvf;
    Real minv1, minv2, maxv1, maxv2;
    Real dqdx1, dqdy1, dqdz1;
    Real dqdx2, dqdy2, dqdz2;
    Real lim1, lim2;
    Real qmin, qmax;
};

typedef bool ( * CheckFun )( RealField & );

class LimField
{
public:
    LimField();
    virtual ~LimField();
public:
    virtual void Init(){};
    Real ModifyLimiter( Real phil, Real phir );
    void CalcFaceValue();
    void CalcFaceValueWeighted();
    void GetQlQr();
    virtual void BcQlQrFix();
public:
    int nEqu;
    MRField * q;
    MRField * dqdx, * dqdy, * dqdz;
    MRField * limiter;

    MRField * qf1, * qf2;
    CheckFun ckfun;
};


class Limiter
{
public:
    Limiter();
    virtual ~Limiter();
public:
    Lim * lim;
    LimField * limf;
    int limflag;
public:
    void Alloc();
    void DeAlloc();
    void SetInitValue();
    void CalcLimiter();
    void CalcLimiterScalar();

    void CalcZeroLimiter();
    void CalcNoLimiter();
    void CalcBarthLimiter();
    void CalcVencatLimiter();
    void CalcLocalBarthLimiter();
    void CalcLocalVencatLimiter();
    void PrepareData();
    void CalcMinMaxDiff();
};

bool NoCheck( RealField & q );

EndNameSpace