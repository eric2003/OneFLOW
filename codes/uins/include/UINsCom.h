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
#include "INsCom.h"

BeginNameSpace( ONEFLOW )

class BcRecord;

class UINsField
{
public:
    UINsField();
    ~UINsField();
public:
    void Init();
public:
    MRField * q, * q1, * q2;
    MRField * dq;
    MRField * rhs, * drhs;
    MRField * gama, * gama1, * gama2;
    MRField * dqdx, * dqdy, * dqdz;
    MRField * dtdx, * dtdy, * dtdz;
    MRField * bc_q;
    MRField * bcdqdx, * bcdqdy, * bcdqdz;
    MRField * invsr; //inviscid spectrum radius;
    MRField * vissr; //viscous  spectrum radius;
    MRField * impsr; //implicit spectrum radius;
    MRField * visl;
    MRField * vist;
    MRField * timestep;
    MRField * tempr;
    MRField * invflux;
    MRField * visflux;
    MRField * limiter;
    MRField * res;
    MRField * res1;
    MRField * res2;
    RealField * spectrum;
};

extern UINsField uinsf;

EndNameSpace