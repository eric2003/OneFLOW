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

class UTurbField
{
public:
    UTurbField();
    ~UTurbField();
public:
    void Init();
    MRField * q, * q1, * q2;
    MRField * dq;
    MRField * rhs, * drhs;
    MRField * visl, * vist;
    MRField * dqdx, * dqdy, * dqdz;
    MRField * q_ns, * dqdx_ns, * dqdy_ns, * dqdz_ns;
    MRField * bld;
    MRField * res;
    MRField * impsr;
    MRField * cross;
    MRField * len_scale;
    MRField * timestep;
    MRField * matrix_l, * matrix_r;
    RealField * dist;
};

extern UTurbField uturbf;


EndNameSpace