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
#include "Configure.h"
#include "HXDefine.h"
BeginNameSpace( ONEFLOW )

typedef void ( * TIME_INTEGRAL )( void );

const int MULTI_STAGE = 1;
const int LUSGS = 2;
const int SIMPLE = 3;


class SweepState
{
public:
    SweepState ();
    ~SweepState();
public:
    static int nSweeps;
};

class TimeIntegral
{
public:
    TimeIntegral ();
    ~TimeIntegral();
public:
    static TIME_INTEGRAL timeIntegral;
public:
    static void Init();
    static void Relaxation( int nCycles );
public:
    static void RungeKutta();
    static void Lusgs();
	static void Simple();
};

EndNameSpace