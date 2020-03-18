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
#include <map>
#include <string>
using namespace std;

BeginNameSpace( ONEFLOW )

const int CFD_SOLVER         = 0;
const int NS_SOLVER          = 1;
const int TURB_SOLVER        = 2;
const int INC_NS_SOLVER      = 3;
const int HIGH_ORDER_SOLVER  = 4;
const int SCALAR_SOLVER      = 5;
const int GRID_SOLVER        = 6;
const int SIMU_SOLVER        = 7;
const int SOLVER_CAPACITY    = 8;

const int INTERFACE_DATA          = 0;
const int INTERFACE_DQ_DATA       = 1;
const int INTERFACE_GRADIENT_DATA = 2;
const int INTERFACE_OVERSET_DATA  = 3;

const int SOLVER_BASED = 0;
const int GRID_BASED   = 1;

const int NO_DATA   = 0;
const int WITH_DATA = 1;

const int FIELD_FLOW = 0;
const int FIELD_RHS  = 1;

const int GREAT_SEND = 0;
const int GREAT_RECV = 1;

const int SEND_STORAGE = 0;
const int RECV_STORAGE = 1;

const int GREAT_READ  = 0;
const int GREAT_WRITE = 1;
const int GREAT_ZERO  = 2;

int GetOppositeSendRecv( int iSr );

extern map< string, int > * solverTypeMap;
extern map< string, int > * interfaceMap;
extern map< string, int > * sendRecvMap;
extern map< string, int > * fieldIdMap;

void CreateSysMap();

EndNameSpace