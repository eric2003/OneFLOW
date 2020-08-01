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
#include "Com.h"

BeginNameSpace( ONEFLOW )

const int _INVISCID  = 0;
const int _LAMINAR   = 1;
const int _ALGEBRAIC = 2;
const int _ONE_EQU   = 3;
const int _TWO_EQU   = 4;

class INsCom
{
public:
    INsCom();
    ~INsCom();
public:
    void Init();
public:
    bool init_flag;
    int nBEqu;
    int nEqu;
    int nTEqu;
    int nSpecies;
    int visSRModel;
    int timestepModel;
    int visTimestepModel;
    int chemModel;
    int nTModel;
    int bctype;
    int bcdtkey;
    int icmpInv;
    int nProbe;
    int ischeme;
    int ivischeme;
public:
    int faceOuterNormal;
public:
    Real vis;
    Real visl, visl1, visl2;
    Real vist, vist1, vist2;
    Real prl, prt;
    Real oprl, oprt;
    Real const_cp, kcp;
    Real reynolds, oreynolds;
public:
    Real mach_ref;
    Real gama_ref;
    Real gama;
    Real gama1;
    Real gama2;
public:
    Real twall_dim;
    Real twall;
    Real statecoef;

    Real elevation;
    Real tref_dim;
    Real pref_dim;
    Real dref_dim;
    Real vref_dim;
    Real cref_dim;
    Real tref;
    Real pref;
    Real dref;
    Real vref;
    Real cref;
    int gasInfoStrategy;
    int machStrategy;
    string gasModelFile;
    Real schmidtl;
    Real schmidtt;
    Real reylref_dim;
    Real aoa, aos;
    RealField refns;
    Real visref_dim;
    Real csuth;
    Real csuth_dim;
    Real dim_amw; //dimensional average molecular weight
    Real amw; //average molecular weight
public:
    RealField q1;
    RealField q2;

    RealField q;
    RealField q0;
    RealField dq;

    RealField prim;
    RealField prim0;

    RealField t;
    RealField t0;

    RealField prims1;
    RealField prims2;

    RealField primt1;
    RealField primt2;

    RealField ts1;
    RealField ts2;

    RealField tt1;
    RealField tt2;

    RealField inflow;
    RealField * bcflow;
public:
    Real invsr; //inviscid spectrum radius;
    Real vissr; //viscous  spectrum radius;
    Real dt;
    Real timestep;
    Real minTimestep;
    Real max_time_ratio;
    Real physicalTimestep;
};

extern INsCom inscom;

void INsExtract( RealField & prim, Real & rm, Real & um, Real & vm, Real & wm, Real & pm );

bool INsCheckFunction( RealField & q );

EndNameSpace