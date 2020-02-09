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

BeginNameSpace( ONEFLOW )

class Ctrl
{
public:
    Ctrl();
    ~Ctrl();
public:
    int startStrategy;
    int ireadwdst;
    int iConv;
    Real lhscoef;
    Real pdt, pdt1;
    int linearTwoStepMethods;
    int showfield;
    int idualtime;
    int time_integral;
    int ieigenfix;
    Real centropy1;
    Real centropy2;
    int idump;
    int inflowType;
    RealField initplane;
    RealField initflow1;
    RealField initflow2;
    int rk_stage;
    RealField rk_coef;
    Real maxTime;
    Real currTime;
    int iexitflag;
    int ilim;
    Real vencat_coef;
    int nrokplus;
    int ivischeme;
    string heatfluxFile;
public:
    void Init();
};

class VisModel
{
public:
    VisModel();
    ~VisModel();
public:
    void Init();
public:
    int vismodel;
    string visname;
};

void CopyArray( RealField & arr, const string & name );

extern Ctrl ctrl;
extern VisModel vis_model;

EndNameSpace