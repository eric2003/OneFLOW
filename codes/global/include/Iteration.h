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

class Iteration
{
public:
    Iteration();
    ~Iteration();
public:
    static int innerSteps;
    static int outerSteps;
    static int maxSteps;
    static int dualtime;
    static int nFieldSave;
    static int nVisualSave;
    static int nResSave;
    static int nForceSave;
    static Real cfl;
    static Real cflst;
    static Real cfled;
    static int ncfl;
public:
    static void Init();
    static bool InnerOk();
    static bool ResOk();
    static bool ForceOk();
};

class SimuIterState
{
public:
    SimuIterState();
    ~SimuIterState();
public:
    static bool Running();
    static bool InnerRunning();
};

class SaveState
{
public:
    SaveState();
    ~SaveState();
public:
    static int nFileSaveSteps;
};

EndNameSpace