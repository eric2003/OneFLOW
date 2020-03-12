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
#include "SimuCtrl.h"
#include "FileUtil.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

bool SimuCtrl::hx_debug = false;
bool SimuCtrl::run_from_ide = false;
string SimuCtrl::system_root = "";

SimuCtrl::SimuCtrl()
{

}

SimuCtrl::~SimuCtrl()
{
}

void SimuCtrl::Init()
{
    string exePath = HX_GetExePath();
    cout << " exe path = " << exePath << "\n";
    string local_root = "/system/";
    if ( SimuCtrl::run_from_ide )
    {
        string curr_dir = HX_GetCurrentDir();
        cout << " curr_dir = " << curr_dir << "\n";
        SimuCtrl::system_root = curr_dir + local_root;
    }
    else
    {
        SimuCtrl::system_root = exePath + local_root;
    }
}

EndNameSpace
