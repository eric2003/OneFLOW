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
string SimuCtrl::execute_dir = "";
string SimuCtrl::current_dir = "";

SimuCtrl::SimuCtrl()
{

}

SimuCtrl::~SimuCtrl()
{
}

void SimuCtrl::Init()
{
    SimuCtrl::execute_dir = HX_GetExePath();
    SimuCtrl::current_dir = HX_GetCurrentDir();

    cout << " SimuCtrl::execute_dir = " << SimuCtrl::execute_dir << "\n";
    cout << " SimuCtrl::current_dir = " << SimuCtrl::current_dir << "\n";

    string local_root = "/system/";
    if ( SimuCtrl::run_from_ide )
    {
        string current_dir_now = RemoveEndSlash( SimuCtrl::current_dir );
        SimuCtrl::system_root = current_dir_now + local_root;
    }
    else
    {
        string execute_dir = RemoveEndSlash( SimuCtrl::execute_dir );
        SimuCtrl::system_root = SimuCtrl::execute_dir + local_root;
    }
    cout << "SimuCtrl::system_root = " << SimuCtrl::system_root << "\n";
}

EndNameSpace
