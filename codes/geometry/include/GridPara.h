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
#include <string>
using namespace std;

BeginNameSpace( ONEFLOW )

//网格参数
class GridPara
{
public:
    GridPara();
    ~GridPara();
public:
    string topo;
    string filetype; //plot3d, cgns...
    string target_filetype;
    string format; //binary, ascii
    string gridFile;
    string bcFile;
    string targetFile;

    //对网格执行的转换操作
    int gridObj;    

    //是否是多块网格
    int multiBlock;

    int axis_dir;
    //网格缩放因子
    Real gridScale;
    RealField gridTrans;
public:
    //初始化网格参数
    void Init();
};

extern GridPara grid_para;

int GetGridTopoType();

EndNameSpace