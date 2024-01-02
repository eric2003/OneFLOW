/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2024 He Xin and the OneFLOW contributors.
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


BeginNameSpace( ONEFLOW )

//Mesh parameters
class GridPara
{
public:
    GridPara();
    ~GridPara();
public:
    std::string topo;
    std::string filetype; //plot3d, cgns...
    std::string target_filetype;
    std::string format; //binary, ascii
    std::string gridFile;
    std::string bcFile;
    std::string targetFile;

    //Conversion operations performed on the grid
    int gridObj;    

    //Is it a multiblock mesh
    int multiBlock;

    int axis_dir;
    //Mesh scaling factor
    Real gridScale;
    RealField gridTrans;
public:
    //Initialize mesh parameters
    void Init();
};

extern GridPara grid_para;

int GetGridTopoType();

EndNameSpace
