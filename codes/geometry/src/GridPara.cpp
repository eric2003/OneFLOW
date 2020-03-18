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

#include "GridPara.h"
#include "Ctrl.h"
#include "DataBase.h"
#include <string>
using namespace std;

BeginNameSpace( ONEFLOW )

GridPara grid_para;

GridPara::GridPara()
{
    ;
}

GridPara::~GridPara()
{
    ;
}

void GridPara::Init()
{
    //设置原始网格的文件名称
    this->gridFile = GetDataValue< string >("sourceGridFileName");
    this->bcFile   = GetDataValue< string >("sourceGridBcName");
    //set target grid file name
    this->targetFile = ONEFLOW::GetDataValue< string >( "targetGridFileName" );
    //设置原始网格格式
    this->filetype = GetDataValue< string >("sourceGridType");
    //设置target grid type
    this->target_filetype = GetDataValue< string >("targetGridType");
    //设置原始网格的拓扑形式
    this->topo = GetDataValue< string >("topoType");

    this->multiBlock = GetDataValue< int >( "multiBlock" );
    //设置要进行的网格操作
    this->gridObj = GetDataValue< int >("gridObj");
    //设置网格缩放比例
    this->gridScale =  GetDataValue< Real >( "gridScale" );
    //设置网格平移量
    this->gridTrans.resize( 3 );
    CopyArray( this->gridTrans, "gridTrans" );

    this->axis_dir = GetDataValue< int >( "axis_dir" );
}

int GetGridTopoType()
{
    return 0;
}

EndNameSpace