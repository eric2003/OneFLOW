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
    //Sets the file name of the original grid
    this->gridFile = GetDataValue< string >("sourceGridFileName");
    this->bcFile   = GetDataValue< string >("sourceGridBcName");
    //set target grid file name
    this->targetFile = ONEFLOW::GetDataValue< string >( "targetGridFileName" );
    //Format original grid
    this->filetype = GetDataValue< string >("sourceGridType");
    //Set target grid type
    this->target_filetype = GetDataValue< string >("targetGridType");
    //Sets the topology of the original mesh
    this->topo = GetDataValue< string >("topoType");

    this->multiBlock = GetDataValue< int >( "multiBlock" );
    //Set the grid operation to be performed
    this->gridObj = GetDataValue< int >("gridObj");
    //Set mesh scale
    this->gridScale =  GetDataValue< Real >( "gridScale" );
    //Set mesh translation amount
    this->gridTrans.resize( 3 );
    CopyArray( this->gridTrans, "gridTrans" );

    this->axis_dir = GetDataValue< int >( "axis_dir" );
}

int GetGridTopoType()
{
    return 0;
}

EndNameSpace