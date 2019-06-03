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
    this->multiBlock = GetDataValue< int >( "multiBlock" );
    this->gridObj = ONEFLOW::GetGridObj();
	this->gridScale =  GetDataValue< Real >( "gridScale" );
	this->gridTrans.resize( 3 );
	CopyArray( this->gridTrans, "gridTrans" );
	this->axis_dir = GetDataValue< int >( "axis_dir" );
}

int GetGridObj()
{
    return ONEFLOW::GetDataValue< int >( "gridObj" );
}

string GetSourceGridType()
{
    return ONEFLOW::GetDataValue< string >( "sourceGridType" );
}

string GetTopoType()
{
    return ONEFLOW::GetDataValue< string >( "topoType" );
}

int GetGridTopoType()
{
    return 0;
}

EndNameSpace