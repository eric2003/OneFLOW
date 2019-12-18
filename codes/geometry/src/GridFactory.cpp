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

#include "GridFactory.h"
#include "CgnsFactory.h"
#include "GridPara.h"
#include "GridMediator.h"
#include "DomainInp.h"
#include "Su2Grid.h"
#include "DataBase.h"
#include "ClassicGrid.h"
#include "StrGrid.h"
#include "PointSearch.h"
#include "BcRecord.h"
#include "HXMath.h"
#include "Partition.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

void GenerateGrid()
{
    GridFactory * gf = new GridFactory();
    gf->Run();
    delete gf;
}

GridFactory::GridFactory()
{
}

GridFactory::~GridFactory()
{
}

void GridFactory::Run()
{
    grid_para.Init();

    switch (grid_para.gridObj)
    {
    case 0: //����һЩ�������ε������緽ǻ��Բ����RAE2822���͵ȵ�
        this->DataBaseGrid();
        break;
    case 1:    //ת������
        this->ConvertGrid();
        break;
    case 2:
        this->GeneInp();
        break;
    case 3:    //�������
        this->PartGrid();
        break;
    default:
        break;
    }
}

void GridFactory::GeneInp()
{
    DomainInp * domainInp = new DomainInp();
    domainInp->Run();
    delete domainInp;
}

void GridFactory::PartGrid()
{
    Partition * part = new Partition();
    part->Run();
    delete part;
}

void GridFactory::ConvertGrid()
{
    string sourceGridType = grid_para.filetype; 
    if ( sourceGridType == "plot3d" )
    {
        this->Plot3DProcess();
    }
    else if ( sourceGridType == "su2" )
    {
        this->SU2Process();
    }
    else if ( sourceGridType == "cgns" )
    {
        this->CGNSProcess();
    }
}

void GridFactory::DataBaseGrid()
{
    ClassicGrid * classicGrid = new ClassicGrid();
    classicGrid->Run();
    delete classicGrid;
}

void GridFactory::Plot3DProcess()
{
    CgnsFactory * cgnsFactory = new CgnsFactory();

    cgnsFactory->CommonToOneFlowGrid();

    delete cgnsFactory;
}

void GridFactory::SU2Process()
{
    Su2Grid * su2Grid = new Su2Grid();
    su2Grid->Su2ToOneFlowGrid();
    delete su2Grid;
}

void GridFactory::CGNSProcess()
{
    CgnsFactory * cgnsFactory = new CgnsFactory();

    cgnsFactory->GenerateGrid();

    delete cgnsFactory;
}

EndNameSpace