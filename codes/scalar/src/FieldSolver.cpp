/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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

#include "FieldSolver.h"
#include "FieldPara.h"
#include "ScalarGrid.h"
#include "MetisGrid.h"
#include "ScalarField.h"
#include <iostream>
#include <vector>
using namespace std;

BeginNameSpace( ONEFLOW )

FieldSolver::FieldSolver()
{
    this->grid = new ScalarGrid();
    this->field = new ScalarField();
    this->para = new FieldPara();
}

FieldSolver::~FieldSolver()
{
    delete this->grid;
    delete this->field;
    delete this->para;
}

void FieldSolver::Run()
{
    this->Init();

    this->SolveFlowField();
}

void FieldSolver::Init()
{
    this->InitCtrlParameter();

    this->InitGrid();

    this->InitFlowField();
}

void FieldSolver::InitCtrlParameter()
{
    this->para->Init();
}

void FieldSolver::InitGrid()
{
    this->grid->GenerateGrid( this->para->nx, 0, this->para->len );
    this->grid->CalcTopology();
    this->grid->CalcMetrics1D();

    Part part;
    NetGrid netGrid;
    part.PartitionGrid( this->grid, 4, & netGrid );
    int kkk = 1;
}

void FieldSolver::InitFlowField()
{
    field->InitFlowField( this->grid );
}

void FieldSolver::SolveFlowField()
{
    field->SolveFlowField( this->para );
}

EndNameSpace