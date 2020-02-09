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

#include "UNsSolver.h"
#include "Mesh.h"
#include "FaceMesh.h"
#include "CellMesh.h"
#include "SolverInfo.h"
#include "SolverState.h"
#include "FaceTopo.h"
#include "Boundary.h"
#include "BcRecord.h"
#include "DataBase.h"
#include "NsIdx.h"
#include "HXMath.h"
#include "UNsLusgs.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

REGISTER_SOLVER( UNsSolver )

UNsSolver::UNsSolver()
{
}

UNsSolver::~UNsSolver()
{
}

void UNsSolver::StaticInit()
{
    NsSolver::StaticInit();
    LusgsState::AddSolver( this->sid, this->gridType, new UNsLusgs() );
}

void UNsSolver::Init()
{
}

void UNsSolver::Run()
{
}

EndNameSpace