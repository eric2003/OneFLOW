/*---------------------------------------------------------------------------*\
OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
Copyright (C) 2017-2022 He Xin and the OneFLOW contributors.
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
#include "Simu.h"
#include "Geom.h"
#include "Solver.h"
#include "Cmpi.h"
#include "CfdPara.h"
#include "Project.h"
#include <iostream>

Simu::Simu(int argc, char **argv)
{
    Project::Init( argc, argv );
    Cmpi::Init( argc, argv );
}

Simu::~Simu()
{
    Cmpi::Finalize();
}

void Simu::Init(int argc, char **argv)
{
}

void Simu::Run()
{
    Geom_t::Init();
    Geom * geom = new Geom{};
    geom->Init();
    geom->GenerateGrid();
    geom->ComputeGeom();

    //cfd parameter
    CfdPara * cfd_para = new CfdPara{};
    cfd_para->Init( geom );

    Solver * solver = new Solver{};
    solver->Run( cfd_para, geom );
    delete cfd_para;
    delete geom;
    delete solver;
    Geom_t::Finalize();
}