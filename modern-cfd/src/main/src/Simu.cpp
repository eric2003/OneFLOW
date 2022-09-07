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
#include "tools.h"
#include <json/json.h>
#include <iostream>
#include <fstream>


FieldSolver_t::FieldSolver_t()
{
    this->cfd_para = new CfdPara{};
}

FieldSolver_t::~FieldSolver_t()
{
    delete this->cfd_para;
}


void FieldSolver_t::Init( Json::Value & root )
{
    this->cfd_para->Init( root );
}

void FieldSolver_t::Run()
{
    Geom_t::Init();
    Geom_t::LoadGrid();
    Geom * geom = new Geom{};
    geom->Init();
    geom->GenerateGrid();
    geom->ComputeGeom();

    Solver * solver = new Solver{};
    solver->Run( this->cfd_para, geom );
    delete cfd_para;
    delete geom;
    delete solver;
    Geom_t::Finalize();
}

void GridSolver_t::Init( Json::Value & root )
{
    Json::Value item = root[ "gridgen" ];
    std::string gridobj_str = item[ "gridobj" ].asString();
    Geom_t::gridName = item[ "gridfile" ].asString();
    if ( gridobj_str == "read" )
    {
        Geom_t::gridobj = 1;
    }
    else if ( gridobj_str == "modify" )
    {
        Geom_t::gridobj = 2;
    }
    else
    {
        Geom_t::gridobj = 0;
    }
}

void GridSolver_t::Run()
{
    Geom_t::Init();
    Geom_t::LoadGrid();
    Geom_t::Finalize();
}

Simu::Simu(int argc, char **argv)
{
    this->solver_t = 0;
    Cmpi::Init( argc, argv );
    Project::Init( argc, argv );
}

Simu::~Simu()
{
    Cmpi::Finalize();
    delete this->solver_t;
}

void Simu::Init(int argc, char **argv)
{
}

void Simu::ReadControlParameter()
{
    Cmpi::server_out << " Simu::ReadControlParameter() Project::prj_grid_dir = " << Project::prj_grid_dir << "\n";
    Cmpi::server_out << " Simu::ReadControlParameter() Project::prj_script_dir = " << Project::prj_script_dir << "\n";
    Cmpi::server_out << " Simu::ReadControlParameter() Project::prj_log_dir = " << Project::prj_log_dir << "\n";

    Json::Value root;
    std::ifstream ifs;
    std::string filename = ::add_string( Project::prj_script_dir, "/cfd.json" );
    ifs.open( filename );

    Json::CharReaderBuilder builder;
    builder[ "collectComments" ] = true;
    JSONCPP_STRING errs;
    if ( !Json::parseFromStream(builder, ifs, &root, &errs) )
    {
        Cmpi::server_out << errs << "\n";
        exit( 1 );
    }

    ifs.close();

    if ( root.isObject() )
    {
        std::cout << "root is object " << std::endl;
        Project::simu_task_name = root[ "simutask" ].asString();
        this->Process( root );
    }
    else
    {
        Cmpi::server_out << "root is not object " << "\n";
    }

}

void Simu::Process( Json::Value & root )
{
    TaskLineEnum simu_task = this->GetTaskLine();
    switch ( simu_task )
    {
    case TaskLineEnum::SOLVE_FIELD:
        this->solver_t = new FieldSolver_t();
        this->solver_t->Init( root );
        break;
    case TaskLineEnum::GRID_GEN:
        this->solver_t = new GridSolver_t();
        this->solver_t->Init( root );
        break;
    default:
    {
        std::cerr << "unknown simu_task value!" << std::endl;
        std::exit( EXIT_FAILURE );
    }
    break;
    }
}

void Simu::Run()
{
    this->ReadControlParameter();
    this->solver_t->Run();
}

const std::map<std::string, TaskLineEnum> taskFilter = 
{
    {"undefined",TaskLineEnum::UNDEFINED},
    {"solve_field",TaskLineEnum::SOLVE_FIELD},
    {"grid_gen",TaskLineEnum::GRID_GEN}
};

TaskLineEnum Simu::GetTaskLine()
{
    TaskLineEnum simu_task = TaskLineEnum::UNDEFINED;
    const std::string& taskStr = Project::simu_task_name;
    if ( taskFilter.find(taskStr) != taskFilter.end() )
    {
        simu_task = taskFilter.at(taskStr);
    }
    return simu_task;
}