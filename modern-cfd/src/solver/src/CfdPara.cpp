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
#include "CfdPara.h"
#include "Geom.h"
#include "Project.h"
#include "tools.h"
#include <json/json.h>
#include <fstream>
#include <iostream>

CfdPara::CfdPara()
{
    //this->Init();
}

CfdPara::~CfdPara()
{
    ;
}

void CfdPara::Init( Geom * geom )
{
    this->ReadJsonCfdFile();
    //this->irestart = 0;
    //this->cfl = 0.5;
    //this->simu_time = 0.625;
    //this->cspeed = 1.0;
    this->dt = Geom_t::dx * cfl / cspeed;
    this->fnt = ( simu_time + SMALL ) / dt;
    this->nt = fnt;
}

void CfdPara::ReadJsonCfdFile()
{
    std::cout << " ReadJsonCfdFile() Project::prj_grid_dir = " << Project::prj_grid_dir << "\n";
    std::cout << " ReadJsonCfdFile() Project::prj_script_dir = " << Project::prj_script_dir << "\n";
    std::cout << " ReadJsonCfdFile() Project::prj_log_dir = " << Project::prj_log_dir << "\n";

    Json::Value root;
    std::ifstream ifs;
    std::string filename = ::add_string( Project::prj_script_dir, "/cfd.json" );
    ifs.open( filename );

    Json::CharReaderBuilder builder;
    builder[ "collectComments" ] = true;
    JSONCPP_STRING errs;
    if ( !Json::parseFromStream(builder, ifs, &root, &errs) )
    {
        std::cout << errs << std::endl;
        exit( 1 );
    }

    if ( root.isObject() )
    {
        std::cout << "root is object " << std::endl;
        this->irestart = root[ "irestart" ].asInt();
        this->cspeed = root[ "cspeed" ].asFloat();

        Json::Value item = root[ "timestep" ];
        this->cfl = item[ "cfl" ].asFloat();
        this->simu_time = item[ "simu_time" ].asFloat();

        std::cout << "cfl = " << this->cfl << std::endl;
        std::cout << "simu_time = " << this->simu_time << std::endl;
        std::cout << "cspeed = " << this->cspeed  << std::endl;
    }
    else
    {
        std::cout << "root is not object " << std::endl;
    }

    std::cout << root << std::endl;
    std::cout << "--------------------------------" << std::endl;
    std::string myJsonString = root.toStyledString();
    std::cout << myJsonString << std::endl;

}

