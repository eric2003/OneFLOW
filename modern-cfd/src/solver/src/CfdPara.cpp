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
#include "Cmpi.h"
#include "Geom.h"
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

void CfdPara::Init( Json::Value & root )
{
    this->irestart = root[ "irestart" ].asInt();
    this->cspeed = root[ "cspeed" ].asFloat();

    Json::Value item = root[ "timestep" ];
    this->cfl = item[ "cfl" ].asFloat();
    this->simu_time = item[ "simu_time" ].asFloat();

    item = root[ "solver" ];

    std::string gridobj_str = item[ "gridobj" ].asString();
    Geom_t::gridName = item[ "gridfile" ].asString();
    if ( gridobj_str == "read" )
    {
        Geom_t::gridobj = 1;
    }
    else
    {
        Geom_t::gridobj = 0;
    }

    Cmpi::server_out << "cfl = " << this->cfl << "\n";
    Cmpi::server_out << "simu_time = " << this->simu_time << "\n";
    Cmpi::server_out << "cspeed = " << this->cspeed << "\n";
}


