/*---------------------------------------------------------------------------*\
OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
Copyright (C) 2017-2023 He Xin and the OneFLOW contributors.
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
#include "Project.h"
#include "tools.h"
#include "Cmpi.h"
#include <json/json.h>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <map>

std::string Project::system_root;
std::string Project::current_dir;
std::string Project::execute_dir;
std::string Project::prj_rootdir;
std::string Project::prj_dir;
std::string Project::prj_grid_dir;
std::string Project::prj_script_dir;
std::string Project::prj_log_dir;
std::string Project::prj_results_dir;
std::string Project::simu_task_name;

Project::Project()
{
}

Project::~Project()
{
}

std::string get_exe_fullpath( const std::string &exeName )
{
    std::filesystem::path full_path;
    full_path = std::filesystem::current_path();
    Cmpi::server_out << "std::filesystem::absolute(exeName) = " << std::filesystem::absolute( exeName ) << "\n";
    return std::filesystem::absolute( exeName ).parent_path().string();
}


void Project::Init(int argc, char **argv)
{
    Cmpi::server_out << "argv[0]=" << argv[ 0 ] << "\n";

    Project::execute_dir = get_exe_fullpath( argv[0] );
    Project::current_dir = std::filesystem::current_path().string();

    Cmpi::server_out << " Project::execute_dir = " << Project::execute_dir << "\n";
    Cmpi::server_out << " Project::current_dir = " << Project::current_dir << "\n";

    Project::prj_rootdir = Project::current_dir;
    Project::prj_grid_dir = "grid";
    Project::prj_script_dir = "script";
    Project::prj_log_dir = "log";

    //Project::SetProjectRootDir("workdir/oneflow1d");
    Cmpi::server_out << " Project::prj_rootdir = " << Project::prj_rootdir << "\n";
    Project::prj_rootdir = std::filesystem::current_path().parent_path().parent_path().string() + "/project_tmp";
    Cmpi::server_out << " Project::prj_rootdir = " << Project::prj_rootdir << "\n";
    std::string prjName = "myprj";
    if ( argc > 1 )
    {
        Cmpi::server_out << "argv[1] = " << argv[1] << "\n";
        prjName = argv[ 1 ];
    }
    Project::SetProjectDir( prjName );
    Cmpi::server_out << " Project::prj_dir = " << Project::prj_dir << "\n";
}

void Project::SetProjectDir( const std::string & prjName )
{
    Project::prj_dir = add_string( Project::prj_rootdir, "/", prjName );
    Project::prj_grid_dir = add_string( Project::prj_dir, "/", "grid" );
    Project::prj_script_dir = add_string( Project::prj_dir, "/", "script" );
    Project::prj_log_dir = add_string( Project::prj_dir, "/", "log" );
    Project::prj_results_dir = add_string( Project::prj_dir, "/", "results" );
    Cmpi::server_out << " Project::prj_grid_dir = " << Project::prj_grid_dir << "\n";
    Cmpi::server_out << " Project::prj_script_dir = " << Project::prj_script_dir << "\n";
    Cmpi::server_out << " Project::prj_log_dir = " << Project::prj_log_dir << "\n";
}

void Project::SetProjectRootDir( const std::string & prjName )
{
    std::filesystem::path myprjpath = prjName;
    Cmpi::server_out << " prjName = " << prjName << "\n";

    if ( myprjpath.is_relative() )
    {
        std::cout << " prjName is relative path\n";
        Project::prj_rootdir = add_string( Project::current_dir, "/", prjName );
    }
}
