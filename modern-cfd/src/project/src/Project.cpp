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
#include "Project.h"
#include <filesystem>
#include <iostream>

template <class... Ts>
void print_all(std::ostream& os, Ts const&... args) {
    ((os << args << "\n" ), ... );
}

template <class... Ts>
std::string add_string( Ts const&... args )
{
    std::ostringstream oss;
    ((oss << args), ... );
    return oss.str();
}

std::string Project::system_root;
std::string Project::current_dir;
std::string Project::execute_dir;
std::string Project::prj_rootdir;


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
    std::cout << "std::filesystem::absolute(exeName) = " << std::filesystem::absolute(exeName) << std::endl;
    return std::filesystem::absolute( exeName ).parent_path().string();
}


void Project::Init(int argc, char **argv)
{
    std::cout << "argv[0]=" << argv[0] << "\n";
    Project::execute_dir = get_exe_fullpath( argv[0] );
    Project::current_dir = std::filesystem::current_path().string();

    std::cout << " Project::execute_dir = " << Project::execute_dir << "\n";
    std::cout << " Project::current_dir = " << Project::current_dir << "\n";

    //Project::SetProjectRootDir( prjName );

    //print_all(std::cout, 1, 2, 3);
    //print_all(std::cout, 1, "Project::execute_dir", 3);

    //std::cout << add_string( 1, 2, 3 ) << std::endl;
    //std::cout << add_string( 1, "2system ", 3 ) << std::endl;
    //std::cout << "no value:" << add_string() << std::endl;

    Project::prj_rootdir = Project::current_dir;
}

void Project::SetProjectRootDir( const std::string & prjName )
{
    Project::prj_rootdir = add_string( Project::current_dir, "/", prjName );
}