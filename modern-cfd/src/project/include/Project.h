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
#pragma once
#include <string>

class Project
{
public:
    Project();
    ~Project();
public:
    static std::string system_root;
    static std::string current_dir;
    static std::string execute_dir;
    static std::string prj_rootdir;
    static std::string prj_dir;
    static std::string prj_grid_dir;
    static std::string prj_script_dir;
    static std::string prj_log_dir;
    static std::string prj_results_dir;
public:
    static std::string simu_task_name;
public:
    static void Init( int argc, char **argv );
    static void SetProjectRootDir( const std::string & prjName );
    static void SetProjectDir( const std::string & prjName );
    //static void ReadControlParameter();
public:
    //static TaskLineEnum GetTaskLine();
};

