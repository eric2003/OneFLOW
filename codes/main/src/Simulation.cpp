/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
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
#include "Simulation.h"
#include "SimuImp.h"
#include "SimpleSimu.h"
#include "MpiTest.h"
#include "JsonTest.h"
#include "CgnsTest.h"
#include "HybridParallel.h"
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>


BeginNameSpace( ONEFLOW )

Simulation::Simulation( int argc, char ** argv )
{
    this->ProcessCmdLineArgs( argc, argv );
}

Simulation::~Simulation()
{
    if ( ! args.empty() ) 
    {
        args.clear();
        args.shrink_to_fit();
    };
}

void Simulation::ProcessCmdLineArgs( int argc, char ** argv )
{
    args.resize( argc );

    std::cout << "nPara = " << args.size() << "\n";

    for ( int i = 0; i < argc; ++ i )
    {   
        args[ i ] = argv[ i ];
        //std::cout << "arguments[ " << i << " ] = " << args[ i ] << std::endl;
        std::cout << "argv[" << i << "] = " << args[ i ] << std::endl;
    }
}

std::unique_ptr<SimuBase> Simulation::MakeDefaultSimulation()
{
    const char* envName = std::getenv("ONEFLOW_DEFAULT_TEST");
    auto& reg = TestRegistry::Instance();

    auto nameList = reg.GetAllRegisteredNames();
    std::sort(nameList.begin(), nameList.end());

    std::cerr << "\n===== OneFLOW lightweight test selector =====\n";
    std::cerr << "Available test case values:\n";
    for (const auto& name : nameList)
    {
        std::cerr << "    - " << name << "\n";
    }

    const std::string defaultCase = "hybrid_parallel";
    std::string selectedCase;
    std::unique_ptr<SimuBase> ptr;

    if (envName != nullptr && std::strlen(envName) > 0)
    {
        selectedCase = envName;
        ptr = reg.Create(selectedCase);
        if (!ptr)
        {
            std::cerr << "\n[WARNING] Test case \"" << selectedCase << "\" is not registered.\n";
            std::cerr << "Fallback to built-in default test case.\n";
            selectedCase = defaultCase;
            ptr = reg.Create(selectedCase);
        }
    }
    else
    {
        // environment variable is not set, use built-in default
        std::cerr << "\n[INFO] Environment variable ONEFLOW_DEFAULT_TEST is NOT set.\n";
        selectedCase = defaultCase;
        ptr = reg.Create(selectedCase);

        std::cerr << "\nHow to set ONEFLOW_DEFAULT_TEST:\n";
        std::cerr << "  Linux / macOS (bash/zsh):\n";
        std::cerr << "      export ONEFLOW_DEFAULT_TEST=\"mpi_test\"\n";
        std::cerr << "  Windows Command Prompt (cmd):\n";
        std::cerr << "      set ONEFLOW_DEFAULT_TEST=mpi_test\n";
        std::cerr << "  Windows PowerShell:\n";
        std::cerr << "      $env:ONEFLOW_DEFAULT_TEST=\"mpi_test\"\n";
    }

    std::cerr << "\n[INFO] Currently using test case: \"" << selectedCase << "\"\n";
    std::cerr << "============================================\n\n";

    return ptr;
}

void Simulation::RunImpl()
{
    int nPara = static_cast<int>( args.size() );

    if ( nPara == 1 )
    {
        std::cout << "\n===== ONEFLOW Light-weight Test Mode =====\n";
        auto simu = MakeDefaultSimulation();
        if ( ! simu )
        {
            // Prefer throwing so it is handled uniformly
            throw std::runtime_error( "[Error] No available default test case!" );
        }
        simu->Run();
    }
    else if ( nPara == 2 )
    {
        // Unified error handling
        throw std::runtime_error( "wrong argument number !" );
    }
    else // nPara >= 3
    {
        std::cout << "\n===== ONEFLOW Full Simulation Mode =====\n";
        auto simu = std::make_unique<SimuImp>( args );
        simu->Run();
    }
}


int Simulation::Run()
{
    try
    {
        RunImpl();
        return 0;
    }
    catch ( const std::exception & e )
    {
        std::cerr << "\n========== Fatal Error ==========\n"
            << e.what() << "\n"
            << "=================================\n";
        return EXIT_FAILURE;
    }
    catch ( ... )
    {
        std::cerr << "\n========== Unknown Fatal Error ==========\n";
        return EXIT_FAILURE;
    }
}

EndNameSpace
