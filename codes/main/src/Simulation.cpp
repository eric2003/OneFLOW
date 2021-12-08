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
#include "Simulation.h"
#include "SimuImp.h"
#include "SimpleSimu.h"
#include "MpiTest.h"
#include "JsonTest.h"
#include "CgnsTest.h"
#include "HybridParallel.h"
#include <iostream>


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

    //if ( args.size() <= 2 )
    //{
    //    std::cout << " argument number should be 3\n";
    //    exit( 0 );
    //}
}

void Simulation::Run()
{
    int nPara = args.size();
    if ( nPara == 1 )
    {
        this->RunDefaultSimu();
    }
    else if (  nPara == 2 )
    {
        std::cout << " wrong argument number !\n";
        exit( 0 );
    }
    else // nPara >= 3
    { 
        SimuImp * simu = new SimuImp( args );
        simu->Run();
        delete simu;
    }
}

void Simulation::RunDefaultSimu()
{
    //MpiTest * mpiTest = new MpiTest();
    //mpiTest->Run();
    //delete mpiTest;

    //JsonTest * jsonTest = new JsonTest();
    //jsonTest->Run();
    //delete jsonTest;

    //CgnsTest * cgnsTest = new CgnsTest();
    //cgnsTest->Run();
    //delete cgnsTest;

    // CUDA, OpenMP, and MPI parallel
    HybridParallel * hybridParallel = new HybridParallel();
    hybridParallel->Run();
    delete hybridParallel;
}


EndNameSpace
