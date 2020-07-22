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
#include "Simulation.h"
#include "SimuImp.h"
#include "SimpleSimu.h"
#include "MpiTest.h"
#include "JsonTest.h"
#include <iostream>
using namespace std;

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

    cout << "nPara = " << args.size() << "\n";

    for ( int i = 0; i < argc; ++ i )
    {   
        args[ i ] = argv[ i ];
        //cout << "arguments[ " << i << " ] = " << args[ i ] << endl;
        cout << "argv[" << i << "] = " << args[ i ] << endl;
    }

    //if ( args.size() <= 2 )
    //{
    //    cout << " argument number should be 3\n";
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
    else if ( nPara == 2 )
    {
        SimpleSimu * simu = new SimpleSimu( args );
        simu->Run();
        delete simu;
    }
    else
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

    JsonTest * jsonTest = new JsonTest();
    jsonTest->Run();
    delete jsonTest;
}


EndNameSpace
