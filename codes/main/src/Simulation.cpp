/*---------------------------------------------------------------------------*\
OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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
#include "Prj.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

Simulation::Simulation( int argc, char ** argv )
{
	this->argc = argc;
	this->argv = argv;
	cout << " argc = " << argc << "\n";

    for ( int i = 0; i < argc; ++ i )
    {    
        cout<<"argument["<<i<<"] is: " << argv[ i ] << endl;
    }

	if ( this->argc > 1 )
	{
		string prjName = argv[ 1 ];
		PrjStatus::SetPrjBaseDir( prjName );
	}
}

Simulation::~Simulation()
{
	;
}

void Simulation::Run()
{
	this->PreProcess();
	this->MainProcess();
	this->PostProcess();
}

void Simulation::PreProcess()
{
    InitSimu();
}

void Simulation::MainProcess()
{
    RunSimu();
}

void Simulation::PostProcess()
{
}


EndNameSpace
