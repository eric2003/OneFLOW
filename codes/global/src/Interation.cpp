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
#include "Iteration.h"
#include "DataBase.h"
#include "Ctrl.h"

BeginNameSpace( ONEFLOW )

int Iteration::innerSteps = 0;
int Iteration::outerSteps = 0;
int Iteration::maxSteps = 1000;
int Iteration::dualtime = 0;
int Iteration::nFieldSave = 100;
int Iteration::nVisualSave = 100;
int Iteration::nResSave = 1;
int Iteration::nForceSave = 1;

Real Iteration::cfl = 1.0;
Real Iteration::cflst = 1.0;
Real Iteration::cfled = 1.0;
int  Iteration::ncfl = 100;

Iteration::Iteration()
{
    ;
}

Iteration::~Iteration()
{
    ;
}

void Iteration::Init()
{
    Iteration::innerSteps = 0;
    Iteration::outerSteps = 0;

    Iteration::cflst = GetDataValue< Real >( "cflst" );
    Iteration::cfled = GetDataValue< Real >( "cfled" );
    Iteration::ncfl  = GetDataValue< int >( "ncfl" );

    Iteration::maxSteps    = GetDataValue< int >( "maxSteps" );
    Iteration::nForceSave  = GetDataValue< int >( "nForceSave" );
    Iteration::nFieldSave  = GetDataValue< int >( "nFieldSave" );
    Iteration::nVisualSave = GetDataValue< int >( "nVisualSave" );
    Iteration::nResSave    = GetDataValue< int >( "nResSave" );
}

bool Iteration::InnerOk()
{
    return Iteration::innerSteps == 1;
}

bool Iteration::ResOk()
{
	int startStrategy = ONEFLOW::GetDataValue< int >("startStrategy");
	if (startStrategy == 2)
	{
		if (Iteration::innerSteps % nResSave == 0)
		{
			return true;
		}
		return false;
	}
	else
	{
		if (Iteration::outerSteps % nResSave == 0)
		{
			return Iteration::innerSteps == 1;
		}
		return false;
	}
}

bool Iteration::ForceOk()
{
	int startStrategy = ONEFLOW::GetDataValue< int >("startStrategy");
	if (startStrategy == 2)
	{
		if (Iteration::innerSteps % nForceSave == 0)
		{
			return true;
		}
		return false;
	}
	else
	{
		if (Iteration::outerSteps % nForceSave == 0)
		{
			return Iteration::innerSteps == 1;
		}
		return false;
	}
}

SimuIterState::SimuIterState()
{
    ;
}

SimuIterState::~SimuIterState()
{
    ;
}

bool SimuIterState::Running()
{
    if ( ctrl.iexitflag == 0 )
    {
        if ( Iteration::outerSteps >= Iteration::maxSteps )
        {
            return false;
        }
        return true;
    }
    else if ( ctrl.currTime >= ctrl.maxTime )
    {
        return false;
    }
    else
    {
        return true;
    }
}

bool SimuIterState::InnerRunning()
{
    return true;
}


int SaveState::nFileSaveSteps = 1;

SaveState::SaveState()
{
    ;
}

SaveState::~SaveState()
{
    ;
}


EndNameSpace