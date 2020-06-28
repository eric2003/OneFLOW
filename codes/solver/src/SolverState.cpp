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

#include "SolverState.h"
#include "SolverMap.h"
#include "CmxTask.h"
#include "Iteration.h"
#include "Zone.h"
#include "ZoneState.h"
#include "GridState.h"
#include "Ctrl.h"

BeginNameSpace( ONEFLOW )

HXVector< LusgsSolver * > LusgsState::str;
HXVector< LusgsSolver * > LusgsState::uns;

LusgsState::LusgsState()
{
    ;
}

LusgsState::~LusgsState()
{
    ;
}

void LusgsState::Init( int nSolver )
{
    LusgsState::str.resize( nSolver );
    LusgsState::uns.resize( nSolver );
}

void LusgsState::AddSolver( int sid, int gridType, LusgsSolver * solver )
{
    if ( gridType == ONEFLOW::UMESH )
    {
        if ( LusgsState::uns[ sid ] ) return;
        LusgsState::uns[ sid ] = solver;
    }
    else
    {
        if ( LusgsState::str[ sid ] ) return;
        LusgsState::str[ sid ] = solver;
    }
}

LusgsSolver * LusgsState::GetLusgsSolver()
{
    int gridType = ZoneState::zoneType[ ZoneState::zid ];
    if ( gridType == ONEFLOW::UMESH )
    {
        return LusgsState::uns[ SolverState::id ];
    }
    else
    {
        return LusgsState::str[ SolverState::id ];
    }
}

int SolverState::id = 0;
int SolverState::tid = 0;
int SolverState::nSolver = 0;
int SolverState::msgId = -1;
IntField SolverState::convergeFlag;

SolverState::SolverState()
{
    ;
}

SolverState::~SolverState()
{
    ;
}

void SolverState::Init( int nSolver )
{
    SolverState::nSolver = nSolver;
    SolverState::convergeFlag.resize( nSolver );
}

void SolverState::SetTid( int tid )
{
    SolverState::tid = tid;
}

void SolverState::SetTidById( int id )
{
    SolverState::id  = id;
    SolverState::tid = SolverMap::GetTid( id );
}

Solver * SolverState::GetSolver()
{
    int gridType = ZoneState::zoneType[ ZoneState::zid ];
    return SolverMap::GetSolver( SolverState::id, gridType );
}

bool SolverState::Converge()
{
    if ( Iteration::innerSteps == 0 ) return false;
    if ( ctrl.idualtime == 0 ) return true;
    bool flag = true;
    for ( int iSolver = 0; iSolver < SolverState::nSolver; ++ iSolver )
    {
        SolverState::id = iSolver;
        ONEFLOW::SsSgTask( "CALC_UNSTEADY_CRITERION" );
    }
    
    return flag;
}


EndNameSpace