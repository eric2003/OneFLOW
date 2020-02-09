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

#include "SolverInfo.h"
#include <map>
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

SolverInfo::SolverInfo()
{
    ;
}

SolverInfo::~SolverInfo()
{
    ;
}

map< int, SolverInfo * > * SolverInfoFactory::data = 0;

SolverInfoFactory::SolverInfoFactory()
{
}

SolverInfoFactory::~SolverInfoFactory()
{
}

void SolverInfoFactory::Init()
{
    if ( ! SolverInfoFactory::data )
    {
        SolverInfoFactory::data = new map< int, SolverInfo * >();
    }
}

void SolverInfoFactory::AddSolverInfo( int sTid )
{
    SolverInfoFactory::Init();

    map< int, SolverInfo * >::iterator iter;

    iter = SolverInfoFactory::data->find( sTid );
    if ( iter == SolverInfoFactory::data->end() )
    {
        SolverInfo * solverInfo = new SolverInfo();
        ( * SolverInfoFactory::data )[ sTid ] = solverInfo;
    }
}

SolverInfo * SolverInfoFactory::GetSolverInfo( int sTid )
{
    map< int, SolverInfo * >::iterator iter;
    iter = SolverInfoFactory::data->find( sTid );
    return iter->second;
}

void SolverInfoFactory::Free()
{
    if ( ! SolverInfoFactory::data ) return;
    map< int, SolverInfo * >::iterator iter;
    for ( iter = SolverInfoFactory::data->begin(); iter != SolverInfoFactory::data->end(); ++ iter )
    {
        delete iter->second;
    }

    SolverInfoFactory::data->clear();

    delete SolverInfoFactory::data;
    SolverInfoFactory::data = 0;
}


EndNameSpace