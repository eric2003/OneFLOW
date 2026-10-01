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

#include "SolverInfo.h"
#include <map>
#include <memory>
#include <iostream>


BeginNameSpace( ONEFLOW )

SolverInfo::SolverInfo()
{
    ;
}

SolverInfo::~SolverInfo()
{
    ;
}

std::unique_ptr< std::map< int, std::unique_ptr< SolverInfo > > > SolverInfoFactory::data;

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
        SolverInfoFactory::data = std::make_unique< std::map< int, std::unique_ptr< SolverInfo > > >();
    }
}

void SolverInfoFactory::AddSolverInfo( int solverType )
{
    SolverInfoFactory::Init();

    auto iter = SolverInfoFactory::data->find( solverType );
    if ( iter == SolverInfoFactory::data->end() )
    {
        ( * SolverInfoFactory::data )[ solverType ] = std::make_unique< SolverInfo >();
    }
}

SolverInfo * SolverInfoFactory::GetSolverInfo( int solverType )
{
    auto iter = SolverInfoFactory::data->find( solverType );
    return iter->second.get();
}

void SolverInfoFactory::Free()
{
    if ( ! SolverInfoFactory::data ) return;
    SolverInfoFactory::data->clear();
    SolverInfoFactory::data.reset();
}


EndNameSpace
