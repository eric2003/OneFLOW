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

#include "RegisterUtil.h"

BeginNameSpace( ONEFLOW )

VarNameSolver::VarNameSolver()
{
}

VarNameSolver::~VarNameSolver()
{
}

void VarNameSolver::AddFieldName( const string & fieldName )
{
    this->data.push_back( fieldName );
}

map< int, VarNameSolver * > * VarNameFactory::data = 0;
MapIntInt * VarNameFactory::mapData = 0;

VarNameFactory::VarNameFactory()
{
}

VarNameFactory::~VarNameFactory()
{
}

void VarNameFactory::Init()
{
    if ( ! VarNameFactory::data )
    {
        VarNameFactory::data = new map< int, VarNameSolver * >();
        VarNameFactory::mapData = new MapIntInt();
    }
}

void VarNameFactory::AddVarNameSolver( int a, int b )
{
    VarNameFactory::Init();

    VarNameFactory::mapData->AddData( a, b );
    int solverPos = VarNameFactory::mapData->GetId( a, b );

    map< int, VarNameSolver * >::iterator iter;

    iter = VarNameFactory::data->find( solverPos );
    if ( iter == VarNameFactory::data->end() )
    {
        VarNameSolver * varNameSolver = new VarNameSolver();
        ( * VarNameFactory::data )[ solverPos ] = varNameSolver;
    }
}

VarNameSolver * VarNameFactory::GetVarNameSolver( int a, int b )
{
    VarNameFactory::Init();

    int solverId = VarNameFactory::mapData->GetId( a, b );

    map< int, VarNameSolver * >::iterator iter;
    iter = VarNameFactory::data->find( solverId );
    return iter->second;
}

void VarNameFactory::FreeVarNameSolver()
{
    if ( ! VarNameFactory::data ) return;
    map< int, VarNameSolver * >::iterator iter;
    for ( iter = VarNameFactory::data->begin(); iter != VarNameFactory::data->end(); ++ iter )
    {
        delete iter->second;
    }

    VarNameFactory::data->clear();

    delete VarNameFactory::data;
    VarNameFactory::data = 0;

    delete VarNameFactory::mapData;
    VarNameFactory::mapData = 0;
}

bool CmpDataAB::operator()( const DataAB & k1, const DataAB & k2 ) const
{
    if ( k1.a != k2.a )
    {
        return k1.a < k2.a;
    }

    return k1.b < k2.b;
}


MapIntInt::MapIntInt()
{
}

MapIntInt::~MapIntInt()
{
}

void MapIntInt::AddData( int a, int b )
{
    map< DataAB, int, CmpDataAB >::iterator iter;
    DataAB ab;
    ab.a = a;
    ab.b = b;
    iter = this->data.find( ab );
    int index = this->data.size();
    if ( iter == this->data.end() )
    {
        this->data[ ab ] = index;
    }
}

int MapIntInt::GetId( int a, int b )
{
    map< DataAB, int, CmpDataAB >::iterator iter;
    DataAB ab;
    ab.a = a;
    ab.b = b;
    iter = this->data.find( ab );
    return iter->second;
}

EndNameSpace