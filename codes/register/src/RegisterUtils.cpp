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

#include "RegisterUtils.h"
#include "SolverDef.h"
#include "Fatal.h"

BeginNameSpace( ONEFLOW )

VarNameSolver::VarNameSolver()
{
}

VarNameSolver::~VarNameSolver()
{
}

void VarNameSolver::AddFieldName( const std::string & fieldName )
{
    this->data.push_back( fieldName );
}

std::map< int, VarNameSolver * > * VarNameFactory::data = 0;
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
        VarNameFactory::data = new std::map< int, VarNameSolver * >();
        VarNameFactory::mapData = new MapIntInt();
    }
}

void VarNameFactory::AddVarNameSolver( int a, int b )
{
    VarNameFactory::Init();

    VarNameFactory::mapData->AddData( a, b );
    int solverPos = VarNameFactory::mapData->GetId( a, b );

    std::map< int, VarNameSolver * >::iterator iter;

    iter = VarNameFactory::data->find( solverPos );
    if ( iter == VarNameFactory::data->end() )
    {
        VarNameSolver * varNameSolver = new VarNameSolver();
        ( * VarNameFactory::data )[ solverPos ] = varNameSolver;
    }
}

VarNameSolver * VarNameFactory::GetVarNameSolver( int a, int b )
{
    // Required lookup: caller must have registered via AddVarNameSolver.
    // Optional lookup should use FindVarNameSolver instead.
    VarNameSolver * solver =
        VarNameFactory::FindVarNameSolver( a, b );

    if ( solver == nullptr )
    {
        Fatal(
            "VarNameSolver is not registered for the given "
            "solverType / fieldType pair" );
    }

    return solver;
}

void VarNameFactory::FreeVarNameSolver()
{
    if ( ! VarNameFactory::data ) return;
    std::map< int, VarNameSolver * >::iterator iter;
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

VarNameSolver * VarNameFactory::FindVarNameSolver(
    int a,
    int b )
{
    if ( ! VarNameFactory::data ||
        ! VarNameFactory::mapData )
    {
        return nullptr;
    }

    DataAB key;
    key.a = a;
    key.b = b;

    std::map< DataAB, int, CmpDataAB >::iterator iter =
        VarNameFactory::mapData->data.find( key );

    if ( iter == VarNameFactory::mapData->data.end() )
    {
        return nullptr;
    }

    std::map< int, VarNameSolver * >::iterator solverIter =
        VarNameFactory::data->find( iter->second );

    if ( solverIter == VarNameFactory::data->end() )
    {
        return nullptr;
    }

    return solverIter->second;
}

void VarNameFactory::Dump(
    std::ostream & output,
    int solverType )
{
    const int interfaceTypes[] =
    {
        INTERFACE_DATA,
        INTERFACE_DQ_DATA,
        INTERFACE_GRADIENT_DATA,
        INTERFACE_OVERSET_DATA
    };

    const char * interfaceNames[] =
    {
        "INTERFACE_DATA",
        "INTERFACE_DQ",
        "INTERFACE_GRADIENT",
        "INTERFACE_OVERSET"
    };

    output
        << "[Communication Groups]\n";

    const int interfaceTypeCount =
        sizeof( interfaceTypes ) / sizeof( interfaceTypes[ 0 ] );

    for ( int iType = 0; iType < interfaceTypeCount; ++ iType )
    {
        VarNameSolver * varNameSolver =
            VarNameFactory::FindVarNameSolver(
                solverType,
                interfaceTypes[ iType ] );

        output
            << "  "
            << interfaceNames[ iType ]
            << ":\n";

        if ( varNameSolver == nullptr )
        {
            output
                << "    <not registered>\n";
            continue;
        }

        if ( varNameSolver->data.empty() )
        {
            output
                << "    <empty>\n";
            continue;
        }

        for ( int iField = 0;
            iField < varNameSolver->data.size();
            ++ iField )
        {
            output
                << "    "
                << varNameSolver->data[ iField ]
                << '\n';
        }

        if ( interfaceTypes[ iType ] ==
            INTERFACE_OVERSET_DATA )
        {
            output
                << "    Status: Reserved / Not Implemented\n";
        }
    }
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
    std::map< DataAB, int, CmpDataAB >::iterator iter;
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
    DataAB ab;
    ab.a = a;
    ab.b = b;

    std::map< DataAB, int, CmpDataAB >::iterator iter =
        this->data.find( ab );

    if ( iter == this->data.end() )
    {
        Fatal( "MapIntInt key is not registered" );
    }

    return iter->second;
}

EndNameSpace
