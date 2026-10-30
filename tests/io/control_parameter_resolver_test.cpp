/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "ParaFile.h"
#include "ConfigDatabaseAdapter.h"
#include "ConfigDocument.h"
#include "DataBase.h"
#include <gtest/gtest.h>
#include <stdexcept>
#include <string>

using namespace ONEFLOW;

TEST( ControlParameterResolverTest, ParsesLiteralArrayLength )
{
    EXPECT_EQ( GetParameterArraySize( "0" ), 0 );
    EXPECT_EQ( GetParameterArraySize( "17" ), 17 );
}

TEST( ControlParameterResolverTest, ResolvesNamedLengthFromRuntimeDatabase )
{
    SetDataInt( "nCells", 4 );

    EXPECT_EQ( GetParameterArraySize( "nCells" ), 4 );
}

TEST( ControlParameterResolverTest, MissingNamedLengthPreservesDatabaseError )
{
    EXPECT_THROW(
        GetParameterArraySize( "missingArrayLength" ),
        std::runtime_error );
}

TEST( ConfigDatabaseAdapterTest, CommitsTypedEntriesIntoRuntimeDatabase )
{
    ConfigDocument document;
    document.Add( ParameterEntry{ "dimension", "int", { "2" } } );
    document.Add( ParameterEntry{ "global_dt", "Real", { "0.001" } } );
    document.Add( ParameterEntry{ "simutask", "string", { "Grid" } } );

    ConfigDatabaseAdapter::Commit( document );

    EXPECT_EQ( GetDataValue<int>( "dimension" ), 2 );
    EXPECT_DOUBLE_EQ( GetDataValue<Real>( "global_dt" ), 0.001 );
    EXPECT_EQ( GetDataValue<std::string>( "simutask" ), "Grid" );
}

TEST( ConfigDatabaseAdapterTest, LaterDocumentEntryReplacesEarlierValue )
{
    ConfigDocument document;
    document.Add( ParameterEntry{ "dimension", "int", { "2" } } );
    document.Add( ParameterEntry{ "dimension", "int", { "3" } } );

    ConfigDatabaseAdapter::Commit( document );

    EXPECT_EQ( GetDataValue<int>( "dimension" ), 3 );
}
