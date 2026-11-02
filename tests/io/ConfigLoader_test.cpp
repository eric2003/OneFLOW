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
#include "ConfigLoader.h"
#include "TextFileParser.h"
#include <gtest/gtest.h>
#include <atomic>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

using namespace ONEFLOW;

namespace {

class TemporaryConfigFile {
public:
    explicit TemporaryConfigFile( const std::string& contents )
    {
        static std::atomic<unsigned int> nextId{ 0 };
        const auto id = nextId.fetch_add( 1 );
        const auto timestamp =
            std::chrono::steady_clock::now().time_since_epoch().count();
        path_ = std::filesystem::temp_directory_path()
            / ( "oneflow_config_loader_test_" + std::to_string( timestamp ) + "_"
                + std::to_string( id ) + ".txt" );
        std::ofstream output( path_ );
        if ( ! output )
        {
            throw std::runtime_error( "unable to create temporary configuration file" );
        }
        output << contents;
        if ( ! output )
        {
            throw std::runtime_error( "unable to write temporary configuration file" );
        }
    }

    ~TemporaryConfigFile()
    {
        std::error_code error;
        std::filesystem::remove( path_, error );
    }

    std::string Path() const { return path_.string(); }

private:
    std::filesystem::path path_;
};

} // namespace

TEST( ConfigLoaderTest, PreservesScalarTypesValuesAndDeclarationOrder )
{
    TemporaryConfigFile file(
        "# full line comment\n"
        "int maxIter = 42;\n"
        "Real cfl = 0.5; // inline comment\n"
        "string simutask = \"Solve\";\n" );

    ConfigLoader loader;
    loader.ParseFile( file.Path() );

    const auto& entries = loader.Document().Entries();
    ASSERT_EQ( entries.size(), 3u );
    EXPECT_EQ( entries[ 0 ].typeName, "int" );
    EXPECT_EQ( entries[ 0 ].name, "maxIter" );
    ASSERT_EQ( entries[ 0 ].values.size(), 1u );
    EXPECT_EQ( entries[ 0 ].values[ 0 ], "42" );
    EXPECT_EQ( entries[ 1 ].typeName, "Real" );
    EXPECT_EQ( entries[ 1 ].name, "cfl" );
    EXPECT_EQ( entries[ 1 ].values[ 0 ], "0.5" );
    EXPECT_EQ( entries[ 2 ].typeName, "string" );
    EXPECT_EQ( entries[ 2 ].name, "simutask" );
    EXPECT_EQ( entries[ 2 ].values[ 0 ], "Solve" );
}

TEST( ConfigLoaderTest, ParsesLiteralArrayAcrossLines )
{
    TemporaryConfigFile file(
        "Real gridTrans[3] = 1.0, 2.0,\n"
        "    3.0;\n" );

    ConfigLoader loader;
    loader.ParseFile( file.Path() );

    const auto& entries = loader.Document().Entries();
    ASSERT_EQ( entries.size(), 1u );
    EXPECT_EQ( entries[ 0 ].typeName, "Real" );
    EXPECT_EQ( entries[ 0 ].name, "gridTrans" );
    ASSERT_EQ( entries[ 0 ].values.size(), 3u );
    EXPECT_EQ( entries[ 0 ].values[ 0 ], "1.0" );
    EXPECT_EQ( entries[ 0 ].values[ 1 ], "2.0" );
    EXPECT_EQ( entries[ 0 ].values[ 2 ], "3.0" );
}

TEST( ConfigLoaderTest, ResolvesNamedArrayLengthThroughInjectedCallback )
{
    TemporaryConfigFile file( "int cellOwners[nCells] = 4, 5;\n" );
    std::string requestedLengthName;
    ConfigLoader loader( [&]( const std::string& name ) {
        requestedLengthName = name;
        return 2;
    } );

    loader.ParseFile( file.Path() );

    EXPECT_EQ( requestedLengthName, "nCells" );
    const auto& entries = loader.Document().Entries();
    ASSERT_EQ( entries.size(), 1u );
    EXPECT_EQ( entries[ 0 ].name, "cellOwners" );
    ASSERT_EQ( entries[ 0 ].values.size(), 2u );
    EXPECT_EQ( entries[ 0 ].values[ 0 ], "4" );
    EXPECT_EQ( entries[ 0 ].values[ 1 ], "5" );
}

TEST( ConfigLoaderTest, RejectsNamedArrayLengthWithoutResolver )
{
    TemporaryConfigFile file( "int cellOwners[nCells] = 4, 5;\n" );
    ConfigLoader loader;

    EXPECT_THROW( loader.ParseFile( file.Path() ), std::runtime_error );
}

TEST( ConfigLoaderTest, ParsingAnotherFileReplacesPreviousDocument )
{
    TemporaryConfigFile firstFile( "int first = 1;\n" );
    TemporaryConfigFile secondFile( "int second = 2;\n" );
    ConfigLoader loader;

    loader.ParseFile( firstFile.Path() );
    loader.ParseFile( secondFile.Path() );

    const auto& entries = loader.Document().Entries();
    ASSERT_EQ( entries.size(), 1u );
    EXPECT_EQ( entries[ 0 ].name, "second" );
}
