#include "StageProfiler.h"

#include <gtest/gtest.h>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

namespace
{
void SetEnv( const char * key, const std::string & value )
{
#ifdef _WIN32
    _putenv_s( key, value.c_str() );
#else
    setenv( key, value.c_str(), 1 );
#endif
}

struct Row
{
    unsigned long long count = 0;
    double inclusive = 0.0;
    double exclusive = 0.0;
};

class StageProfiler : public testing::Test
{
protected:
    std::filesystem::path directory;
    std::string path;

    void SetUp() override
    {
        directory = std::filesystem::temp_directory_path() /
            ( "oneflow-host-accounting-" + std::to_string(
                std::chrono::steady_clock::now().time_since_epoch().count() ) );
        std::filesystem::create_directories( directory );
        path = ( directory / "rank-%r.tsv" ).string();
        SetEnv( "ONEFLOW_STAGE_BREAKDOWN_FILE", path );
        SetEnv( "OMPI_COMM_WORLD_RANK", "3" );
        if ( std::getenv( "ONEFLOW_HOST_ACCOUNTING" ) == nullptr )
            SetEnv( "ONEFLOW_HOST_ACCOUNTING", "1" );
    }

    void TearDown() override
    {
        std::filesystem::remove_all( directory );
    }

    std::filesystem::path HostPath() const
    {
        return directory / "rank-3.tsv.host.tsv";
    }

    std::map< std::string, Row > Read()
    {
        ONEFLOW::StageProfiler::Flush();
        std::ifstream input( HostPath() );
        EXPECT_TRUE( input.good() );
        std::string header;
        std::getline( input, header );
        EXPECT_EQ( header, "scope\tcount\tinclusive_ms\texclusive_ms" );
        std::map< std::string, Row > rows;
        std::string name;
        Row row;
        while ( input >> name >> row.count >> row.inclusive >> row.exclusive )
            rows[ name ] = row;
        return rows;
    }
};

TEST_F( StageProfiler, NestedScopesDoNotDoubleCount )
{
    {
        ONEFLOW::ScopedHostTimer root( "nested.root", true );
        {
            ONEFLOW::ScopedHostTimer parent( "nested.parent" );
            ONEFLOW::ScopedHostTimer leaf( "nested.leaf" );
        }
        ONEFLOW::ScopedHostTimer sibling( "nested.sibling" );
    }
    const auto rows = Read();
    const Row root = rows.at( "nested.root" );
    const Row parent = rows.at( "nested.parent" );
    const Row leaf = rows.at( "nested.leaf" );
    const Row sibling = rows.at( "nested.sibling" );
    EXPECT_EQ( root.count, 1u );
    EXPECT_NEAR( parent.inclusive, parent.exclusive + leaf.inclusive, 1e-9 );
    EXPECT_NEAR( root.inclusive, root.exclusive + parent.exclusive +
                 leaf.exclusive + sibling.exclusive, 1e-9 );
    EXPECT_GE( root.exclusive, 0.0 );
}

TEST_F( StageProfiler, ExceptionUnwindsAndRestoresParent )
{
    {
        ONEFLOW::ScopedHostTimer root( "exception.root", true );
        try
        {
            ONEFLOW::ScopedHostTimer child( "exception.child" );
            throw std::runtime_error( "test" );
        }
        catch ( const std::runtime_error & ) {}
        ONEFLOW::ScopedHostTimer after( "exception.after" );
    }
    const auto rows = Read();
    EXPECT_NEAR( rows.at( "exception.root" ).inclusive,
        rows.at( "exception.root" ).exclusive + rows.at( "exception.child" ).exclusive +
        rows.at( "exception.after" ).exclusive, 1e-9 );
}

TEST_F( StageProfiler, IgnoresWorkOutsideRootAndAccumulatesCalls )
{
    { ONEFLOW::ScopedHostTimer outside( "outside.before" ); }
    for ( int i = 0; i < 2; ++i )
    {
        ONEFLOW::ScopedHostTimer root( "repeated.root", true );
        ONEFLOW::ScopedHostTimer child( "repeated.child" );
    }
    { ONEFLOW::ScopedHostTimer outside( "outside.after" ); }
    const auto rows = Read();
    EXPECT_EQ( rows.count( "outside.before" ), 0u );
    EXPECT_EQ( rows.count( "outside.after" ), 0u );
    EXPECT_EQ( rows.at( "repeated.root" ).count, 2u );
    EXPECT_EQ( rows.at( "repeated.child" ).count, 2u );
    EXPECT_NEAR( rows.at( "repeated.root" ).inclusive,
        rows.at( "repeated.root" ).exclusive + rows.at( "repeated.child" ).exclusive, 1e-9 );
}

TEST_F( StageProfiler, SwitchControlsSidecar )
{
    { ONEFLOW::ScopedHostTimer root( "disabled.root", true ); }
    ONEFLOW::StageProfiler::Flush();
    EXPECT_EQ( std::filesystem::exists( HostPath() ),
        std::string( std::getenv( "ONEFLOW_HOST_ACCOUNTING" ) ) == "1" );
    EXPECT_TRUE( std::filesystem::exists( directory / "rank-3.tsv" ) );
}

TEST_F( StageProfiler, RegistersResidentProfilingCategories )
{
    ONEFLOW::StageProfiler::Record( "residual_zero", 1.25 );
    ONEFLOW::StageProfiler::Record( "host_mirror_refresh", 2.5 );
    ONEFLOW::StageProfiler::Flush();
    std::ifstream input( directory / "rank-3.tsv" );
    std::string line;
    std::getline( input, line );
    std::map< std::string, double > totals;
    while ( std::getline( input, line ) )
    {
        const auto tab = line.find( '\t' );
        ASSERT_NE( tab, std::string::npos );
        totals[ line.substr( 0, tab ) ] = std::stod( line.substr( tab + 1 ) );
    }
    EXPECT_EQ( totals.count( "residual_zero" ), 1u );
    EXPECT_EQ( totals.count( "host_mirror_refresh" ), 1u );
    EXPECT_DOUBLE_EQ( totals[ "residual_zero" ], 1.25 );
    EXPECT_DOUBLE_EQ( totals[ "host_mirror_refresh" ], 2.5 );
}
}
