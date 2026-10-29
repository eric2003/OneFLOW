#include <gtest/gtest.h>

#include "DataBaseType.h"
#include "DataObject.h"
#include "DataPara.h"
#include "SimuContextParameter.h"
#include "SimuContextIO.h"
#include "Prj.h"

#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace ONEFLOW;

namespace
{

void AddStringParameter(
    DataPara & parameters, const std::string & name, const std::string & value )
{
    auto dataObject = std::make_unique< TDataObject< std::string > >( 1 );
    dataObject->CopyValue( &value, 1 );
    parameters.SetDataEntry( std::make_unique< DataEntry >(
        name, HX_STRING, 1, std::move( dataObject ) ) );
}

void AddIntegerParameter(
    DataPara & parameters, const std::string & name, int value )
{
    auto dataObject = std::make_unique< TDataObject< int > >( 1 );
    dataObject->CopyValue( &value, 1 );
    parameters.SetDataEntry( std::make_unique< DataEntry >(
        name, HX_INT, 1, std::move( dataObject ) ) );
}

std::string ReadTextFile( const std::filesystem::path & path )
{
    std::ifstream file( path );
    return std::string(
        std::istreambuf_iterator<char>( file ),
        std::istreambuf_iterator<char>() );
}

}

TEST( SimuContextParameterTest, MissingOptionalStringReturnsEmpty )
{
    DataPara parameters;

    EXPECT_EQ(
        SimuContextDetail::ReadOptionalStringParameter(
            parameters, "inputRoot" ),
        "" );
}

TEST( SimuContextParameterTest, ReadsStringParameter )
{
    DataPara parameters;
    AddStringParameter( parameters, "inputRoot", "mesh/input" );

    EXPECT_EQ(
        SimuContextDetail::ReadOptionalStringParameter(
            parameters, "inputRoot" ),
        "mesh/input" );
}

TEST( SimuContextParameterTest, RejectsNonStringAndNamesParameter )
{
    DataPara parameters;
    AddIntegerParameter( parameters, "inputRoot", 12 );

    try
    {
        SimuContextDetail::ReadOptionalStringParameter(
            parameters, "inputRoot" );
        FAIL() << "Expected a type validation exception";
    }
    catch ( const std::invalid_argument & error )
    {
        EXPECT_NE(
            std::string( error.what() ).find( "inputRoot" ),
            std::string::npos );
        EXPECT_NE(
            std::string( error.what() ).find( "must be a string" ),
            std::string::npos );
    }
}

TEST( SimuContextParameterTest, RejectsNonStringForEveryCaseIOSetting )
{
    const std::vector< std::string > parameterNames = {
        "projectRoot",
        "inputRoot",
        "outputRoot",
        "meshInputRoot",
        "flowInputRoot",
        "meshOutputRoot",
        "flowOutputRoot",
        "ioPathInfoFileName",
        "sourceGridCaseDir"
    };

    for ( const std::string & name : parameterNames )
    {
        SCOPED_TRACE( name );
        DataPara parameters;
        AddIntegerParameter( parameters, name, 7 );

        try
        {
            SimuContextDetail::ReadCaseIOSettings( parameters );
            FAIL() << "Expected a type validation exception";
        }
        catch ( const std::invalid_argument & error )
        {
            EXPECT_NE(
                std::string( error.what() ).find( name ),
                std::string::npos );
            EXPECT_NE(
                std::string( error.what() ).find( "must be a string" ),
                std::string::npos );
        }
    }
}

TEST( SimuContextParameterTest, EmptyMeshInputAllowsLegacyFallback )
{
    EXPECT_TRUE(
        SimuContextDetail::ShouldUseLegacySourceGridCaseDir( "" ) );
}

TEST( SimuContextParameterTest, DefaultMeshInputAllowsLegacyFallback )
{
    EXPECT_TRUE(
        SimuContextDetail::ShouldUseLegacySourceGridCaseDir( "default" ) );
}

TEST( SimuContextParameterTest, ExplicitMeshInputTakesPrecedence )
{
    EXPECT_FALSE(
        SimuContextDetail::ShouldUseLegacySourceGridCaseDir( "grid-input" ) );
}

TEST( SimuContextParameterTest, CaseRootDoesNotUseLegacyFallback )
{
    EXPECT_FALSE(
        SimuContextDetail::ShouldUseLegacySourceGridCaseDir( "case" ) );
}

TEST( SimuContextParameterTest, RelativeMeshRootDoesNotUseLegacyFallback )
{
    EXPECT_FALSE(
        SimuContextDetail::ShouldUseLegacySourceGridCaseDir( "mesh/input" ) );
}

TEST( SimuContextParameterTest, LegacySourceGridReplacesEmptyMeshRoot )
{
    EXPECT_EQ(
        SimuContextDetail::ApplyLegacySourceGridCaseDir( "", "../grid-case" ),
        "../grid-case" );
}

TEST( SimuContextParameterTest, LegacySourceGridReplacesDefaultMeshRoot )
{
    EXPECT_EQ(
        SimuContextDetail::ApplyLegacySourceGridCaseDir(
            "default", "../grid-case" ),
        "../grid-case" );
}

TEST( SimuContextParameterTest, ReadsGenericAndDomainRootSettings )
{
    DataPara parameters;
    AddStringParameter( parameters, "projectRoot", "../project" );
    AddStringParameter( parameters, "inputRoot", "shared-input" );
    AddStringParameter( parameters, "outputRoot", "shared-output" );
    AddStringParameter( parameters, "meshInputRoot", "../mesh-source" );
    AddStringParameter( parameters, "flowInputRoot", "../flow-source" );
    AddStringParameter( parameters, "meshOutputRoot", "mesh-results" );
    AddStringParameter( parameters, "flowOutputRoot", "flow-results" );
    AddStringParameter( parameters, "ioPathInfoFileName", "diagnostics/paths.txt" );
    AddStringParameter( parameters, "sourceGridCaseDir", "../legacy-grid" );

    const auto settings = SimuContextDetail::ReadCaseIOSettings( parameters );

    EXPECT_EQ( settings.projectRoot, "../project" );
    EXPECT_EQ( settings.inputRoot, "shared-input" );
    EXPECT_EQ( settings.outputRoot, "shared-output" );
    EXPECT_EQ( settings.meshInputRoot, "../mesh-source" );
    EXPECT_EQ( settings.flowInputRoot, "../flow-source" );
    EXPECT_EQ( settings.meshOutputRoot, "mesh-results" );
    EXPECT_EQ( settings.flowOutputRoot, "flow-results" );
    EXPECT_EQ( settings.ioPathInfoFileName, "diagnostics/paths.txt" );
}

TEST( SimuContextParameterTest, ReadsLegacyGridRootOnlyWhenMeshRootAllowsFallback )
{
    DataPara parameters;
    AddStringParameter( parameters, "meshInputRoot", "default" );
    AddStringParameter( parameters, "sourceGridCaseDir", "../producer" );

    const auto settings = SimuContextDetail::ReadCaseIOSettings( parameters );

    EXPECT_EQ( settings.meshInputRoot, "../producer" );
}

TEST( SimuContextParameterTest, ExplicitMeshRootDoesNotReadLegacyGridParameter )
{
    DataPara parameters;
    AddStringParameter( parameters, "meshInputRoot", "../explicit-grid" );
    AddIntegerParameter( parameters, "sourceGridCaseDir", 42 );

    const auto settings = SimuContextDetail::ReadCaseIOSettings( parameters );

    EXPECT_EQ( settings.meshInputRoot, "../explicit-grid" );
}

TEST( SimuContextParameterTest, EmptyLegacySourceGridPreservesMeshRoot )
{
    EXPECT_EQ(
        SimuContextDetail::ApplyLegacySourceGridCaseDir( "", "" ),
        "" );
}

TEST( SimuContextParameterTest, DefaultLegacySourceGridPreservesMeshRoot )
{
    EXPECT_EQ(
        SimuContextDetail::ApplyLegacySourceGridCaseDir(
            "default", "default" ),
        "default" );
}

TEST( SimuContextParameterTest, ExplicitMeshRootWinsOverLegacySourceGrid )
{
    EXPECT_EQ(
        SimuContextDetail::ApplyLegacySourceGridCaseDir(
            "mesh/input", "../grid-case" ),
        "mesh/input" );
}

TEST( SimuContextParameterTest, CaseMeshRootWinsOverLegacySourceGrid )
{
    EXPECT_EQ(
        SimuContextDetail::ApplyLegacySourceGridCaseDir(
            "case", "../grid-case" ),
        "case" );
}


TEST( SimuContextIOTest, DefaultRootsWriteDiagnosticUnderResults )
{
    const std::filesystem::path caseDir =
        std::filesystem::temp_directory_path() / "OneFLOW_CaseIOSettings_Default";
    std::filesystem::remove_all( caseDir );
    std::filesystem::create_directories( caseDir );

    Prj::SetPrjBaseDir( caseDir.string() );
    SimuContextDetail::CaseIOSettings settings;
    SimuContextDetail::ApplyCaseIOSettings( settings );

    EXPECT_TRUE( std::filesystem::exists( caseDir / "results" / "io_paths.txt" ) );

    Prj::ClearPrjBaseDir();
    std::filesystem::remove_all( caseDir );
}

TEST( SimuContextIOTest, ExplicitDiagnosticNameUsesMeshAndFlowOutputRoots )
{
    const std::filesystem::path caseDir =
        std::filesystem::temp_directory_path() / "OneFLOW_CaseIOSettings_DomainRoots";
    std::filesystem::remove_all( caseDir );
    std::filesystem::create_directories( caseDir / "mesh-output" / "diagnostics" );
    std::filesystem::create_directories( caseDir / "flow-output" / "diagnostics" );

    Prj::SetPrjBaseDir( caseDir.string() );
    SimuContextDetail::CaseIOSettings settings;
    settings.meshOutputRoot = "mesh-output";
    settings.flowOutputRoot = "flow-output";
    settings.ioPathInfoFileName = "diagnostics/paths.txt";
    SimuContextDetail::ApplyCaseIOSettings( settings );

    EXPECT_TRUE(
        std::filesystem::exists(
            caseDir / "mesh-output" / "diagnostics" / "paths.txt" ) );
    EXPECT_TRUE(
        std::filesystem::exists(
            caseDir / "flow-output" / "diagnostics" / "paths.txt" ) );

    const std::string meshDiagnostics = ReadTextFile(
        caseDir / "mesh-output" / "diagnostics" / "paths.txt" );
    const std::string flowDiagnostics = ReadTextFile(
        caseDir / "flow-output" / "diagnostics" / "paths.txt" );
    EXPECT_NE(
        meshDiagnostics.find(
            "mesh_output_root=" + Prj::GetMeshOutputBaseDir() ),
        std::string::npos );
    EXPECT_NE(
        flowDiagnostics.find(
            "flow_output_root=" + Prj::GetFlowOutputBaseDir() ),
        std::string::npos );

    Prj::ClearPrjBaseDir();
    std::filesystem::remove_all( caseDir );
}
