#include "RuntimeEnvironment.h"

#include <gtest/gtest.h>
#include <nlohmann/json.hpp>

#include <sstream>
#include <string>

TEST(RuntimeEnvironmentCommandTest, EmitsMachineReadableDependencyInformation)
{
    char executable[] = "OneFLOW";
    char environment[] = "--environment";
    char format[] = "--format";
    char jsonFormat[] = "json";
    char * arguments[] = {
        executable, environment, format, jsonFormat
    };

    std::ostringstream output;
    std::ostringstream error;
    const int status = ONEFLOW::HandleRuntimeEnvironmentCommand(
        4, arguments, output, error);

    ASSERT_EQ(status, 0) << error.str();
    const nlohmann::json document = nlohmann::json::parse(output.str());

    EXPECT_EQ(document.at("schema_version"), 1);
    EXPECT_TRUE(document.at("oneflow").at("version").is_string());
    EXPECT_TRUE(document.at("build").at("compiler").is_string());
    ASSERT_TRUE(document.at("dependencies").is_array());
    ASSERT_EQ(document.at("dependencies").size(), 2);

    const auto & cgns = document.at("dependencies").at(0);
    EXPECT_EQ(cgns.at("name"), "CGNS");
    EXPECT_TRUE(cgns.at("enabled").is_boolean());
    EXPECT_TRUE(cgns.contains("runtime_version"));
    EXPECT_TRUE(cgns.contains("header_version"));
    EXPECT_TRUE(cgns.contains("runtime_library"));

    const auto & hdf5 = document.at("dependencies").at(1);
    EXPECT_EQ(hdf5.at("name"), "HDF5");
    EXPECT_TRUE(hdf5.at("enabled").is_boolean());
    EXPECT_TRUE(hdf5.at("status").is_string());
    EXPECT_TRUE(hdf5.contains("runtime_version"));
    EXPECT_TRUE(hdf5.contains("runtime_version_status"));
    EXPECT_TRUE(hdf5.contains("header_version"));
    EXPECT_TRUE(hdf5.contains("runtime_library"));
}

TEST(RuntimeEnvironmentCommandTest, ReportsUnavailableHdf5RuntimeVersionInText)
{
    char executable[] = "OneFLOW";
    char environment[] = "--environment";
    char * arguments[] = { executable, environment };

    std::ostringstream output;
    std::ostringstream error;
    const int status = ONEFLOW::HandleRuntimeEnvironmentCommand(
        2, arguments, output, error);

    ASSERT_EQ(status, 0) << error.str();
    const std::string report = output.str();
    EXPECT_NE(report.find("HDF5"), std::string::npos);
    EXPECT_NE(report.find("Runtime version : unavailable"), std::string::npos);
    EXPECT_NE(report.find("runtime_api_symbol_not_found"), std::string::npos);
}

TEST(RuntimeEnvironmentCommandTest, PrintsOneFlowVersion)
{
    char executable[] = "OneFLOW";
    char version[] = "--version";
    char * arguments[] = { executable, version };

    std::ostringstream output;
    std::ostringstream error;
    const int status = ONEFLOW::HandleRuntimeEnvironmentCommand(
        2, arguments, output, error);

    EXPECT_EQ(status, 0);
    EXPECT_NE(output.str().find("OneFLOW"), std::string::npos);
    EXPECT_TRUE(error.str().empty());
}

TEST(RuntimeEnvironmentCommandTest, RejectsUnsupportedEnvironmentFormat)
{
    char executable[] = "OneFLOW";
    char environment[] = "--environment";
    char format[] = "--format";
    char unsupported[] = "xml";
    char * arguments[] = {
        executable, environment, format, unsupported
    };

    std::ostringstream output;
    std::ostringstream error;
    const int status = ONEFLOW::HandleRuntimeEnvironmentCommand(
        4, arguments, output, error);

    EXPECT_EQ(status, 2);
    EXPECT_NE(error.str().find("Usage:"), std::string::npos);
    EXPECT_TRUE(output.str().empty());
}

TEST(RuntimeEnvironmentCommandTest, LeavesSimulationArgumentsUnchanged)
{
    char executable[] = "OneFLOW";
    char runMode[] = "r";
    char caseDirectory[] = "case";
    char * arguments[] = { executable, runMode, caseDirectory };

    std::ostringstream output;
    std::ostringstream error;
    const int status = ONEFLOW::HandleRuntimeEnvironmentCommand(
        3, arguments, output, error);

    EXPECT_EQ(status, -1);
    EXPECT_TRUE(output.str().empty());
    EXPECT_TRUE(error.str().empty());
}

TEST(RuntimeEnvironmentCommandTest, RejectsExtraVersionArguments)
{
    char executable[] = "OneFLOW";
    char version[] = "--version";
    char extra[] = "unexpected";
    char * arguments[] = { executable, version, extra };

    std::ostringstream output;
    std::ostringstream error;
    const int status = ONEFLOW::HandleRuntimeEnvironmentCommand(
        3, arguments, output, error);

    EXPECT_EQ(status, 2);
    EXPECT_NE(error.str().find("Usage:"), std::string::npos);
}

