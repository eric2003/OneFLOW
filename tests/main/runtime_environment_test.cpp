#include "RuntimeEnvironment.h"
#include "RuntimeVersionCompatibility.h"

#include <gtest/gtest.h>
#include <nlohmann/json.hpp>

#include <filesystem>
#include <fstream>
#include <iterator>
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
    ASSERT_TRUE(cgns.contains("version_compatibility"));
    const auto & cgnsCompatibility = cgns.at("version_compatibility");
    EXPECT_TRUE(cgnsCompatibility.at("status").is_string());
    EXPECT_TRUE(cgnsCompatibility.at("basis").is_string());
    EXPECT_TRUE(cgnsCompatibility.at("detail").is_string());
    const auto & cgnsRuntimeVersion = cgns.at("runtime_version");
    const auto & cgnsHeaderVersion = cgns.at("header_version");
    if (cgnsHeaderVersion.is_object())
    {
        EXPECT_TRUE(
            cgnsHeaderVersion.at("number").is_number_integer() ||
            cgnsHeaderVersion.at("number").is_null());
        EXPECT_TRUE(cgnsHeaderVersion.at("source").is_string());
    }
    const bool hasCgnsHeader =
        cgnsHeaderVersion.is_object() &&
        cgnsHeaderVersion.contains("value") &&
        cgnsHeaderVersion.at("value") != "unknown";
    if (!hasCgnsHeader)
    {
        EXPECT_EQ(cgnsCompatibility.at("status"), "unknown");
    }
    else if (cgnsRuntimeVersion.is_object() &&
             cgnsRuntimeVersion.at("value").is_null())
    {
        EXPECT_EQ(cgnsCompatibility.at("status"), "assumed_match");
        EXPECT_EQ(cgnsCompatibility.at("basis"), "header_version");
    }
    else if (cgnsRuntimeVersion.is_object() &&
             !cgnsRuntimeVersion.at("value").is_null())
    {
        EXPECT_TRUE(
            cgnsCompatibility.at("status") == "verified_match" ||
            cgnsCompatibility.at("status") == "mismatch");
    }

    const auto & hdf5 = document.at("dependencies").at(1);
    const auto & hdf5Version = hdf5.at("runtime_version");
    EXPECT_EQ(hdf5.at("name"), "HDF5");
    EXPECT_TRUE(hdf5.at("enabled").is_boolean());
    EXPECT_TRUE(hdf5.at("status").is_string());
    EXPECT_TRUE(hdf5.contains("runtime_version"));
    EXPECT_TRUE(hdf5.contains("runtime_version_status"));
    EXPECT_TRUE(hdf5.contains("header_version"));
    EXPECT_TRUE(hdf5.contains("runtime_library"));
    ASSERT_TRUE(hdf5.contains("version_compatibility"));
    const auto & hdf5Compatibility = hdf5.at("version_compatibility");
    EXPECT_TRUE(hdf5Compatibility.at("status").is_string());
    EXPECT_TRUE(hdf5Compatibility.at("basis").is_string());
    EXPECT_TRUE(hdf5Compatibility.at("detail").is_string());
    EXPECT_TRUE(hdf5.at("header_version").is_null());
    EXPECT_EQ(hdf5Compatibility.at("status"), "unknown");
    if (hdf5.at("enabled").get<bool>())
    {
        EXPECT_TRUE(
            hdf5.at("status") == "available" ||
            hdf5.at("status") == "query_failed");
        if (hdf5.at("status") == "available")
        {
            EXPECT_TRUE(hdf5Version.at("value").is_string());
            EXPECT_EQ(hdf5Version.at("source"), "H5get_libversion");
        }
        else
        {
            EXPECT_TRUE(hdf5Version.at("value").is_null());
        }
    }
    else
    {
        EXPECT_EQ(hdf5.at("status"), "not_detected");
        EXPECT_TRUE(hdf5Version.at("value").is_null());
    }
}

TEST(RuntimeVersionCompatibilityTest, VerifiesMatchingIndependentVersions)
{
    const auto result =
        ONEFLOW::AssessVersionCompatibility("4.5.2", "4.5.2");
    EXPECT_EQ(result.status, "verified_match");
    EXPECT_EQ(result.basis, "header_and_runtime_versions");
}

TEST(RuntimeVersionCompatibilityTest, DetectsDifferentIndependentVersions)
{
    const auto result =
        ONEFLOW::AssessVersionCompatibility("4.5.2", "4.5.1");
    EXPECT_EQ(result.status, "mismatch");
    EXPECT_EQ(result.basis, "header_and_runtime_versions");
}

TEST(RuntimeVersionCompatibilityTest, AssumesMatchWhenOnlyHeaderVersionExists)
{
    const auto result =
        ONEFLOW::AssessVersionCompatibility("4.5.2", "");
    EXPECT_EQ(result.status, "assumed_match");
    EXPECT_EQ(result.basis, "header_version");
    EXPECT_NE(result.detail.find("ABI compatibility"), std::string::npos);
}

TEST(RuntimeVersionCompatibilityTest, RuntimeVersionAloneIsNotHeaderEvidence)
{
    const auto result =
        ONEFLOW::AssessVersionCompatibility("", "1.14.6");
    EXPECT_EQ(result.status, "unknown");
    EXPECT_EQ(result.basis, "insufficient_version_evidence");
}

TEST(RuntimeVersionCompatibilityTest, MissingAndUnknownStringsAreNotEvidence)
{
    EXPECT_EQ(
        ONEFLOW::AssessVersionCompatibility("", "").status,
        "unknown");
    EXPECT_EQ(
        ONEFLOW::AssessVersionCompatibility("unknown", "1.2.3").status,
        "unknown");
    EXPECT_EQ(
        ONEFLOW::AssessVersionCompatibility("1.2.3", "unknown").status,
        "assumed_match");
}

TEST(RuntimeEnvironmentSnapshotTest, WritesValidJsonAndCreatesParentDirectory)
{
    const std::filesystem::path root =
        std::filesystem::temp_directory_path()
        / "oneflow_runtime_environment_snapshot_test";
    const std::filesystem::path report =
        root / "results" / "runtime_environment.json";
    std::error_code cleanupError;
    std::filesystem::remove_all(root, cleanupError);

    std::string error;
    ASSERT_TRUE(ONEFLOW::WriteRuntimeEnvironmentJsonFile(
        report.string(), error)) << error;
    ASSERT_TRUE(std::filesystem::is_regular_file(report));

    std::ifstream input(report);
    ASSERT_TRUE(input.good());
    const std::string contents{
        std::istreambuf_iterator<char>(input),
        std::istreambuf_iterator<char>()};
    const nlohmann::json document = nlohmann::json::parse(contents);
    EXPECT_EQ(document.at("schema_version"), 1);
    EXPECT_TRUE(document.at("oneflow").contains("version"));
    EXPECT_TRUE(document.at("dependencies").is_array());

    std::filesystem::remove_all(root, cleanupError);
}

TEST(RuntimeEnvironmentCommandTest, ReportsHdf5RuntimeVersionInText)
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
    EXPECT_NE(report.find("Runtime version :"), std::string::npos);
    EXPECT_NE(report.find("Version consistency:"), std::string::npos);
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

