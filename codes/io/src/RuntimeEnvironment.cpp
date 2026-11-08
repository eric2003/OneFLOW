/*---------------------------------------------------------------------------*\\
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
----------------------------------------------------------------------------*/

#include "RuntimeEnvironment.h"
#include "RuntimeVersionCompatibility.h"

#include <nlohmann/json.hpp>

#include <cstdint>
#include <filesystem>
#include <exception>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>
#elif defined(__unix__) || defined(__APPLE__)
#include <dlfcn.h>
#endif

#if defined(ENABLE_CGNS)
#include <cgnslib.h>
#endif

#ifndef ONEFLOW_PROJECT_VERSION
#define ONEFLOW_PROJECT_VERSION "unknown"
#endif

#ifndef ONEFLOW_COMPILER_ID
#define ONEFLOW_COMPILER_ID "unknown"
#endif

#ifndef ONEFLOW_COMPILER_VERSION
#define ONEFLOW_COMPILER_VERSION "unknown"
#endif

#ifndef ONEFLOW_BUILD_CONFIGURATION
#define ONEFLOW_BUILD_CONFIGURATION "unknown"
#endif

namespace ONEFLOW
{
namespace
{

using Json = nlohmann::json;

std::string GetVersionValue(const Json & version)
{
    if (!version.is_object() ||
        !version.contains("value") ||
        !version["value"].is_string())
    {
        return {};
    }
    return version["value"].get<std::string>();
}

Json BuildVersionCompatibility(
    const Json & headerVersion,
    const Json & runtimeVersion)
{
    const VersionCompatibilityAssessment assessment =
        AssessVersionCompatibility(
            GetVersionValue(headerVersion),
            GetVersionValue(runtimeVersion));
    return {
        {"status", assessment.status},
        {"basis", assessment.basis},
        {"detail", assessment.detail}
    };
}

#if defined(ENABLE_CGNS)

struct CgnsModuleInfo
{
    std::string path;
    std::string linkage = "unresolved";
};

std::string FormatCgnsVersionNumber(int value)
{
    if (value <= 0)
    {
        return {};
    }

    const int major = value / 1000;
    const int minor = (value / 100) % 10;
    const int patch = (value / 10) % 10;
    const int revision = value % 10;

    std::ostringstream stream;
    stream << major << '.' << minor << '.' << patch;
    if (revision != 0)
    {
        stream << '.' << revision;
    }
    return stream.str();
}

CgnsModuleInfo QueryCgnsModule()
{
    CgnsModuleInfo result;

#if defined(_WIN32)
    HMODULE module = nullptr;
    const auto functionAddress = reinterpret_cast<LPCSTR>(
        reinterpret_cast<std::uintptr_t>(&cg_open));
    const DWORD flags =
        GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS |
        GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT;

    if (GetModuleHandleExA(flags, functionAddress, &module) != 0)
    {
        char pathBuffer[32768] = {};
        const DWORD length = GetModuleFileNameA(
            module, pathBuffer, static_cast<DWORD>(sizeof(pathBuffer)));
        if (length > 0 && length < sizeof(pathBuffer))
        {
            result.path.assign(pathBuffer, length);
            result.linkage = module == GetModuleHandleA(nullptr)
                ? "linked_into_executable"
                : "shared";
        }
    }
#elif defined(__unix__) || defined(__APPLE__)
    Dl_info moduleInfo{};
    const auto functionAddress = reinterpret_cast<void *>(
        reinterpret_cast<std::uintptr_t>(&cg_open));
    if (dladdr(functionAddress, &moduleInfo) != 0 &&
        moduleInfo.dli_fname != nullptr)
    {
        result.path = moduleInfo.dli_fname;
        const std::string::size_type sharedObject = result.path.find(".so");
        const std::string::size_type dynamicLibrary = result.path.find(".dylib");
        result.linkage =
            sharedObject != std::string::npos ||
            dynamicLibrary != std::string::npos
            ? "shared"
            : "linked_into_executable";
    }
#endif

    return result;
}

Json QueryCgnsDependency()
{
    const CgnsModuleInfo module = QueryCgnsModule();

    Json runtimeVersion = {
        {"value", nullptr},
        {"source", "unavailable"}
    };
    std::string runtimeStatus = "runtime_version_api_unavailable";

#if defined(ONEFLOW_CGNS_HAS_RUNTIME_VERSION_API)
    int versionNumber = 0;
    char * versionString = nullptr;
    cg_get_cgns_version(&versionNumber);
    cg_get_cgns_version_string(&versionString);

    std::string value;
    if (versionString != nullptr && versionString[0] != '\0')
    {
        value = versionString;
    }
    else
    {
        value = FormatCgnsVersionNumber(versionNumber);
    }

    if (!value.empty())
    {
        runtimeVersion = {
            {"value", value},
            {"source", "runtime_api"}
        };
        runtimeStatus = "available";
    }
    else
    {
        runtimeStatus = "runtime_version_query_failed";
    }
#endif

    Json headerVersion = {
        {"value", nullptr},
        {"number", nullptr},
        {"source", "unavailable"}
    };
#ifdef CGNS_VERSION
    headerVersion = {
        {"value", FormatCgnsVersionNumber(CGNS_VERSION)},
        {"number", CGNS_VERSION},
        {"source", "CGNS_VERSION_compile_time_macro"}
    };
#endif

    return {
        {"name", "CGNS"},
        {"enabled", true},
        {"status", "available"},
        {"runtime_version", runtimeVersion},
        {"runtime_version_status", runtimeStatus},
        {"header_version", headerVersion},
        {"runtime_library", {
            {"module_path", module.path},
            {"linkage", module.linkage}
        }},
        {"version_compatibility",
            BuildVersionCompatibility(headerVersion, runtimeVersion)}
    };
}

#else

Json QueryCgnsDependency()
{
    return {
        {"name", "CGNS"},
        {"enabled", false},
        {"status", "disabled"},
        {"runtime_version", nullptr},
        {"header_version", nullptr},
        {"version_compatibility",
            BuildVersionCompatibility(nullptr, nullptr)},
        {"runtime_library", nullptr}
    };
}

#endif

struct RuntimeSymbolInfo
{
    std::string modulePath;
    std::string linkage = "unavailable";
};

using Hdf5GetLibVersion = int (*)(unsigned int *, unsigned int *, unsigned int *);

RuntimeSymbolInfo FindHdf5VersionSymbol(Hdf5GetLibVersion & function)
{
    RuntimeSymbolInfo result;
    function = nullptr;

#if defined(_WIN32)
    HANDLE process = GetCurrentProcess();
    std::vector<HMODULE> modules(256);
    DWORD bytesNeeded = 0;
    if (!EnumProcessModules(
            process, modules.data(),
            static_cast<DWORD>(modules.size() * sizeof(HMODULE)),
            &bytesNeeded))
    {
        return result;
    }

    if (bytesNeeded > modules.size() * sizeof(HMODULE))
    {
        modules.resize(bytesNeeded / sizeof(HMODULE));
        if (!EnumProcessModules(
                process, modules.data(),
                static_cast<DWORD>(modules.size() * sizeof(HMODULE)),
                &bytesNeeded))
        {
            return result;
        }
    }

    const std::size_t moduleCount =
        bytesNeeded / sizeof(HMODULE);
    for (std::size_t i = 0; i < moduleCount; ++i)
    {
        HMODULE module = modules[i];
        FARPROC symbol = GetProcAddress(module, "H5get_libversion");
        if (symbol == nullptr)
        {
            continue;
        }

        function = reinterpret_cast<Hdf5GetLibVersion>(symbol);
        char pathBuffer[32768] = {};
        const DWORD length = GetModuleFileNameA(
            module, pathBuffer, static_cast<DWORD>(sizeof(pathBuffer)));
        if (length > 0 && length < sizeof(pathBuffer))
        {
            result.modulePath.assign(pathBuffer, length);
        }
        result.linkage = module == GetModuleHandleA(nullptr)
            ? "linked_into_executable"
            : "shared";
        break;
    }
#elif defined(__unix__) || defined(__APPLE__)
    void * symbol = dlsym(RTLD_DEFAULT, "H5get_libversion");
    if (symbol != nullptr)
    {
        function = reinterpret_cast<Hdf5GetLibVersion>(symbol);
        Dl_info moduleInfo{};
        if (dladdr(symbol, &moduleInfo) != 0 && moduleInfo.dli_fname != nullptr)
        {
            result.modulePath = moduleInfo.dli_fname;
            result.linkage = "shared";
        }
    }
#endif

    return result;
}

Json QueryHdf5Dependency()
{
    const Json headerVersion = nullptr;
    Hdf5GetLibVersion getVersion = nullptr;
    const RuntimeSymbolInfo module = FindHdf5VersionSymbol(getVersion);
    if (getVersion == nullptr)
    {
        return {
            {"name", "HDF5"},
            {"enabled", false},
            {"status", "not_detected"},
            {"runtime_version", {
                {"value", nullptr},
                {"source", "unavailable"}
            }},
            {"runtime_version_status", "runtime_api_symbol_not_found"},
            {"header_version", headerVersion},
            {"version_compatibility",
                BuildVersionCompatibility(headerVersion, nullptr)},
            {"runtime_library", {
                {"module_path", module.modulePath},
                {"linkage", module.linkage}
            }},
            {"detail", "H5get_libversion was not exported by a loaded process module."}
        };
    }

    unsigned int major = 0;
    unsigned int minor = 0;
    unsigned int release = 0;
    const int status = getVersion(&major, &minor, &release);
    if (status < 0)
    {
        return {
            {"name", "HDF5"},
            {"enabled", true},
            {"status", "query_failed"},
            {"runtime_version", {
                {"value", nullptr},
                {"source", "unavailable"}
            }},
            {"runtime_version_status", "runtime_api_call_failed"},
            {"header_version", headerVersion},
            {"version_compatibility",
                BuildVersionCompatibility(headerVersion, nullptr)},
            {"runtime_library", {
                {"module_path", module.modulePath},
                {"linkage", module.linkage}
            }}
        };
    }

    std::ostringstream version;
    version << major << '.' << minor << '.' << release;
    return {
        {"name", "HDF5"},
        {"enabled", true},
        {"status", "available"},
        {"runtime_version", {
            {"value", version.str()},
            {"source", "H5get_libversion"}
        }},
        {"runtime_version_status", "available"},
        {"header_version", headerVersion},
        {"version_compatibility",
            BuildVersionCompatibility(nullptr, Json{{"value", version.str()}})},
        {"runtime_library", {
            {"module_path", module.modulePath},
            {"linkage", module.linkage}
        }}
    };
}

Json QueryRuntimeEnvironment()
{
    Json dependencies = Json::array();
    dependencies.push_back(QueryCgnsDependency());
    dependencies.push_back(QueryHdf5Dependency());

    return {
        {"schema_version", 1},
        {"oneflow", {
            {"version", ONEFLOW_PROJECT_VERSION}
        }},
        {"build", {
            {"compiler", ONEFLOW_COMPILER_ID},
            {"compiler_version", ONEFLOW_COMPILER_VERSION},
            {"configuration", ONEFLOW_BUILD_CONFIGURATION}
        }},
        {"dependencies", dependencies}
    };
}

void WriteRuntimeEnvironmentText(
    const Json & environment,
    std::ostream & output)
{
    output << "OneFLOW Runtime Environment\n";
    output << "----------------------------------------\n";
    output << "OneFLOW       : "
           << environment["oneflow"]["version"].get<std::string>() << '\n';
    output << "Compiler      : "
           << environment["build"]["compiler"].get<std::string>() << ' '
           << environment["build"]["compiler_version"].get<std::string>() << '\n';
    output << "Configuration : "
           << environment["build"]["configuration"].get<std::string>() << '\n';
    output << "\nDependencies\n";
    output << "----------------------------------------\n";

    for (const Json & dependency : environment["dependencies"])
    {
        output << dependency["name"].get<std::string>() << '\n';
        output << "  Status          : "
               << dependency["status"].get<std::string>() << '\n';

        const Json & runtimeVersion = dependency["runtime_version"];
        if (!runtimeVersion.is_null() &&
            !runtimeVersion["value"].is_null())
        {
            output << "  Runtime version : "
                   << runtimeVersion["value"].get<std::string>()
                   << " (" << runtimeVersion["source"].get<std::string>()
                   << ")\n";
        }
        else
        {
            output << "  Runtime version : unavailable\n";
        }

        const Json & headerVersion = dependency["header_version"];
        if (!headerVersion.is_null())
        {
            output << "  Header version  : "
                   << headerVersion["value"].get<std::string>()
                   << " (" << headerVersion.value("source", "compile-time")
                   << ")\n";
            if (headerVersion.contains("header_path"))
            {
                output << "  Header path     : "
                       << headerVersion["header_path"].get<std::string>()
                       << '\n';
            }
        }

        if (dependency.contains("version_compatibility"))
        {
            const Json & compatibility = dependency["version_compatibility"];
            output << "  Version consistency: "
                   << compatibility["status"].get<std::string>()
                   << " (" << compatibility["basis"].get<std::string>()
                   << ")\n";
            output << "  Version detail     : "
                   << compatibility["detail"].get<std::string>() << '\n';
        }

        const Json & runtimeLibrary = dependency["runtime_library"];
        if (!runtimeLibrary.is_null())
        {
            output << "  Linkage         : "
                   << runtimeLibrary["linkage"].get<std::string>() << '\n';
            output << "  Loaded module   : "
                   << runtimeLibrary["module_path"].get<std::string>() << '\n';
        }

        if (dependency.contains("runtime_version_status"))
        {
            output << "  Query status    : "
                   << dependency["runtime_version_status"].get<std::string>()
                   << '\n';
        }
        if (dependency.contains("detail"))
        {
            output << "  Detail          : "
                   << dependency["detail"].get<std::string>() << '\n';
        }
    }
}

} // namespace

bool WriteRuntimeEnvironmentJsonFile(
    const std::string & fileName,
    std::string & error)
{
    error.clear();

    try
    {
        const std::filesystem::path outputPath(fileName);
        if (outputPath.empty())
        {
            error = "output file path is empty";
            return false;
        }

        const std::filesystem::path parentPath = outputPath.parent_path();
        if (!parentPath.empty())
        {
            std::error_code filesystemError;
            std::filesystem::create_directories(parentPath, filesystemError);
            if (filesystemError)
            {
                error = "could not create output directory: " +
                    filesystemError.message();
                return false;
            }
        }

        std::ofstream output(
            outputPath, std::ios_base::out | std::ios_base::trunc);
        if (!output)
        {
            error = "could not open runtime environment report";
            return false;
        }

        output << QueryRuntimeEnvironment().dump(2) << '\n';
        output.flush();
        if (!output)
        {
            error = "failed to write runtime environment report";
            return false;
        }

        output.close();
        if (!output)
        {
            error = "failed to close runtime environment report";
            return false;
        }
        return true;
    }
    catch (const std::exception & exception)
    {
        error = exception.what();
        return false;
    }
    catch (...)
    {
        error = "unknown error while writing runtime environment report";
        return false;
    }
}

int HandleRuntimeEnvironmentCommand(
    int argc,
    char ** argv,
    std::ostream & output,
    std::ostream & error)
{
    if (argc < 2 || argv == nullptr || argv[1] == nullptr)
    {
        return -1;
    }

    const std::string command = argv[1];
    if (command == "--version")
    {
        if (argc != 2)
        {
            error << "Usage: OneFLOW --version\n";
            return 2;
        }
        output << "OneFLOW " << ONEFLOW_PROJECT_VERSION << '\n';
        return 0;
    }

    if (command != "--environment")
    {
        return -1;
    }

    bool jsonFormat = false;
    if (argc == 4 &&
        argv[2] != nullptr &&
        argv[3] != nullptr &&
        std::string(argv[2]) == "--format" &&
        std::string(argv[3]) == "json")
    {
        jsonFormat = true;
    }
    else if (argc != 2)
    {
        error << "Usage: OneFLOW --environment [--format json]\n";
        return 2;
    }

    const Json environment = QueryRuntimeEnvironment();
    if (jsonFormat)
    {
        output << environment.dump(2) << '\n';
    }
    else
    {
        WriteRuntimeEnvironmentText(environment, output);
    }
    return 0;
}

} // namespace ONEFLOW
