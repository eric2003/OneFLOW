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

#include "FileUtils.h"
#include "Fatal.h"

#ifdef _WINDOWS
#include <windows.h>
#include <direct.h>
#include <io.h>
#else
    #ifdef WIN_GNU
        #include <windows.h>
        #include <direct.h>
        #include <io.h>
    #else
        #include <sys/stat.h>
        #include <unistd.h>
    #endif
#endif

#include <vector>
#include <climits>
#include <iostream>
#include <filesystem>

BeginNameSpace( ONEFLOW )

bool HX_IsDirectory(const std::string& dirName)
{
    std::error_code ec;
    return std::filesystem::is_directory(dirName, ec) && !ec;
}

bool HX_CreateDirectory( const std::string & dirName )
{
    try
    {
        if ( std::filesystem::exists( dirName ) )
        {
            if ( std::filesystem::is_directory( dirName ) )
            {
                std::cout << "Directory already exists: "
                    << dirName << std::endl;
                return true;
            }
            else
            {
                std::cerr << "Path exists but is not a directory: "
                    << dirName << std::endl;
                return false;
            }
        }

        // Create the directory recursively with default permissions.
        if ( std::filesystem::create_directories( dirName ) )
        {
            std::cout << "Directory created successfully: "
                << dirName << std::endl;
            return true;
        }
        else
        {
            std::cerr << "Failed to create directory: "
                << dirName << std::endl;
            return false;
        }
    }
    catch ( const std::filesystem::filesystem_error & e )
    {
        std::cerr << "Filesystem error: "
            << e.what() << std::endl;
        return false;
    }
}

std::string HX_GetExeDirectory()
{
#ifdef _WIN32
    std::wstring wbuf(MAX_PATH, L'\0');
    DWORD len = 0;
    while (true) {
        len = GetModuleFileNameW(nullptr, wbuf.data(), static_cast<DWORD>(wbuf.size()));
        if (len == 0) return {};
        if (len < wbuf.size()) break;
        wbuf.resize(wbuf.size() * 2);
    }
    return std::filesystem::path(wbuf.substr(0, len)).parent_path().string();
#else
    std::vector<char> buf(PATH_MAX);
    ssize_t count = -1;
    while (true) {
        count = readlink("/proc/self/exe", buf.data(), buf.size());
        if (count < 0) return {};
        if (static_cast<size_t>(count) < buf.size()) break;
        buf.resize(buf.size() * 2);
    }
    buf[count] = '\0';
    return std::filesystem::path(std::string(buf.data(), count)).parent_path().string();
#endif
}


std::string HX_GetCurrentDirectory()
{
    try {
        return std::filesystem::current_path().string();
    } catch (const std::filesystem::filesystem_error&) {
        return {};  // or rethrow / log the error
    }
}

bool EndWithBackwardSlash( const std::string & fileName )
{
    return !fileName.empty() && fileName.back() == '\\';
}

bool EndWithForwardSlash( const std::string & fileName )
{
    return !fileName.empty() && fileName.back() == '/';
}

bool EndWithSlash( const std::string & fileName )
{
    return EndWithForwardSlash( fileName ) ||
        EndWithBackwardSlash( fileName );
}

bool StartWithForwardSlash( const std::string & fileName )
{
    return !fileName.empty() && fileName.front() == '/';
}

std::string RemoveFirstSlash( const std::string & fileName )
{
    if ( !fileName.empty() &&
        ( fileName.front() == '/' || fileName.front() == '\\' ) )
    {
        return fileName.substr( 1 );
    }

    return fileName;
}

std::string RemoveEndSlash( const std::string & fileName )
{
    if ( !fileName.empty() &&
        ( fileName.back() == '/' || fileName.back() == '\\' ) )
    {
        return fileName.substr( 0, fileName.size() - 1 );
    }

    return fileName;
}

void GetFileNameExtension( const std::string & fullName, std::string & mainName, std::string & extensionName, const std::string & fileNameSeparator )
{
    std::basic_string <char>::size_type index;

    index         = fullName.find_last_of( fileNameSeparator );
    mainName      = fullName.substr( 0, index );
    extensionName = fullName.substr( index+1, fullName.length() - index - 1 );
}

void ModifyFileMainName( std::string & fileName,  const std::string & newMainName )
{
    std::string mainName, extensionName;
    ONEFLOW::GetFileNameExtension( fileName, mainName, extensionName, "." );

    std::ostringstream oss;
    oss << newMainName << "." << extensionName;

    fileName = oss.str();
}

void ModifyFileExtensionName( std::string & fileName,  const std::string & newExtensionName )
{
    std::string mainName, extensionName;
    ONEFLOW::GetFileNameExtension( fileName, mainName, extensionName, "." );

    std::ostringstream oss;
    oss << mainName << "." << newExtensionName;

    fileName = oss.str();
}

EndNameSpace
