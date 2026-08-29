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

#include "FileUtil.h"
#include "Stop.h"

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

#include <iostream>
#include <filesystem>

BeginNameSpace( ONEFLOW )

bool DirExist(const std::string& dirName)
{
    std::error_code ec;
    return std::filesystem::is_directory(dirName, ec) && !ec;
}

bool MakeDir(const std::string& dirName)
{
    try {
        if (std::filesystem::exists(dirName)) {
            if (std::filesystem::is_directory(dirName)) {
                std::cout << "Directory already exists: " << dirName << std::endl;
                return true;
            } else {
                std::cerr << "Path exists but is not a directory: " << dirName << std::endl;
                return false;
            }
        }

        // 创建目录（默认权限，支持递归）
        if (std::filesystem::create_directories(dirName)) {
            std::cout << "Directory created successfully: " << dirName << std::endl;
            return true;
        } else {
            std::cerr << "Failed to create directory: " << dirName << std::endl;
            return false;
        }
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
        return false;
    }
}

//std::string HX_GetExeDirectory()
//{
//    char buffer[ FILENAME_MAX ] = { 0 };
//#ifdef _WIN32
//    GetModuleFileName( NULL, buffer, FILENAME_MAX );
//#else
//    std::size_t count = readlink( "/proc/self/exe", buffer, FILENAME_MAX );
//#endif
//    std::string::size_type pos = std::string( buffer ).find_last_of( "\\/" );
//    return std::string( buffer ).substr( 0, pos);
//}

std::filesystem::path HX_GetExeDirectory()
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
    return std::filesystem::path(wbuf.substr(0, len)).parent_path();
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
    return std::filesystem::path(std::string(buf.data(), count)).parent_path();
#endif
}


//std::string HX_GetCurrentDir()
//{
//#ifdef _WINDOWS
//    char * cwd = _getcwd( 0, 0 );
//#else
//    char * cwd = getcwd( 0, 0 ); 
//#endif
//    std::string working_dir( cwd ) ;
//    std::free( cwd ) ;
//    return working_dir ;
//}

std::string HX_GetCurrentDir()
{
    try {
        return std::filesystem::current_path().string();
    } catch (const std::filesystem::filesystem_error&) {
        return {};  // or rethrow / log the error
    }
}

bool EndWithSlash( const std::string & fileName )
{
    if ( EndWithForwardSlash( fileName ) ||
        EndWithBackwardSlash( fileName ) )
    {
        return true;
    }
    return false;
}

bool EndWithBackwardSlash( const std::string & fileName )
{
    size_t pos = fileName.find_last_of("\\");
    size_t ss = fileName.size();
    if ( ss == 0 )
    {
        return false;
    }
    else
    {
        bool flag = fileName.substr( ss - 1, 1 ) == "\\";
        return flag;
    }
}

bool EndWithForwardSlash( const std::string & fileName )
{
    size_t pos = fileName.find_last_of("/");
    size_t ss = fileName.size();
    if ( ss == 0 )
    {
        return false;
    }
    else
    {
        bool flag = fileName.substr( ss - 1, 1 ) == "/";
        return flag;
    }
}

bool StartWithForwardSlash( const std::string & fileName )
{
    size_t pos = fileName.find_first_of("/");
    if ( fileName.size() == 0 )
    {
        return false;
    }

    if ( fileName.substr( 0,1 ) == "/" )
    {
        return true;
    }
    return false;
}

std::string RemoveFirstSlash( const std::string & fileName )
{
    if ( StartWithForwardSlash( fileName ) )
    {
        int len = fileName.size();
        return fileName.substr( 1, len - 1 );
    }
    return fileName;
}

std::string RemoveEndSlash( const std::string & fileName )
{
    if ( EndWithSlash( fileName ) )
    {
        int len = fileName.size();
        return fileName.substr( 0, len - 1 );
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
