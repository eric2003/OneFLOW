/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software; you can redistribute it and/or modify it
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
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "ProcessRunner.h"

#include <cerrno>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <spawn.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
extern char ** environ;
#endif

BeginNameSpace( ONEFLOW )

namespace
{
#ifdef _WIN32
    std::wstring Utf8ToWide( const std::string & value )
    {
        if ( value.empty() ) return {};
        const int count = MultiByteToWideChar(
            CP_UTF8, MB_ERR_INVALID_CHARS, value.data(),
            static_cast< int >( value.size() ), nullptr, 0 );
        if ( count <= 0 )
            throw std::runtime_error( "Invalid UTF-8 process argument" );
        std::wstring result( static_cast< std::size_t >( count ), L'\0' );
        if ( MultiByteToWideChar(
                 CP_UTF8, MB_ERR_INVALID_CHARS, value.data(),
                 static_cast< int >( value.size() ), result.data(), count ) != count )
            throw std::runtime_error( "Failed converting UTF-8 process argument" );
        return result;
    }

    std::wstring QuoteWindowsArgument( const std::wstring & value )
    {
        if ( ! value.empty() &&
             value.find_first_of( L" \t\n\v\"" ) == std::wstring::npos )
            return value;

        std::wstring quoted = L"\"";
        std::size_t backslashes = 0;
        for ( wchar_t character : value )
        {
            if ( character == L'\\' )
            {
                ++backslashes;
                continue;
            }
            if ( character == L'\"' )
            {
                quoted.append( backslashes * 2 + 1, L'\\' );
                quoted += L'\"';
                backslashes = 0;
                continue;
            }
            quoted.append( backslashes, L'\\' );
            backslashes = 0;
            quoted += character;
        }
        quoted.append( backslashes * 2, L'\\' );
        quoted += L'\"';
        return quoted;
    }

    std::wstring ResolveWindowsExecutable( const std::wstring & executable )
    {
        if ( executable.find_first_of( L"\\/" ) != std::wstring::npos )
            return executable;
        std::vector< wchar_t > buffer( 32768, L'\0' );
        const DWORD length = SearchPathW(
            nullptr, executable.c_str(), L".exe",
            static_cast< DWORD >( buffer.size() ), buffer.data(), nullptr );
        if ( length == 0 || length >= buffer.size() )
            return executable;
        return std::wstring( buffer.data(), length );
    }
#endif
}

ProcessResult RunProcess(
    const std::string & executable,
    const std::vector< std::string > & arguments,
    const std::string & workingDirectory )
{
    if ( executable.empty() )
        throw std::invalid_argument( "Process executable must not be empty" );

    ProcessResult result;
#ifdef _WIN32
    const std::wstring wideExecutable = ResolveWindowsExecutable( Utf8ToWide( executable ) );
    std::wstring commandLine = QuoteWindowsArgument( wideExecutable );
    for ( const std::string & argument : arguments )
        commandLine += L" " + QuoteWindowsArgument( Utf8ToWide( argument ) );
    const std::wstring wideWorkingDirectory =
        workingDirectory.empty() ? std::wstring{} : Utf8ToWide( workingDirectory );

    SECURITY_ATTRIBUTES security{};
    security.nLength = sizeof( security );
    security.bInheritHandle = TRUE;
    HANDLE readPipe = nullptr;
    HANDLE writePipe = nullptr;
    if ( ! CreatePipe( &readPipe, &writePipe, &security, 0 ) )
    {
        result.errorMessage = "CreatePipe failed with Windows error " +
            std::to_string( GetLastError() );
        return result;
    }
    if ( ! SetHandleInformation( readPipe, HANDLE_FLAG_INHERIT, 0 ) )
    {
        result.errorMessage = "SetHandleInformation failed with Windows error " +
            std::to_string( GetLastError() );
        CloseHandle( readPipe );
        CloseHandle( writePipe );
        return result;
    }

    HANDLE inputHandle = CreateFileW(
        L"NUL", GENERIC_READ, FILE_SHARE_READ | FILE_SHARE_WRITE,
        &security, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr );
    if ( inputHandle == INVALID_HANDLE_VALUE )
    {
        result.errorMessage = "Opening NUL for child standard input failed with Windows error " +
            std::to_string( GetLastError() );
        CloseHandle( readPipe );
        CloseHandle( writePipe );
        return result;
    }

    STARTUPINFOW startup{};
    startup.cb = sizeof( startup );
    startup.dwFlags = STARTF_USESTDHANDLES;
    startup.hStdOutput = writePipe;
    startup.hStdError = writePipe;
    startup.hStdInput = inputHandle;
    PROCESS_INFORMATION process{};
    const BOOL created = CreateProcessW(
        wideExecutable.c_str(), commandLine.data(), nullptr, nullptr, TRUE,
        CREATE_NO_WINDOW, nullptr,
        wideWorkingDirectory.empty() ? nullptr : wideWorkingDirectory.c_str(),
        &startup, &process );
    CloseHandle( inputHandle );
    CloseHandle( writePipe );
    if ( ! created )
    {
        result.errorMessage = "CreateProcessW failed with Windows error " +
            std::to_string( GetLastError() ) + " for executable " + executable;
        CloseHandle( readPipe );
        return result;
    }

    result.started = true;
    char buffer[ 4096 ];
    DWORD bytesRead = 0;
    while ( ReadFile( readPipe, buffer, sizeof( buffer ), &bytesRead, nullptr ) &&
            bytesRead > 0 )
        result.output.append( buffer, static_cast< std::size_t >( bytesRead ) );
    CloseHandle( readPipe );
    WaitForSingleObject( process.hProcess, INFINITE );
    DWORD exitCode = 0;
    if ( GetExitCodeProcess( process.hProcess, &exitCode ) )
        result.exitCode = static_cast< int >( exitCode );
    else
    {
        result.exitCode = -1;
        result.errorMessage = "GetExitCodeProcess failed with Windows error " +
            std::to_string( GetLastError() );
    }
    CloseHandle( process.hThread );
    CloseHandle( process.hProcess );
#else
    int outputPipe[ 2 ] = { -1, -1 };
    if ( pipe( outputPipe ) != 0 )
    {
        result.errorMessage = std::string( "pipe failed: " ) + std::strerror( errno );
        return result;
    }

    posix_spawn_file_actions_t actions;
    int spawnError = posix_spawn_file_actions_init( &actions );
    const bool actionsInitialized = spawnError == 0;
    if ( spawnError == 0 )
        spawnError = posix_spawn_file_actions_adddup2(
            &actions, outputPipe[ 1 ], STDOUT_FILENO );
    if ( spawnError == 0 )
        spawnError = posix_spawn_file_actions_adddup2(
            &actions, outputPipe[ 1 ], STDERR_FILENO );
    if ( spawnError == 0 )
        spawnError = posix_spawn_file_actions_addclose( &actions, outputPipe[ 0 ] );
    if ( spawnError == 0 )
        spawnError = posix_spawn_file_actions_addclose( &actions, outputPipe[ 1 ] );
#if defined(__GLIBC__) || defined(__APPLE__)
    if ( spawnError == 0 && ! workingDirectory.empty() )
        spawnError = posix_spawn_file_actions_addchdir_np(
            &actions, workingDirectory.c_str() );
#else
    if ( spawnError == 0 && ! workingDirectory.empty() )
        spawnError = ENOTSUP;
#endif

    std::vector< char * > argv;
    argv.reserve( arguments.size() + 2 );
    argv.push_back( const_cast< char * >( executable.c_str() ) );
    for ( const std::string & argument : arguments )
        argv.push_back( const_cast< char * >( argument.c_str() ) );
    argv.push_back( nullptr );

    pid_t processId = -1;
    if ( spawnError == 0 )
        spawnError = posix_spawnp(
            &processId, executable.c_str(), &actions, nullptr, argv.data(), environ );
    if ( actionsInitialized )
        posix_spawn_file_actions_destroy( &actions );
    close( outputPipe[ 1 ] );

    if ( spawnError != 0 )
    {
        result.errorMessage = std::string( "posix_spawnp failed: " ) +
            std::strerror( spawnError );
        close( outputPipe[ 0 ] );
        return result;
    }

    result.started = true;
    char buffer[ 4096 ];
    ssize_t bytesRead = 0;
    while ( ( bytesRead = read( outputPipe[ 0 ], buffer, sizeof( buffer ) ) ) > 0 )
        result.output.append( buffer, static_cast< std::size_t >( bytesRead ) );
    close( outputPipe[ 0 ] );

    int status = 0;
    while ( waitpid( processId, &status, 0 ) < 0 )
    {
        if ( errno == EINTR ) continue;
        result.errorMessage = std::string( "waitpid failed: " ) + std::strerror( errno );
        return result;
    }
    if ( WIFEXITED( status ) )
        result.exitCode = WEXITSTATUS( status );
    else if ( WIFSIGNALED( status ) )
        result.exitCode = 128 + WTERMSIG( status );
    else
        result.exitCode = -1;
#endif
    return result;
}

EndNameSpace
