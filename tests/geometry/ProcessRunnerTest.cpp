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
#include "ProcessRunner.h"

#include <gtest/gtest.h>

#include <string>
#include <vector>

namespace
{
    TEST( ProcessRunnerTest, CapturesOutputAndPreservesArgumentBoundaries )
    {
#ifdef _WIN32
        const ONEFLOW::ProcessResult result = ONEFLOW::RunProcess(
            "cmd.exe", { "/d", "/c", "echo ProcessRunner argument with spaces" } );
#else
        const ONEFLOW::ProcessResult result = ONEFLOW::RunProcess(
            "/bin/echo", { "ProcessRunner argument with spaces" } );
#endif
        ASSERT_TRUE( result.started ) << result.errorMessage;
        EXPECT_EQ( result.exitCode, 0 ) << result.errorMessage << "\\n" << result.output;
        EXPECT_NE( result.output.find( "ProcessRunner argument with spaces" ),
                   std::string::npos );
    }

    TEST( ProcessRunnerTest, ReturnsNonzeroExitCodeAndCapturesDiagnostics )
    {
#ifdef _WIN32
        const ONEFLOW::ProcessResult result = ONEFLOW::RunProcess(
            "cmd.exe", { "/d", "/c", "echo ProcessRunner diagnostic 1>&2 & exit /b 7" } );
#else
        const ONEFLOW::ProcessResult result = ONEFLOW::RunProcess(
            "/bin/sh", { "-c", "printf 'ProcessRunner diagnostic\\n' >&2; exit 7" } );
#endif
        ASSERT_TRUE( result.started ) << result.errorMessage;
        EXPECT_EQ( result.exitCode, 7 ) << result.errorMessage << "\\n" << result.output;
        EXPECT_NE( result.output.find( "ProcessRunner diagnostic" ),
                   std::string::npos ) << result.errorMessage << "\\n" << result.output;
    }

    TEST( ProcessRunnerTest, ReportsMissingExecutableWithoutInvokingAShell )
    {
        const ONEFLOW::ProcessResult result = ONEFLOW::RunProcess(
            "oneflow_process_runner_executable_that_does_not_exist", {} );
        EXPECT_FALSE( result.started );
        EXPECT_FALSE( result.errorMessage.empty() );
    }
}
