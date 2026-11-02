/*---------------------------------------------------------------------------*\\
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
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.
----------------------------------------------------------------------------*/

#pragma once

#include <exception>
#include <utility>

namespace ONEFLOW
{
namespace SimuCaseLifecycleDetail
{

template< typename Execute, typename Cleanup >
void ExecuteWithCleanup( Execute && execute, Cleanup && cleanup )
{
    try
    {
        std::forward< Execute >( execute )();
    }
    catch ( ... )
    {
        try
        {
            std::forward< Cleanup >( cleanup )();
        }
        catch ( ... )
        {
            // Keep the original execution failure.
        }
        throw;
    }

    std::forward< Cleanup >( cleanup )();
}

// Run every independent teardown step, then propagate the first failure.
// This prevents one failed release from skipping later resource cleanup.
template< typename... Cleanup >
void ExecuteAllCleanupSteps( Cleanup &&... cleanup )
{
    std::exception_ptr firstFailure;
    auto runCleanup = [&firstFailure]( auto && step )
    {
        try
        {
            std::forward< decltype( step ) >( step )();
        }
        catch ( ... )
        {
            if ( ! firstFailure )
            {
                firstFailure = std::current_exception();
            }
        }
    };

    ( runCleanup( std::forward< Cleanup >( cleanup ) ), ... );

    if ( firstFailure )
    {
        std::rethrow_exception( firstFailure );
    }
}

template< typename Setup, typename Rollback >
void ExecuteWithRollback( Setup && setup, Rollback && rollback )
{
    try
    {
        std::forward< Setup >( setup )();
    }
    catch ( ... )
    {
        try
        {
            std::forward< Rollback >( rollback )();
        }
        catch ( ... )
        {
            // Keep the original setup failure.
        }
        throw;
    }
}

template< typename Setup, typename Execute, typename Finalize >
void ExecuteWithProcessLifecycle(
    Setup && setup,
    Execute && execute,
    Finalize && finalize )
{
    ExecuteWithCleanup(
        [&]()
        {
            std::forward< Setup >( setup )();
            std::forward< Execute >( execute )();
        },
        [&]() { std::forward< Finalize >( finalize )(); } );
}

template< typename CaseDirs, typename RunCase >
void RunCasesSequentially( const CaseDirs & caseDirs, RunCase && runCase )
{
    for ( const auto & caseDir : caseDirs )
    {
        runCase( caseDir );
    }
}

template< typename Setup, typename CaseDirs, typename RunCase, typename Finalize >
void RunCasesWithProcessLifecycle(
    Setup && setup,
    const CaseDirs & caseDirs,
    RunCase && runCase,
    Finalize && finalize )
{
    ExecuteWithProcessLifecycle(
        std::forward< Setup >( setup ),
        [&]()
        {
            RunCasesSequentially(
                caseDirs, std::forward< RunCase >( runCase ) );
        },
        std::forward< Finalize >( finalize ) );
}

}
}