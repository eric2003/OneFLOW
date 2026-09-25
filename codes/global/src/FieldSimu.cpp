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
#include "FieldSimu.h"
#include "SimuContext.h"
#include "EulerDomainStateSync.h"
#include "CmxTaskNames.h"
#include "Iteration.h"
#include "Ctrl.h"
#include "NsCom.h"
#include "UsdData.h"
#include "MultiBlock.h"
#include "SolverMap.h"
#include "SolverCatalog.h"
#include "CmxTask.h"
#include "Multigrid.h"
#include "BcData.h"
#include "GridState.h"
#include "FieldManager.h"
#include "SolverState.h"
#include "SolverDef.h"
#include "Parallel.h"
#include "RegisterUtils.h"
#include <iostream>
#include <stdexcept>

BeginNameSpace( ONEFLOW )

void FieldSimuSetupGlobals()
{
    InitFlowSimuGlobal();
}

void FieldSimuLoadGrid()
{
    MultiBlock::LoadGridAndBuildLink();
}

void FieldSimuPrepareWallDist()
{
    MultiBlock::ProcessFlowWallDist();
}


void FieldSimuCreateSolvers()
{
    // Prefer SolverCatalog name at the pipeline boundary (owns via SolverMap today).
    SolverCatalog::CreateDefault();
}

void FieldSimuCreateSolvers( const SimuContext & ctx )
{
    if ( ctx.HasExpandedSolverNames() )
    {
        SolverCatalog::CreateDefault(
            ONEFLOW::UMESH,
            &ctx.ExpandedSolverNames() );
    }
    else
    {
        SolverCatalog::CreateDefault();
    }
}

void FieldSimuInitFlowField()
{
    // Stage entry: task name enters CmxTask here (no numerical change)
    ONEFLOW::MultiSolverMultiGridTask( kInitFlowFieldTaskName );
}

void DumpFieldEnvironments()
{
    if ( Parallel::GetPid() != Parallel::GetServerid() )
    {
        return;
    }

    const int savedSolverIndex =
        SolverState::solverIndex;

    const int savedSolverType =
        SolverState::solverType;

    std::cout
        << "\n"
        << "========================================\n"
        << "       Field Environment Summary\n"
        << "========================================\n";

    for ( int solverIndex = 0;
        solverIndex < SolverState::nSolver;
        ++ solverIndex )
    {
        SolverState::SetSolverTypeBySolverIndex(
            solverIndex );

        const int solverType =
            SolverState::solverType;

        FieldManager * fieldManager =
            FieldManagerRegistry::GetFieldManager(
                solverType );

        std::cout
            << "\n"
            << "[Solver "
            << solverIndex
            << ", type "
            << solverType
            << "]\n";

        if ( fieldManager == nullptr )
        {
            std::cout
                << "  <FieldManager not found>\n";
            continue;
        }

        fieldManager->DumpFieldEnvironment(
            std::cout );
    }

    SolverState::solverIndex =
        savedSolverIndex;

    SolverState::solverType =
        savedSolverType;

    //std::cout
    //    << "========================================\n"
    //    << std::endl;
}

void DumpCommunicationEnvironments()
{
    if ( Parallel::GetPid() != Parallel::GetServerid() )
    {
        return;
    }

    const int savedSolverIndex =
        SolverState::solverIndex;

    const int savedSolverType =
        SolverState::solverType;

    std::cout
        << "\n"
        << "========================================\n"
        << "    Communication Environment Summary\n"
        << "========================================\n";

    for ( int solverIndex = 0;
        solverIndex < SolverState::nSolver;
        ++ solverIndex )
    {
        SolverState::SetSolverTypeBySolverIndex(
            solverIndex );

        const int solverType =
            SolverState::solverType;

        std::cout
            << "\n"
            << "[Solver "
            << solverIndex
            << ", type "
            << solverType
            << "]\n";

        VarNameFactory::Dump(
            std::cout,
            solverType );
    }

    SolverState::solverIndex =
        savedSolverIndex;

    SolverState::solverType =
        savedSolverType;
}

namespace
{
    bool CheckInterfaceStorageContainsCommunicationField(
        std::ostream & output,
        const InterfaceFieldProperty & interfaceFieldProperty,
        VarNameSolver * varNameSolver )
    {
        if ( varNameSolver == nullptr )
        {
            return true;
        }

        const FieldDefinitionTable::Data & interfaceFields =
            interfaceFieldProperty.GetData();

        bool consistent = true;

        for ( int fieldId = 0;
            fieldId < varNameSolver->data.size();
            ++ fieldId )
        {
            const std::string & fieldName =
                varNameSolver->data[ fieldId ];

            if ( interfaceFields.find( fieldName ) ==
                interfaceFields.end() )
            {
                consistent = false;

                output
                    << "    ERROR: communication field \""
                    << fieldName
                    << "\" is not registered in Interface Storage\n";
            }
        }

        return consistent;
    }
}

void CheckCommunicationInterfaceConsistency()
{
    if ( Parallel::GetPid() != Parallel::GetServerid() )
    {
        return;
    }

    const int savedSolverIndex =
        SolverState::solverIndex;

    const int savedSolverType =
        SolverState::solverType;

    std::cout
        << "\n"
        << "========================================\n"
        << " Communication / Interface Consistency\n"
        << "========================================\n";

    const int interfaceTypes[] =
    {
        ONEFLOW::INTERFACE_DATA,
        ONEFLOW::INTERFACE_DQ_DATA,
        ONEFLOW::INTERFACE_GRADIENT_DATA,
        ONEFLOW::INTERFACE_OVERSET_DATA
    };

    const char * interfaceNames[] =
    {
        "INTERFACE_DATA",
        "INTERFACE_DQ",
        "INTERFACE_GRADIENT",
        "INTERFACE_OVERSET"
    };

    for ( int solverIndex = 0;
        solverIndex < SolverState::nSolver;
        ++ solverIndex )
    {
        SolverState::SetSolverTypeBySolverIndex(
            solverIndex );

        const int solverType =
            SolverState::solverType;

        FieldManager * fieldManager =
            FieldManagerRegistry::GetFieldManager(
                solverType );

        std::cout
            << "\n"
            << "[Solver "
            << solverIndex
            << ", type "
            << solverType
            << "]\n";

        if ( fieldManager == nullptr )
        {
            std::cout
                << "  <FieldManager not found>\n";
            continue;
        }

        const InterfaceFieldProperty & interfaceFieldProperty =
            fieldManager->GetInterfaceFieldProperty();

        bool consistent = true;

        for ( int iType = 0; iType < 4; ++ iType )
        {
            VarNameSolver * varNameSolver =
                VarNameFactory::GetVarNameSolver(
                    solverType,
                    interfaceTypes[ iType ] );

            bool groupConsistent =
                CheckInterfaceStorageContainsCommunicationField(
                    std::cout,
                    interfaceFieldProperty,
                    varNameSolver );

            if ( ! groupConsistent )
            {
                consistent = false;

                std::cout
                    << "  "
                    << interfaceNames[ iType ]
                    << ": INCONSISTENT\n";
            }
        }

        if ( consistent )
        {
            std::cout
                << "  Result: OK\n";
        }
        else
        {
            std::cout
                << "  Result: INCONSISTENT\n";
        }
    }

    SolverState::solverIndex =
        savedSolverIndex;

    SolverState::solverType =
        savedSolverType;
}

void FieldSimuRun()
{
    MultigridSolve();
}

void FieldSimuRun( SimuContext & context )
{
    MultigridSolve( context );
}

void FieldPipeline::Run()
{
    FieldSimuSetupGlobals();
    FieldSimuLoadGrid();
    FieldSimuPrepareWallDist();
    FieldSimuCreateSolvers();
    FieldSimuInitFlowField();
    //DumpFieldEnvironments();
    //DumpCommunicationEnvironments();
    //CheckCommunicationInterfaceConsistency();
    FieldSimuRun();
}

void FieldPipeline::Run( SimuContext & ctx )
{
    FieldSimuSetupGlobals();
    FieldSimuLoadGrid();
    FieldSimuPrepareWallDist();
    FieldSimuCreateSolvers( ctx );
    FieldSimuInitFlowField();
    //DumpFieldEnvironments();
    //DumpCommunicationEnvironments();
    //CheckCommunicationInterfaceConsistency();
    SyncAllEulerDomainStates( ctx );
    FieldSimuRun( ctx );
}

void FieldSimuRunPipeline()
{
    FieldPipeline::Run();
}

void FieldSimuRunPipeline( SimuContext & ctx )
{
    FieldPipeline::Run( ctx );
}

void FieldSimu()
{
    FieldPipeline::Run();
}

void FieldSimu( SimuContext & context )
{
    FieldSimuRunPipeline( context );
}

void InitFlowSimuGlobal()
{
    vis_model.Init();
    ctrl.Init();
    Iteration::Init();
    usd.InitBasic();
}

void InitializeSolver()
{
    // Compatibility alias for older call sites
    FieldSimuInitFlowField();
}

EndNameSpace
