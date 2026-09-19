#include "EulerDomainStateSync.h"

#include "CpuEulerDomainBackend.h"
#include "EulerDomainMrFieldAdapter.h"
#include "EulerDomainStateLifecycle.h"
#include "SimuContext.h"
#include "Ctrl.h"
#include "DataBase.h"
#include "FieldImp.h"
#include "Grid.h"
#include "GridState.h"
#include "NsCom.h"
#include "SolverState.h"
#include "Zone.h"
#include "ZoneState.h"

#include <cmath>
#include <memory>
#include <stdexcept>

BeginNameSpace( ONEFLOW )

namespace
{

bool BuildCurrentEulerDomainState(
    EulerDomainProblem & problem,
    EulerDomainStateKey & key,
    std::unique_ptr< EulerDomainMrFieldSnapshot > & snapshot )
{
    Grid * grid = Zone::GetGrid();
    if ( grid == nullptr || grid->nCells <= 0 ) return false;

    MRField * q = GetFieldPointer< MRField >( grid, "q" );
    if ( q == nullptr ) return false;

    const int nEquations = static_cast< int >( q->GetNEqu() );
    if ( nEquations != 3 && nEquations != 5 ) return false;

    Real minDistance = 0.0;
    Real maxDistance = 0.0;
    grid->GetMinMaxDistance( minDistance, maxDistance );
    if ( ! std::isfinite( minDistance ) || minDistance <= 0.0 )
    {
        throw std::runtime_error(
            "cannot derive a positive Euler domain length scale" );
    }
    if ( ! std::isfinite( ctrl.pdt ) || ctrl.pdt <= 0.0 )
    {
        throw std::runtime_error(
            "cannot initialize Euler domain state with invalid dt" );
    }
    if ( ! std::isfinite( nscom.gama_ref ) || nscom.gama_ref <= 1.0 )
    {
        throw std::runtime_error(
            "cannot initialize Euler domain state with invalid gamma" );
    }

    problem.nCells = grid->nCells;
    problem.nGhostCells = grid->nBFaces;
    problem.nEquations = nEquations;
    problem.gamma = nscom.gama_ref;
    problem.dt = ctrl.pdt;
    problem.dx = minDistance;

    key = {
        SolverState::solverIndex,
        ZoneState::zid,
        GridState::gridLevel,
        AccelBackendKind::CPU };
    snapshot = std::make_unique< EulerDomainMrFieldSnapshot >(
        *q, grid->nCells );
    return true;
}

}

void SyncCurrentEulerDomainState(
    SimuContext & context,
    CpuEulerDomainBackend & backend,
    bool restart )
{
    EulerDomainProblem problem;
    EulerDomainStateKey key;
    std::unique_ptr< EulerDomainMrFieldSnapshot > snapshot;
    if ( ! BuildCurrentEulerDomainState( problem, key, snapshot ) ) return;

    const EulerDomainConstFieldView field = snapshot->View();
    if ( restart )
    {
        context.RestartAccelState( backend, problem, key, field );
    }
    else
    {
        context.InitializeAccelState( backend, problem, key, field );
    }
}

void UploadCurrentEulerDomainState(
    SimuContext & context,
    CpuEulerDomainBackend & backend )
{
    EulerDomainProblem problem;
    EulerDomainStateKey key;
    std::unique_ptr< EulerDomainMrFieldSnapshot > snapshot;
    if ( ! BuildCurrentEulerDomainState( problem, key, snapshot ) ) return;

    const EulerDomainConstFieldView field = snapshot->View();
    if ( context.AccelStates().Contains( key ) )
    {
        backend.Upload( context.AccelStates().Get( key ), field );
    }
    else
    {
        context.InitializeAccelState( backend, problem, key, field );
    }
}

void SyncAllEulerDomainStates( SimuContext & context )
{
    CpuEulerDomainBackend backend;
    const bool restart = ctrl.startStrategy == 1 || ctrl.startStrategy == 3;
    for ( int solverIndex = 0;
        solverIndex < SolverState::nSolver; ++ solverIndex )
    {
        SolverState::SetSolverTypeBySolverIndex( solverIndex );
        for ( int gridLevel = 0;
            gridLevel < GridState::nGrids; ++ gridLevel )
        {
            GridState::SetGridLevel( gridLevel );
            SyncCurrentEulerDomainState( context, backend, restart );
        }
    }
}

EndNameSpace
