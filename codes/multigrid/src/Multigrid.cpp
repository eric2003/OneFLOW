/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2024 He Xin and the OneFLOW contributors.
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

#include "Multigrid.h"
#include "INsInvterm.h"
#include "Mesh.h"
#include "Ctrl.h"
#include "Solver.h"
#include "Iteration.h"
#include "DataBase.h"
#include "UNsSolver.h"
#include "Zone.h"
#include "TimeIntegral.h"
#include "Parallel.h"
#include "SolverState.h"
#include "GridState.h"
#include "CmxTask.h"
#include "BgField.h"
#include "TimeSpan.h"
#include <iostream>



BeginNameSpace( ONEFLOW )

int MG::nPre;
int MG::nPost;
int MG::nMulti;
int MG::iterMode;
int MG::cycleType;
int MG::currentGridLevel = 0;

MG::MG()
{
    ;
}

MG::~MG()
{
    ;
}

void MG::Init()
{
    MG::nPre = ONEFLOW::GetDataValue< int >( "nmgPreSmoothingSteps" );
    MG::nPost = ONEFLOW::GetDataValue< int >( "nmgPostSmoothingSteps" );
    MG::nMulti = ONEFLOW::GetDataValue< int >( "nmg" );
    MG::iterMode = ONEFLOW::GetDataValue< int >( "mgIterMode" );
    MG::cycleType = ONEFLOW::GetDataValue< int >( "mgCycleType" );

    if ( MG::nMulti <= 1 )
    {
        MG::nPre = 1;
        MG::nPost = 0;
    }
}

void MG::Allocate()
{
    MG::Init();
    BgField::Init();
}

void MG::Deallocate()
{
    BgField::Free();
}

void MG::MWrap( FunctionPointer multigridPointer, int gridLevel )
{
    if ( this->iterMode == 0 )
    {
        ( this->*multigridPointer )( gridLevel );
    }
    else
    {
        for ( int sid = 0; sid < SolverState::nSolver; ++ sid )
        {
            SolverState::SetTidById( sid );
            ( this->*multigridPointer )( gridLevel );
        }
    }
}

void MG::MultigridSolve()
{
    this->Allocate();
    this->Run();
    this->Deallocate();
}

void MG::Run()
{
	int startStrategy = ONEFLOW::GetDataValue< int >("startStrategy");
	if (startStrategy == 2|| startStrategy == 3)
	{
		double rhs_V = 1e-8;
		double rhs_u = 1e-8;
		double rhs_v = 1e-8;
		double rhs_w = 1e-8;

		int maxIterSteps = GetDataValue< int >("maxIterSteps");

		iinv.remax_V = 1;
		iinv.remax_up = 1;
		iinv.remax_vp = 1;
		iinv.remax_wp = 1;

		TimeSpan * timeSpan = new TimeSpan();
		while (SimuIterState::Running())
		{

			while (iinv.remax_up > rhs_u || iinv.remax_vp > rhs_v || iinv.remax_wp > rhs_w)
			{

				if (Iteration::innerSteps >= maxIterSteps) break;

				ctrl.currTime += ctrl.pdt;

				Iteration::outerSteps++;
				Iteration::innerSteps++;

				this->SolveInnerIter();

			}
			this->OuterProcess(timeSpan);
		}
		delete timeSpan;
	}
	else
	{
		TimeSpan * timeSpan = new TimeSpan();
		while ( SimuIterState::Running() )
		{
			Iteration::outerSteps ++;
			ctrl.currTime += ctrl.pdt;

			//Inner loop
			Iteration::innerSteps = 0;
			while ( !SolverState::Converge() )
			{
				Iteration::innerSteps ++;

				this->SolveInnerIter();
			}
			this->OuterProcess( timeSpan );
		}
		delete timeSpan;
	}
}

void MG::InnerProcess()
{
    ONEFLOW::MsMgTask( "POST_PROCESS" );
}

void MG::OuterProcess( TimeSpan * timeSpan )
{
    if ( Iteration::outerSteps % Iteration::nFieldSave == 0 )
    {
        if ( Parallel::IsServer() )
        {
            std::cout << "dumping field...";
            std::cout << "  finished " << std::endl;
            timeSpan->ShowTimeSpan();
        }
    }
}

void MG::PreprocessMultigridFlowField( int gl )
{
    GridState::SetGridLevel( gl );

    ONEFLOW::SsSgTask( "STORE_RHS" );
}

void MG::InitializeCoarseGridFlowFieldByRestrictFineGridFlowField( int fgl )
{
    //+ Restrict from fine to coarse grid for all q
    //+ Corresponding nsmb: wsav = restr (W)
    //Call the following statement to update the Q value (CQ) on the thin grid
    GridState::SetGridLevel( fgl );

    ONEFLOW::SsSgTask( "RESTRICT_ALL_Q" );
}

void MG::StoreCoarseGridFlowFieldToTemporaryStorage( int fgl )
{
    int cgl = GridState::GetCGridLevel( fgl );
    //+ After loadq (coat grid, cqsav), cqsav is actually wsav (corresponding to nsmb notation)
    //Loadq takes the Q value (CQ) from the sparse grid and assigns it to cqsav
    GridState::SetGridLevel( cgl );

    ONEFLOW::SsSgTask( "LOAD_Q" );
}

void MG::PrepareFineGridResiduals( int fgl )
{
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+ Where (* residual) = - (* RHS) corresponds to - F
    //+ Note that residual is the value on the dense grid, and RHS should be f. check this
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //This should be Ni, NJ, NK instead of CNI, CNJ, CNK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    //Initialize residual to the last driver (F is 0 for the first time, and there will be a value later)
    //Residual = ( - f ) == ( - rhs );
    GridState::SetGridLevel( fgl );

    ONEFLOW::SsSgTask( "LOAD_RESIDUALS" );

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+After updateresiduals, generalresidualfield = Rl (W) - F
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ONEFLOW::SsSgTask( "UPDATE_RESIDUALS" );
}

void MG::PrepareCoarseGridResiduals( int fgl )
{
    //+Residual incoarsegrid = - restr (RL (W) - F) after restrict defect
    GridState::SetGridLevel( fgl );
    ONEFLOW::SsSgTask( "RESTRICT_DEFECT" );

    //After updateresiduals, residualincoarsegrid = RL-1 (wsav) - restr (RL (W) - F)
    int cgl = GridState::GetCGridLevel( fgl );
    GridState::SetGridLevel( cgl );
    ONEFLOW::SsSgTask( "UPDATE_RESIDUALS" );
}

void MG::SolveCoarseGridFlowField( int fgl )
{
    int cgl = GridState::GetCGridLevel( fgl );
    //After solvecoarsegridflowfield (CGL), the Q value on the thin grid will be updated
    //In this case, the Q value on the sparse grid is equal to W0 in nsmb

    //cycleType = 1 V cycle
    //cycleType = 2 W cycle
    for ( int iCycle = 0; iCycle < MG::cycleType; ++ iCycle )
    {
        this->SolveMultigridFlowField( cgl );
    }
}

void MG::CorrectFineGridFlowFieldByInterplateCoarseGridFlowField( int fgl )
{
    int cgl = GridState::GetCGridLevel( fgl );

    //w0 - wsav
    GridState::SetGridLevel( cgl );

    ONEFLOW::SsSgTask( "MODIFY_COARSEGRID" );

    //w = w + prol( w0 - wsav )
    GridState::SetGridLevel( fgl );

    ONEFLOW::SsSgTask( "MODIFY_FINEGRID" );

    //The following is actually to restore the Q value on the sparse grid.
    GridState::SetGridLevel( cgl );

    ONEFLOW::SsSgTask( "RECOVER_COARSEGRID" );
}

void MG::PreRelaxationCycle( int gl )
{
    GridState::SetGridLevel( gl );

    TimeIntegral::Relaxation( MG::nPre );
}

void MG::PostRelaxationCycle( int gl )
{
    GridState::SetGridLevel( gl );

    TimeIntegral::Relaxation( MG::nPost );
}

void MG::PostprocessMultigridFlowField( int gl )
{
    // This is only to restore the initial value of the general residual field. As for the usefulness of this, let's say something else.
    GridState::SetGridLevel( gl );

    ONEFLOW::SsSgTask( "RECOVER_RESIDUALS" );
}

void MG::FastSolveFlowFieldByMultigridMethod( int gl )
{
    // If grid is in the coarest gridLevelIndex
    if ( ONEFLOW::DoNotNeedMultigridMethod( gl )  ) return;

    this->MWrap( & MG::PreprocessMultigridFlowField                            , gl );
    this->MWrap( & MG::InitializeCoarseGridFlowFieldByRestrictFineGridFlowField, gl );
    this->MWrap( & MG::StoreCoarseGridFlowFieldToTemporaryStorage              , gl );
    this->MWrap( & MG::PrepareFineGridResiduals                                , gl );
    this->MWrap( & MG::PrepareCoarseGridResiduals                              , gl );
    this->SolveCoarseGridFlowField( gl );
    this->MWrap( & MG::CorrectFineGridFlowFieldByInterplateCoarseGridFlowField , gl );
}

void MG::SolveMultigridFlowField( int gl )
{
	int startStrategy = ONEFLOW::GetDataValue< int >("startStrategy");
	if (startStrategy == 2|| startStrategy == 3)
	{
		//this->MWrap(&MG::PreprocessMultigridFlowField, gl);
		this->MWrap(&MG::PreRelaxationCycle, gl);
		//this->FastSolveFlowFieldByMultigridMethod(gl);
		//this->MWrap(&MG::PostRelaxationCycle, gl);
		//this->MWrap(&MG::PostprocessMultigridFlowField, gl);
	}
	else
	{
		this->MWrap(&MG::PreprocessMultigridFlowField, gl);
		this->MWrap(&MG::PreRelaxationCycle, gl);
		this->FastSolveFlowFieldByMultigridMethod(gl);
		this->MWrap(&MG::PostRelaxationCycle, gl);
		this->MWrap(&MG::PostprocessMultigridFlowField, gl);
	}
}

void MG::SolveInnerIter()
{
    if ( MG::iterMode == 0 )
    {
        this->WeakIter();
    }
    else
    {
        this->StrongIter();
    }

    this->InnerProcess();
}

void MG::ZeroResidualsForAllSolvers()
{
    for ( int sId = 0; sId < SolverState::nSolver; ++ sId )
    {
        SolverState::SetTidById( sId );
        ONEFLOW::SsSgTask( "ZERO_RESIDUALS" );
    }
}

void MG::WeakIter()
{
    this->ZeroResidualsForAllSolvers();

    for ( int sid = 0; sid < SolverState::nSolver; ++ sid )
    {
        SolverState::SetTidById( sid );
        this->SolveMultigridFlowField( 0 );
    }
}

void MG::StrongIter()
{
	int startStrategy = ONEFLOW::GetDataValue< int >("startStrategy");
	if (startStrategy == 2|| startStrategy == 3)
	{
		//this->ZeroResidualsForAllSolvers();
		this->SolveMultigridFlowField(0);
	}
	else
	{
		this->ZeroResidualsForAllSolvers();

		this->SolveMultigridFlowField(0);
	}
}


bool DoNotNeedMultigridMethod( int gl )
{
    return true;
}

void MultigridSolve()
{
    MG * mg = new MG();
    mg->MultigridSolve();
    delete mg;
}

EndNameSpace
