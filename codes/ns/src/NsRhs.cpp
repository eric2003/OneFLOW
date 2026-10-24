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

#include "NsRhs.h"
#include "UNsBcSolver.h"
#include "Zone.h"
#include "DataBase.h"
#include "Iteration.h"
#include "NsCom.h"
#include "UCom.h"
#include "UNsCom.h"
#include "UnsGrid.h"
#include "NsIdx.h"
#include "UNsInvFlux.h"
#include "UNsVisFlux.h"
#include "UNsUnsteady.h"
#include "Ctrl.h"
#include "FieldWrap.h"
#include "HipProductionPolicy.h"
#ifdef ONEFLOW_ENABLE_HIP
#include "HipFluxBackend.h"
#endif

#include <iostream>
#include <memory>


BeginNameSpace( ONEFLOW )

void NsRhs::UpdateResiduals()
{
	NsCalcRHS();
}

void NsCalcBc()
{
	auto uNsBcSolver = std::make_unique<UNsBcSolver>();
	uNsBcSolver->Init();
	uNsBcSolver->CalcBc();
}

void NsCalcBcDebug( const std::string & title )
{
	std::cout << title << "\n";
	auto uNsBcSolver = std::make_unique<UNsBcSolver>();
	uNsBcSolver->Init();
	uNsBcSolver->CalcBc();

}

void NsCalcGamaT( int flag )
{
	ug.Init();
	unsf.Init();
	ug.SetStEd(flag);

	if ( nscom.chemModel == 1 )
	{
	}
	else
	{
		Real oamw = one;
		Real * density = ( * unsf.q )[ IDX::IR ].data();
		Real * pressure = ( * unsf.q )[ IDX::IP ].data();
		Real * gama = ( * unsf.gama )[ 0 ].data();
		// Inviscid three-equation configurations set nTModel=0 and
		// consequently allocate no temperature components.  Gamma is still
		// required by the flux path, but there is no temperature field to
		// populate in that case.
		Real * temperature = nullptr;
		if ( nscom.nTModel > IDX::ITT && unsf.tempr != nullptr
		     && unsf.tempr->GetNEqu() > IDX::ITT )
		{
			temperature = ( * unsf.tempr )[ IDX::ITT ].data();
		}
		for ( int cId = ug.ist; cId < ug.ied; ++ cId )
		{
			gama[ cId ] = nscom.gama_ref;
			if ( temperature != nullptr )
			{
				temperature[ cId ] =
					pressure[ cId ] /
					( nscom.statecoef * density[ cId ] * oamw );
			}
		}
	}
}

void NsRefreshHostStateMirror()
{
#ifdef ONEFLOW_ENABLE_HIP
    if ( ! HipProductionBackendSelected() ) return;
    UnsGrid * grid = Zone::GetUnsGrid();
    if ( grid == nullptr ) return;
    MRField * q = GetFieldPointer< MRField >( grid, "q" );
    if ( q == nullptr ) return;
    const int nEquations = static_cast< int >( q->GetNEqu() );
    if ( nEquations <= 0 || nEquations > 5 ) return;
    ug.Init();
    const int nCells = ug.nCells;
    const int nGhostCells = ug.nTCell - nCells;
    Real * components[ 5 ] = {};
    for ( int equation = 0; equation < nEquations; ++ equation )
        components[ equation ] = ( * q )[ equation ].data();
    if ( ! HipFluxBackend::Shared().RefreshHostPrimitiveState(
             grid, nCells, nGhostCells, nEquations, components ) )
    {
        return;
    }
    // The host-only derived fields (constant gamma, temperature) were skipped
    // while the device owned q; rebuild them from the refreshed mirror so the
    // restart/visualization output is self-consistent.
    NsCalcGamaT( F_INNER );
    if ( nGhostCells > 0 ) NsCalcGamaT( F_GHOST );
#else
    return;
#endif
}

void NsMirrorHostBoundaryState()
{
#ifdef ONEFLOW_ENABLE_HIP
    if ( ! HipDeviceAuthoritativeStateEnabled() ) return;
    if ( ! HipBoundaryStateMirrorEnabled() ) return;
    UnsGrid * grid = Zone::GetUnsGrid();
    if ( grid == nullptr || grid->nBFaces <= 0 ) return;
    MRField * q = GetFieldPointer< MRField >( grid, "q" );
    if ( q == nullptr ) return;
    const int nEquations = static_cast< int >( q->GetNEqu() );
    if ( nEquations <= 0 || nEquations > 5 ) return;
    ug.Init();
    if ( ug.lcf == nullptr ) return;
    Real * components[ 5 ] = {};
    for ( int equation = 0; equation < nEquations; ++ equation )
        components[ equation ] = ( * q )[ equation ].data();
    // Every boundary face's owner cell is a superset of the wall-adjacent cells
    // the force/coefficient readers use.  Duplicate ids are harmless (the
    // same value is written twice) and avoid a sort on the per-step path.
    HipFluxBackend::Shared().MirrorPrimitiveCellsToHost(
        grid, nEquations, ug.lcf->data(), grid->nBFaces, components );
#else
    return;
#endif
}

void NsCalcRHS()
{
	NsCalcInvFlux();

	NsCalcVisFlux();

	NsCalcSrcFlux();
}


void NsCalcInvFlux()
{
    // Reuse the host-side limiter/face work buffers across RK stages.  The
    // fields are rebound in Init() for the active zone/grid, while qf1/qf2
    // and invflux retain their capacity when the topology is unchanged.
    static UNsInvFlux uNsInvFlux;
    uNsInvFlux.CalcFlux();
}

void NsCalcVisFlux()
{
	auto uNsVisFlux = std::make_unique<UNsVisFlux>();
	uNsVisFlux->CalcFlux();
}

void NsCalcSrcFlux()
{
	if (nscom.chemModel == 1)
	{
		NsCalcChemSrc();
	}

	NsCalcTurbEnergy();

	//dual time step source
	if ( ctrl.idualtime == 1 )
	{
		NsCalcDualTimeStepSrc();
	}
}

void NsCalcChemSrc()
{
	;
}

void NsCalcTurbEnergy()
{
	;
}

void NsCalcDualTimeStepSrc()
{
	auto unsUnsteady = std::make_unique<UNsUnsteady>();
	unsUnsteady->CalcDualTimeSrc();
}

EndNameSpace
