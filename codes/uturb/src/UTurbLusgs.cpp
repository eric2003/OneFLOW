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

#include "UTurbLusgs.h"
#include "UNsSpectrum.h"
#include "UTurbSpectrum.h"
#include "UTurbCom.h"
#include "UCom.h"
#include "TurbCom.h"
#include "TurbRhs.h"
#include "UNsCom.h"
#include "Com.h"
#include "UnsGrid.h"
#include "CellMesh.h"
#include "CellTopo.h"
#include "Zone.h"
#include "Ctrl.h"
#include "NsIdx.h"
#include "HXMath.h"
#include "Parallel.h"
#include "Iteration.h"
#include "HipProductionPolicy.h"
#ifdef ONEFLOW_ENABLE_HIP
#include "HipTurbulenceBackend.h"
#endif
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <vector>
#include <iostream>


BeginNameSpace( ONEFLOW )

UTurbLusgs::UTurbLusgs()
{
}

UTurbLusgs::~UTurbLusgs()
{
}

void UTurbLusgs::SingleSweep()
{
    this->LowerSweep();
    this->UpperSweep();
}

void UTurbLusgs::Init()
{
    UnsGrid * grid = Zone::GetUnsGrid();
    FaceTopo & faceTopo = grid->GetFaceTopo();
    CellTopo & cellTopo = grid->GetCellMesh().GetCellTopo();
    cellTopo.CalcC2f( faceTopo );
    ug.Init();
    turblu.Init();
    uturbf.Init();
}

void UTurbLusgs::LowerSweep()
{
    if ( this->SweepHip( true ) ) return;

    this->Init();

    for ( int cId = 0; cId < ug.nCells; ++ cId )
    {
        ug.cId = cId;

        gcom.blank = ( * ug.blankf )[ ug.cId ];

        if ( this->IsOversetCell() )
        {
            this->ZeroOversetCell();
        }
        else
        {
            this->PrepareSweep();

            this->ZeroFluxIncrement();
        
            this->SolveLowerCell();

            this->CalcLowerChange();
        }

        this->Update();
    }
    
    //UploadInterfaceValue( grid, dqField, "dqField",  numberOfTotalEquations );
}

void UTurbLusgs::UpperSweep()
{
    if ( this->SweepHip( false ) ) return;

    this->Init();
    //this->LusgsBoundary();
    //DownloadInterfaceValue( grid, dqField, "dqField",  numberOfTotalEquations );
    
    for ( int cId = ug.nCells - 1; cId >= 0; -- cId )
    {
        ug.cId = cId;
        gcom.blank = ( * ug.blankf )[ ug.cId ];

        if ( this->IsOversetCell() )
        {
            this->ZeroOversetCell();
        }
        else
        {
            this->PrepareSweep();

            this->ZeroFluxIncrement();

            this->SolveUpperCell();

            this->CalcUpperChange();
        }

        this->Update();
    }

    //CalcTurbulentViscosity();
}

void UTurbLusgs::SolveLowerCell()
{
    int fn = ( * ug.c2f )[ ug.cId ].size();
    for ( int iFace = 0; iFace < fn; ++ iFace )
    {
        int fId = ( * ug.c2f )[ ug.cId ][ iFace ];
        this->SolveLower( fId );
    }
}

void UTurbLusgs::SolveUpperCell()
{
    int fn = ( * ug.c2f )[ ug.cId ].size();
    for ( int iFace = 0; iFace < fn; ++ iFace )
    {
        int fId = ( * ug.c2f )[ ug.cId ][ iFace ];
        this->SolveUpper( fId );
    }
}

void UTurbLusgs::SolveLower( int fId )
{
    if ( this->CanNotLowerSolve( fId ) ) return;

    this->Solve( fId, - 1 );
}

void UTurbLusgs::SolveUpper( int fId )
{
    if ( this->CanNotUpperSolve( fId ) ) return;

    this->Solve( fId, - 1 );
}

bool UTurbLusgs::CanNotLowerSolve( int fId )
{
    ug.fId = fId;

    ug.lc = ( * ug.lcf )[ ug.fId ];
    ug.rc = ( * ug.rcf )[ ug.fId ];

    // One of lc  and rc must be cell itself.
    // Now its neighboring cell belongs to lower triangular
    return ( ug.lc > ug.cId || ug.rc > ug.cId );
}

bool UTurbLusgs::CanNotUpperSolve( int fId )
{
    ug.fId = fId;

    ug.lc = ( * ug.lcf )[ ug.fId ];
    ug.rc = ( * ug.rcf )[ ug.fId ];

    // One of lc  and rc must be cell itself.
    // Now its neighboring cell belongs to upper triangular
    return ( ug.lc < ug.cId || ug.rc < ug.cId );
}

void UTurbLusgs::Solve( int fId, int signValue )
{
    ug.fId = fId;

    ug.lc = ( * ug.lcf )[ ug.fId ];
    ug.rc = ( * ug.rcf )[ ug.fId ];

    this->SetMeshGeometry();

    this->PrepareData();

    this->GetStandardFluxIncrement( signValue );

    this->AddFluxIncrement();
}

void UTurbLusgs::SetMeshGeometry()
{
    gcom.SetGeometry();
}

void UTurbLusgs::PrepareData()
{
    if ( gcom.swapflag )
    {
        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            turblu.matrix[ iEqu ] = ( * uturbf.matrix_l )[ iEqu ][ ug.fId ];
        }
    }
    else
    {
        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            turblu.matrix[ iEqu ] = ( * uturbf.matrix_r )[ iEqu ][ ug.fId ];
        }
    }

    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turblu.dqj[ iEqu ] = ( * uturbf.dq )[ iEqu ][ ug.rc ];
    }

}

void UTurbLusgs::PrepareSweep()
{
    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        gcom.blank = ( * ug.blankf )[ ug.cId ];

        turblu.dqi[ iEqu ] = ( * uturbf.dq  )[ iEqu ][ ug.cId ]; //The initial value of dqfield is 0 (conserved or original variable)
        turblu.rhs[ iEqu ] = ( * uturbf.rhs )[ iEqu ][ ug.cId ]; //It is better to have RHS in RHS
    }

    if ( turblu.numberOfSweeps > 1 )
    {
        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            turblu.dqi0[ iEqu ] = turblu.dqi[ iEqu ];
            turblu.drhs[ iEqu ] = ( * uturbf.drhs )[ iEqu ][ ug.cId ];

            turblu.dqi[ iEqu ] = ( * uturbf.rhs )[ iEqu ][ ug.cId ] - turblu.drhs[ iEqu ];
            ( * unsf.dq )[ iEqu ][ ug.cId ] = turblu.dqi[ iEqu ];
            turblu.drhs[ iEqu ] = 0.0;
        }
    }

    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turblu.radius[ iEqu ] = ( * uturbf.impsr )[ iEqu ][ ug.cId ];
    }
}

void UTurbLusgs::Update()
{
    if ( turblu.numberOfSweeps > 1 )
    {
        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            ( * uturbf.dq   )[ iEqu ][ ug.cId ]  = turblu.dqi[ iEqu ];
            ( * uturbf.drhs )[ iEqu ][ ug.cId ]  = turblu.drhs[ iEqu ];
        }
    }
    else
    {
        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            ( * uturbf.dq   )[ iEqu ][ ug.cId ]  = turblu.dqi[ iEqu ];
        }
    }
}



namespace
{

// Dependency levels for one sweep direction.
//
// The host sweep is Gauss-Seidel in cell-id order: cell `c` consumes the
// already-updated delta of every neighbour it depends on.  A cell's level is
// one more than the highest level among those dependencies, so all cells on a
// level are mutually independent and launching the levels in order reproduces
// the sequential chain exactly - not approximately.
//
// `lower` selects the direction: a lower sweep depends on neighbours with a
// smaller id, an upper sweep on neighbours with a larger id.  A dependency
// cycle (only possible with a non-monotone ordering, which the face-id-ordered
// host walk does not have) would leave a cell unassigned; the caller treats an
// unassigned cell as an error rather than silently dropping it.
void ComputeSweepLevels(
    const std::vector< int > & offsets,
    const std::vector< int > & ids,
    const int * leftCell,
    const int * rightCell,
    int nCells,
    int nBoundaryFaces,
    bool lower,
    std::vector< int > & levelOffset,
    std::vector< int > & levelCells )
{
    std::vector< int > level( static_cast< std::size_t >( nCells ), -1 );

    // Relaxation to a fixed point.  Each pass assigns every cell whose
    // dependencies are already assigned; the number of passes is the level
    // count, which for a structured 3D mesh is O(n^(1/3)) * 3 and therefore
    // small.  A cycle would stop making progress, which the final check
    // reports.
    int assigned = 0;
    int currentLevel = 0;
    while ( assigned < nCells )
    {
        int newlyAssigned = 0;
        for ( int cell = 0; cell < nCells; ++ cell )
        {
            if ( level[ cell ] >= 0 ) continue;

            bool ready = true;
            int dependencyLevel = -1;
            for ( int slot = offsets[ cell ]; slot < offsets[ cell + 1 ]; ++ slot )
            {
                const int face = ids[ slot ];
                const int left = leftCell[ face ];
                const bool isLeft = left == cell;
                const int neighbour = isLeft ? rightCell[ face ] : left;

                // The host's CanNotLowerSolve / CanNotUpperSolve filters a face
                // out of the triangle when the other cell is on the wrong side,
                // so only these faces are actual dependencies.
                if ( lower )
                {
                    if ( neighbour > cell ) continue;
                }
                else
                {
                    if ( neighbour < cell ) continue;
                }
                // A boundary face's ghost is not swept: the loop runs over
                // owned cells only, and the turbulence BC writes q[ghost], not
                // dq[ghost].  Its delta is therefore constant for the whole
                // sweep, so it contributes a fixed term and imposes no ordering
                // constraint.  The device kernel still reads it, which is why
                // the delta array spans the ghosts; the value there is whatever
                // the host array holds, copied verbatim.
                if ( neighbour >= nCells ) continue;  // constant, not swept
                if ( level[ neighbour ] < 0 )
                {
                    ready = false;
                    break;
                }
                if ( level[ neighbour ] > dependencyLevel )
                    dependencyLevel = level[ neighbour ];
            }
            if ( ! ready ) continue;

            level[ cell ] = dependencyLevel + 1;
            ++newlyAssigned;
            if ( level[ cell ] > currentLevel ) currentLevel = level[ cell ];
        }
        if ( newlyAssigned == 0 )
        {
            throw std::runtime_error(
                "SA LU-SGS sweep dependency graph has a cycle" );
        }
        assigned += newlyAssigned;
    }

    const int nLevels = currentLevel + 1;
    levelOffset.assign( static_cast< std::size_t >( nLevels ) + 1, 0 );
    for ( int cell = 0; cell < nCells; ++ cell )
        ++ levelOffset[ level[ cell ] + 1 ];
    for ( int l = 0; l < nLevels; ++ l )
        levelOffset[ l + 1 ] += levelOffset[ l ];

    levelCells.assign( static_cast< std::size_t >( nCells ), 0 );
    std::vector< int > cursor( levelOffset.begin(), levelOffset.end() - 1 );
    for ( int cell = 0; cell < nCells; ++ cell )
        levelCells[ static_cast< std::size_t >( cursor[ level[ cell ] ] ++ ) ] = cell;
}

// Per-direction schedules, cached alongside the CSR.
struct SweepSchedule
{
    const int * cachedLeft = nullptr;
    const int * cachedRight = nullptr;
    int cachedCells = -1;
    int cachedFaces = -1;
    bool valid = false;
    std::vector< int > offsets;
    std::vector< int > ids;
    std::vector< int > levelOffset[ 2 ];
    std::vector< int > levelCells[ 2 ];
};

SweepSchedule & Schedule()
{
    static SweepSchedule schedule;
    return schedule;
}

}

bool UTurbLusgs::SweepHip( bool lowerSweep )
{
#ifdef ONEFLOW_ENABLE_HIP
    if ( turbcom.nEqu != 1 ) return false;
    if ( ! HipProductionFeatureEnabled( "ONEFLOW_ENABLE_UTURB_HIP_LUSGS" ) )
    {
        return false;
    }
    if ( turblu.numberOfSweeps != 1 )
    {
        // The multi-sweep branches carry a drhs accumulation the device kernel
        // does not reproduce; declining keeps those configurations correct.
        return false;
    }
    UnsGrid * grid = Zone::GetUnsGrid();
    if ( grid == nullptr ) return false;
    if ( uturbf.dq == nullptr || uturbf.rhs == nullptr
         || uturbf.impsr == nullptr || uturbf.matrix_l == nullptr
         || uturbf.matrix_r == nullptr || uturbf.meanFlow.q == nullptr )
    {
        return false;
    }

    ug.Init();
    const int nCells = ug.nCells;
    const int nFaces = ug.nFaces;
    if ( nCells <= 0 || nFaces <= 0 ) return false;

    const int * leftCell = ug.lcf->data();
    const int * rightCell = ug.rcf->data();
    const int nBoundaryFaces = ug.nBFaces;

    SweepSchedule & schedule = Schedule();
    if ( ! schedule.valid
         || schedule.cachedLeft != leftCell
         || schedule.cachedRight != rightCell
         || schedule.cachedCells != nCells
         || schedule.cachedFaces != nFaces )
    {
        // Cell-to-face CSR.  A cell's faces come out in ascending face id, the
        // order the host accumulates them in.
        schedule.offsets.assign( static_cast< std::size_t >( nCells ) + 1, 0 );
        for ( int face = 0; face < nFaces; ++ face )
        {
            const int left = leftCell[ face ];
            if ( left >= 0 && left < nCells ) ++ schedule.offsets[ left + 1 ];
            if ( face >= nBoundaryFaces )
            {
                const int right = rightCell[ face ];
                if ( right >= 0 && right < nCells )
                    ++ schedule.offsets[ right + 1 ];
            }
        }
        for ( int cell = 0; cell < nCells; ++ cell )
            schedule.offsets[ cell + 1 ] += schedule.offsets[ cell ];

        schedule.ids.assign(
            static_cast< std::size_t >( schedule.offsets[ nCells ] ), 0 );
        std::vector< int > cursor(
            schedule.offsets.begin(), schedule.offsets.end() - 1 );
        for ( int face = 0; face < nFaces; ++ face )
        {
            const int left = leftCell[ face ];
            if ( left >= 0 && left < nCells )
                schedule.ids[ cursor[ left ] ++ ] = face;
            if ( face >= nBoundaryFaces )
            {
                const int right = rightCell[ face ];
                if ( right >= 0 && right < nCells )
                    schedule.ids[ cursor[ right ] ++ ] = face;
            }
        }

        ComputeSweepLevels(
            schedule.offsets, schedule.ids, leftCell, rightCell,
            nCells, nBoundaryFaces, true,
            schedule.levelOffset[ 0 ], schedule.levelCells[ 0 ] );
        ComputeSweepLevels(
            schedule.offsets, schedule.ids, leftCell, rightCell,
            nCells, nBoundaryFaces, false,
            schedule.levelOffset[ 1 ], schedule.levelCells[ 1 ] );

        schedule.cachedLeft = leftCell;
        schedule.cachedRight = rightCell;
        schedule.cachedCells = nCells;
        schedule.cachedFaces = nFaces;
        schedule.valid = true;
    }

    const int direction = lowerSweep ? 0 : 1;
    const std::vector< int > & useOffset = schedule.levelOffset[ direction ];
    const std::vector< int > & useCells = schedule.levelCells[ direction ];

    const int nLevels = static_cast< int >( useOffset.size() ) - 1;
    if ( nLevels <= 0 ) return false;

    TurbSaSweepView view;
    view.nCells = nCells;
    view.nFaces = nFaces;
    view.nLevels = nLevels;
    // delta spans the ghost range: the upper sweep multiplies a boundary
    // face's constant ghost delta by matrix_r, as the host does.
    view.nGhostCells = ug.nTCell - nCells;
    view.cellFaceOffsets = schedule.offsets.data();
    view.cellFaceIds = schedule.ids.data();
    view.leftCell = leftCell;
    view.rightCell = rightCell;
    view.levelOffset = useOffset.data();
    view.levelCells = useCells.data();
    view.rhs = ( * uturbf.rhs )[ ISA ].data();
    view.implicitSpectrum = ( * uturbf.impsr )[ 0 ].data();
    view.matrixLeft = ( * uturbf.matrix_l )[ 0 ].data();
    view.matrixRight = ( * uturbf.matrix_r )[ 0 ].data();
    view.delta = ( * uturbf.dq )[ ISA ].data();
    view.lowerSweep = lowerSweep;
    view.numberOfSweeps = turblu.numberOfSweeps;
    view.cacheKey = grid;

    HipTurbulenceBackend::Shared().CalcSweep( view );
    return true;
#else
    ( void ) lowerSweep;
    return false;
#endif
}

EndNameSpace
