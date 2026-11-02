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

#include "Grad.h"
#include "StageProfiler.h"

#include <cstdlib>
#include "UGrad.h"
#include "UCom.h"
#include "FieldManager.h"
#include "InterField.h"
#include "FaceTopo.h"
#include "BcRecord.h"
#include "UnsGrid.h"
#include "FaceMesh.h"
#include "CellMesh.h"
#include "CellTopo.h"
#include "Zone.h"
#include "Iteration.h"
#include "Ctrl.h"
#include "DataBase.h"
#include "StringUtils.h"
#include <iostream>
#include <iomanip>


BeginNameSpace( ONEFLOW )

void Grad::CalcGrad()
{
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        ONEFLOW::CalcGradGGCellWeight( ( * q )[ iEqu ], ( * dqdx )[ iEqu ], ( * dqdy )[ iEqu ], ( * dqdz )[ iEqu ] );
    }

    if ( ctrl.debugFieldDump != 0 )
    {
        Real mindiff = 1.0e-10;
        int idumpface = 1;
        int idumpcell = 0;
        const MRField * qField = q;
        const MRField * dqdxField = dqdx;
        const MRField * dqdyField = dqdy;
        const MRField * dqdzField = dqdz;

        std::string fname = AddString( "grad.", name, ".debug" );
        std::cout << "varname = " << name << "\n";
        HXDebug::DumpField( fname, qField );
        HXDebug::CompareFile( 1.0e-12, idumpcell );

        // Literal "q": debug always dumps NS mean-flow reference, even when
        // Grad::name is turbq/tempr. See FieldWrap.h (cross-solver path).
        const UnsGrid * grid = Zone::GetConstUnsGrid();
        const MRField * qq = GetFieldPointer< MRField >( grid, "q" );
        HXDebug::DumpField( "flowq.debug", qq );
        HXDebug::CompareFile( 1.0e-12, idumpcell );

        std::string fnamedqdx = AddString( "grad.", name, ".dqdx.debug" );
        std::string fnamedqdy = AddString( "grad.", name, ".dqdy.debug" );
        std::string fnamedqdz = AddString( "grad.", name, ".dqdz.debug" );

        HXDebug::DumpField( fnamedqdx, dqdxField );
        HXDebug::CompareFile( mindiff, idumpcell );
        HXDebug::DumpField( fnamedqdy, dqdyField );
        HXDebug::CompareFile( mindiff, idumpcell );
        HXDebug::DumpField( fnamedqdz, dqdzField );
        HXDebug::CompareFile( mindiff, idumpcell );
    }

    this->SwapBcGrad();
}

void Grad::CalcGradDebug()
{
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        ONEFLOW::CalcGradDebug( ( * q )[ iEqu ], ( * dqdx )[ iEqu ], ( * dqdy )[ iEqu ], ( * dqdz )[ iEqu ] );
    }

    this->SwapBcGrad();
}

namespace
{

// StoreBcGrad fills the face-indexed boundary gradient (bdqdx/bdqdy/bdqdz).
// No consumer reads those fields in this tree, and the store costs a random
// gather over 3*neq cell-sized arrays per RK stage (about 4 ms/step on the
// 8-zone DCU case). Keep it behind an explicit switch so a future consumer can
// re-enable it instead of silently reading zeros.
bool StoreBoundaryGradientRequested()
{
    static const bool enabled = []()
    {
        const char * value = std::getenv( "ONEFLOW_ENABLE_BC_GRADIENT_STORE" );
        return value != nullptr && value[ 0 ] == '1';
    }();
    return enabled;
}

}

void Grad::SwapBcGrad()
{
    UnsGrid * grid = Zone::GetUnsGrid();

    {
        ScopedStageTimer timer( "swap_upload" );
    UploadInterfaceValue( grid, dqdx, namex, nEqu );
    UploadInterfaceValue( grid, dqdy, namey, nEqu );
    UploadInterfaceValue( grid, dqdz, namez, nEqu );
    }

    {
        ScopedStageTimer timer( "swap_download" );
    DownloadInterfaceValue( grid, dqdx, namex, nEqu );
    DownloadInterfaceValue( grid, dqdy, namey, nEqu );
    DownloadInterfaceValue( grid, dqdz, namez, nEqu );
    }

    if ( ctrl.debugFieldDump != 0 )
    {
        Real mindiff = 1.0e-10;
        int idumpface = 1;
        int idumpcell = 0;
        const MRField * dqdxField = dqdx;
        const MRField * dqdyField = dqdy;
        const MRField * dqdzField = dqdz;

        HXDebug::DumpField( "grad.dqdxDownload.debug", dqdxField );
        HXDebug::CompareFile( mindiff, idumpcell );
        HXDebug::DumpField( "grad.dqdyDownload.debug", dqdyField );
        HXDebug::CompareFile( mindiff, idumpcell );
        HXDebug::DumpField( "grad.dqdzDownload.debug", dqdzField );
        HXDebug::CompareFile( mindiff, idumpcell );
    }
    if ( istore == 1 && StoreBoundaryGradientRequested() )
    {
        ScopedStageTimer timer( "swap_store" );
        this->StoreBcGrad();
    }
}

void Grad::StoreBcGrad()
{
    UnsGrid * grid = Zone::GetUnsGrid();

    FaceTopo * faceTopo = &grid->GetFaceTopo();

    IntField & bcType = faceTopo->GetBcRecord().bcType;

    int nBFaces = bcType.size();

    for ( int iFace = 0; iFace < nBFaces; ++ iFace )
    {
        int bc_type = bcType[ iFace ];

        int lc = faceTopo->GetLeftCells()[ iFace ];
        int rc = faceTopo->GetRightCells()[ iFace ];

        for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
        {
            ( * bdqdx )[ iEqu ][ iFace ] = ( * dqdx )[ iEqu ][ lc ];
            ( * bdqdy )[ iEqu ][ iFace ] = ( * dqdy )[ iEqu ][ lc ];
            ( * bdqdz )[ iEqu ][ iFace ] = ( * dqdz )[ iEqu ][ lc ];
        }
    }

}

EndNameSpace
