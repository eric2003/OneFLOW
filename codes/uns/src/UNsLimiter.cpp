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

#include "UNsLimiter.h"
#include "UnsGrid.h"
#include "Ctrl.h"
#include "Zone.h"
#include "DataBase.h"
#include "UCom.h"
#include "NsCom.h"
#include "UNsCom.h"
#include "Boundary.h"
#include "BcRecord.h"

BeginNameSpace( ONEFLOW )

NsLimField::NsLimField()
{
    this->nEqu = nscom.nEqu;
}

void NsLimField::Init()
{
    // Literal NS names only (also used as UTurbInvFlux::nslimiter under TURB).
    UnsGrid * grid = Zone::GetUnsGrid();
    BindFieldsAndAllocFaces( grid, "q", "dqdx", "dqdy", "dqdz", "limiter" );
    this->ckfun = & NsCheckFunction;
}

void NsLimField::BcQlQrFix()
{
    for ( int fId = 0; fId < ug.nBFaces; ++ fId )
    {
        int bcType = ug.bcRecord->bcType[ fId ];
        if ( bcType == BC::INTERFACE ) continue;
        if ( bcType == BC::PERIODIC ) continue;

        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        for ( int iEqu = 0; iEqu < this->nEqu; ++ iEqu )
        {
            Real tmp = half * ( ( * this->q )[ iEqu ][ ug.lc ] + ( * this->q )[ iEqu ][ ug.rc ] );

            ( * this->qf1 )[ iEqu ][ ug.fId ] = tmp;
            ( * this->qf2 )[ iEqu ][ ug.fId ] = tmp;
        }

        // Literal bc_q - do not use unsf.bc_q (may be stale under non-NS SolverState).
        if ( bcType == BC::SOLID_SURFACE )
        {
            MRField * bcQ = GetFieldPointer< MRField >( Zone::GetUnsGrid(), "bc_q" );
            if ( bcQ == nullptr ) continue;
            for ( int iEqu = 0; iEqu < this->nEqu; ++ iEqu )
            {
                ( * this->qf1 )[ iEqu ][ ug.fId ] = ( * bcQ )[ iEqu ][ ug.fId ];
                ( * this->qf2 )[ iEqu ][ ug.fId ] = ( * bcQ )[ iEqu ][ ug.fId ];
            }
        }
    }
}

NsLimiter::NsLimiter()
{
    InjectLimField( std::make_unique<NsLimField>() );
    limflag = ctrl.ilim;
}

EndNameSpace
