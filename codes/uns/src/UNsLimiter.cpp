/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2020 He Xin and the OneFLOW contributors.
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
    qf1 = 0;
    qf2 = 0;
    this->nEqu = nscom.nEqu;
}

NsLimField::~NsLimField()
{
    delete qf1;
    delete qf2;
}

void NsLimField::Init()
{
    UnsGrid * grid = Zone::GetUnsGrid();

    q       = GetFieldPointer< MRField > ( grid, "q" );
    dqdx    = GetFieldPointer< MRField > ( grid, "dqdx" );
    dqdy    = GetFieldPointer< MRField > ( grid, "dqdy" );
    dqdz    = GetFieldPointer< MRField > ( grid, "dqdz" );
    limiter = GetFieldPointer< MRField > ( grid, "limiter" );

    this->nEqu = q->GetNEqu();

    qf1 = new MRField( this->nEqu, grid->nFace );
    qf2 = new MRField( this->nEqu, grid->nFace );

    this->ckfun = & NsCheckFunction;
}

void NsLimField::BcQlQrFix()
{
    for ( int fId = 0; fId < ug.nBFace; ++ fId )
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

        if ( bcType == BC::SOLID_SURFACE )
        {
            for ( int iEqu = 0; iEqu < this->nEqu; ++ iEqu )
            {
                ( * this->qf1 )[ iEqu ][ ug.fId ] = ( * unsf.bc_q )[ iEqu ][ ug.fId ];
                ( * this->qf2 )[ iEqu ][ ug.fId ] = ( * unsf.bc_q )[ iEqu ][ ug.fId ];
            }
        }
    }
}

NsLimiter::NsLimiter()
{
    limf = new NsLimField();
    limflag = ctrl.ilim;
}

NsLimiter::~NsLimiter()
{
    delete limf;
}

EndNameSpace