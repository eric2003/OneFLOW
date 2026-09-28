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

    OneFLOW is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\---------------------------------------------------------------------------*/

#include "UNsUnsteady.h"
#include "NsUnsteady.h"
#include "SolverDef.h"
#include "Iteration.h"
#include "UNsCom.h"
#include "UCom.h"
#include "NsCom.h"
#include "NsInvFlux.h"
#include "UnsGrid.h"
#include "Zone.h"
#include "DataBase.h"

BeginNameSpace( ONEFLOW )

UNsUnsteady::UNsUnsteady()
{
    this->solverType = NS_SOLVER;
    field = new UsdField();
    this->SetEquationCount( nscom.nTEqu );
    field->InitBasic( NS_SOLVER );

    this->srcFun = & UNsUnstPrepareSrcData;
    this->criFun = & UNsUnstPrepareCriData;

    ug.Init();
    unsf.Init();
}

UNsUnsteady::~UNsUnsteady()
{
    delete field;
}


void UNsUnstPrepareSrcData( Unsteady * unst )
{
    UUnsteady * unsteady = static_cast< UUnsteady * >( unst );
    UsdField * field = unst->field;

    MRField * q =
        field->GetFlow( UsdField::HistoryLevel::Current );

    MRField * q1 =
        field->GetFlow( UsdField::HistoryLevel::Previous );

    MRField * q2 =
        field->GetFlow( UsdField::HistoryLevel::Old );

    for ( int iEqu = 0; iEqu < unsteady->nEqu; ++ iEqu )
    {
        unsteady->prim[ iEqu ] =
            ( * q )[ iEqu ][ ug.cId ];

        unsteady->prim1[ iEqu ] =
            ( * q1 )[ iEqu ][ ug.cId ];

        unsteady->prim2[ iEqu ] =
            ( * q2 )[ iEqu ][ ug.cId ];
    }
    nscom.gama = ( * unsf.gama  )[ 0 ][ ug.cId ];
    gcom.cvol  = ( * ug.cvol  )[ ug.cId ];
    gcom.cvol1 = ( * ug.cvol1 )[ ug.cId ];
    gcom.cvol2 = ( * ug.cvol2 )[ ug.cId ];

    PrimToQ( unsteady->prim , nscom.gama, unsteady->q  );
    PrimToQ( unsteady->prim1, nscom.gama, unsteady->q1 );
    PrimToQ( unsteady->prim2, nscom.gama, unsteady->q2 );
}

void UNsUnstPrepareCriData( Unsteady * unst )
{
    UUnsteady * unsteady = static_cast< UUnsteady * >( unst );
    UsdField * field = unst->field;

    MRField * q =
        field->GetFlow( UsdField::HistoryLevel::Current );

    MRField * q1 =
        field->GetFlow( UsdField::HistoryLevel::Previous );

    MRField * q2 =
        field->GetFlow( UsdField::HistoryLevel::Old );

    for ( int iEqu = 0; iEqu < unsteady->nEqu; ++ iEqu )
    {
        unsteady->prim [ iEqu ] =
            ( * q )[ iEqu ][ ug.cId ];

        unsteady->prim1[ iEqu ] =
            ( * q1 )[ iEqu ][ ug.cId ];

        unsteady->prim2[ iEqu ] =
            ( * q2 )[ iEqu ][ ug.cId ];
    }

    nscom.gama = ( * unsf.gama  )[ 0 ][ ug.cId ];


    PrimToQ( unsteady->prim , nscom.gama, unsteady->q  );
    PrimToQ( unsteady->prim1, nscom.gama, unsteady->q1 );
    PrimToQ( unsteady->prim2, nscom.gama, unsteady->q2 );
}

EndNameSpace
