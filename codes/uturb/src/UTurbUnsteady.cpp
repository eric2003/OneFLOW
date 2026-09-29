/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\---------------------------------------------------------------------------*/

#include "UTurbUnsteady.h"
#include "SolverDef.h"
#include "Com.h"
#include "UCom.h"
#include "UnsGrid.h"
#include "Zone.h"
#include "DataBase.h"
#include "TurbCom.h"
#include "UTurbCom.h"
#include "NsIdx.h"

BeginNameSpace( ONEFLOW )

UTurbUnsteady::UTurbUnsteady()
{
    this->SetEquationCount( turbcom.nEqu );

    this->srcFun = & UTurbUnstPrepareSrcData;
    this->criFun = & UTurbUnstPrepareCriData;

    ug.Init();
    uturbf.Init();
}
void UTurbUnstPrepareSrcData( UUnsteady * unsteady )
{
    UnsteadyFieldView * field = &unsteady->field;

    MRField * q =
        field->GetFlow( UnsteadyFieldView::HistoryLevel::Current );

    MRField * q1 =
        field->GetFlow( UnsteadyFieldView::HistoryLevel::Previous );

    MRField * q2 =
        field->GetFlow( UnsteadyFieldView::HistoryLevel::Old );

    for ( int iEqu = 0; iEqu < unsteady->GetEquationCount(); ++ iEqu )
    {
        unsteady->prim[ iEqu ] =
            ( * q )[ iEqu ][ ug.cId ];

        unsteady->prim1[ iEqu ] =
            ( * q1 )[ iEqu ][ ug.cId ];

        unsteady->prim2[ iEqu ] =
            ( * q2 )[ iEqu ][ ug.cId ];
    }

    gcom.cvol  = ( * ug.cvol  )[ ug.cId ];
    gcom.cvol1 = ( * ug.cvol1 )[ ug.cId ];
    gcom.cvol2 = ( * ug.cvol2 )[ ug.cId ];

    Real coef = 1.0;

    if ( unsteady->GetEquationCount() >= 2 )
    {
        coef  = ( * uturbf.q_ns )[ IDX::IR ][ ug.cId ];
    }

    for ( int iEqu = 0; iEqu < unsteady->GetEquationCount(); ++ iEqu )
    {
        unsteady->q [ iEqu ] = coef * unsteady->prim [ iEqu ];
        unsteady->q1[ iEqu ] = coef * unsteady->prim1[ iEqu ];
        unsteady->q2[ iEqu ] = coef * unsteady->prim2[ iEqu ];
    }
}

void UTurbUnstPrepareCriData( UUnsteady * unsteady )
{
    UnsteadyFieldView * field = &unsteady->field;

    MRField * q =
        field->GetFlow( UnsteadyFieldView::HistoryLevel::Current );

    MRField * q1 =
        field->GetFlow( UnsteadyFieldView::HistoryLevel::Previous );

    MRField * q2 =
        field->GetFlow( UnsteadyFieldView::HistoryLevel::Old );

    for ( int iEqu = 0; iEqu < unsteady->GetEquationCount(); ++ iEqu )
    {
        unsteady->prim [ iEqu ] =
            ( * q )[ iEqu ][ ug.cId ];

        unsteady->prim1[ iEqu ] =
            ( * q1 )[ iEqu ][ ug.cId ];

        unsteady->prim2[ iEqu ] =
            ( * q2 )[ iEqu ][ ug.cId ];
    }

    gcom.cvol  = ( * ug.cvol  )[ ug.cId ];
    gcom.cvol1 = ( * ug.cvol1 )[ ug.cId ];
    gcom.cvol2 = ( * ug.cvol2 )[ ug.cId ];

    Real coef = 1.0;

    if ( unsteady->GetEquationCount() >= 2 )
    {
        coef  = ( * uturbf.q_ns )[ IDX::IR ][ ug.cId ];
    }

    for ( int iEqu = 0; iEqu < unsteady->GetEquationCount(); ++ iEqu )
    {
        unsteady->q [ iEqu ] = coef * unsteady->prim [ iEqu ];
        unsteady->q1[ iEqu ] = coef * unsteady->prim1[ iEqu ];
        unsteady->q2[ iEqu ] = coef * unsteady->prim2[ iEqu ];
    }

}

EndNameSpace
