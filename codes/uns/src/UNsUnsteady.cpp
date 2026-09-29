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
    this->SetEquationCount( nscom.nTEqu );

    this->SetSourceFunction( & UNsUnstPrepareSrcData );
    this->SetCriterionFunction( & UNsUnstPrepareCriData );

    ug.Init();
    unsf.Init();
}

void UNsUnstPrepareSrcData( UUnsteady * unsteady )
{
    MRField * q =
        unsteady->GetFlow( Unsteady::HistoryLevel::Current );

    MRField * q1 =
        unsteady->GetFlow( Unsteady::HistoryLevel::Previous );

    MRField * q2 =
        unsteady->GetFlow( Unsteady::HistoryLevel::Old );

    RealField & primitive =
        unsteady->GetPrimitive(
            Unsteady::HistoryLevel::Current );

    RealField & primitive1 =
        unsteady->GetPrimitive(
            Unsteady::HistoryLevel::Previous );

    RealField & primitive2 =
        unsteady->GetPrimitive(
            Unsteady::HistoryLevel::Old );

    RealField & conservative =
        unsteady->GetConservative(
            Unsteady::HistoryLevel::Current );

    RealField & conservative1 =
        unsteady->GetConservative(
            Unsteady::HistoryLevel::Previous );

    RealField & conservative2 =
        unsteady->GetConservative(
            Unsteady::HistoryLevel::Old );

    for ( int iEqu = 0; iEqu < unsteady->GetEquationCount(); ++ iEqu )
    {
        primitive[ iEqu ] =
            ( * q )[ iEqu ][ ug.cId ];

        primitive1[ iEqu ] =
            ( * q1 )[ iEqu ][ ug.cId ];

        primitive2[ iEqu ] =
            ( * q2 )[ iEqu ][ ug.cId ];
    }
    nscom.gama = ( * unsf.gama  )[ 0 ][ ug.cId ];
    gcom.cvol  = ( * ug.cvol  )[ ug.cId ];
    gcom.cvol1 = ( * ug.cvol1 )[ ug.cId ];
    gcom.cvol2 = ( * ug.cvol2 )[ ug.cId ];

    PrimToQ( primitive , nscom.gama, conservative  );
    PrimToQ( primitive1, nscom.gama, conservative1 );
    PrimToQ( primitive2, nscom.gama, conservative2 );
}

void UNsUnstPrepareCriData( UUnsteady * unsteady )
{
    MRField * q =
        unsteady->GetFlow( Unsteady::HistoryLevel::Current );

    MRField * q1 =
        unsteady->GetFlow( Unsteady::HistoryLevel::Previous );

    MRField * q2 =
        unsteady->GetFlow( Unsteady::HistoryLevel::Old );

    RealField & primitive =
        unsteady->GetPrimitive(
            Unsteady::HistoryLevel::Current );

    RealField & primitive1 =
        unsteady->GetPrimitive(
            Unsteady::HistoryLevel::Previous );

    RealField & primitive2 =
        unsteady->GetPrimitive(
            Unsteady::HistoryLevel::Old );

    RealField & conservative =
        unsteady->GetConservative(
            Unsteady::HistoryLevel::Current );

    RealField & conservative1 =
        unsteady->GetConservative(
            Unsteady::HistoryLevel::Previous );

    RealField & conservative2 =
        unsteady->GetConservative(
            Unsteady::HistoryLevel::Old );

    for ( int iEqu = 0; iEqu < unsteady->GetEquationCount(); ++ iEqu )
    {
        primitive [ iEqu ] =
            ( * q )[ iEqu ][ ug.cId ];

        primitive1[ iEqu ] =
            ( * q1 )[ iEqu ][ ug.cId ];

        primitive2[ iEqu ] =
            ( * q2 )[ iEqu ][ ug.cId ];
    }

    nscom.gama = ( * unsf.gama  )[ 0 ][ ug.cId ];


    PrimToQ( primitive , nscom.gama, conservative  );
    PrimToQ( primitive1, nscom.gama, conservative1 );
    PrimToQ( primitive2, nscom.gama, conservative2 );
}

EndNameSpace
