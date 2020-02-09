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

UNsUsdField::UNsUsdField()
{
}

UNsUsdField::~UNsUsdField()
{
}

void UNsUsdField::Init()
{
    this->InitBasic( NS_SOLVER );
}

UNsUnsteady::UNsUnsteady()
{
    this->sTid = NS_SOLVER;
    data = new NsUsdData();
    field = new UNsUsdField();
    data->Init();
    field->Init();

    this->srcFun = & UNsUnstPrepareSrcData;
    this->criFun = & UNsUnstPrepareCriData;

    ug.Init();
    unsf.Init();
}

UNsUnsteady::~UNsUnsteady()
{
    delete data;
    delete field;
}


void UNsUnstPrepareSrcData( Unsteady * unst )
{
    UsdData * data = unst->data;
    UsdField * field = unst->field;
    for ( int iEqu = 0; iEqu < data->nEqu; ++ iEqu )
    {
        data->prim [ iEqu ] = ( * field->q  )[ iEqu ][ ug.cId ];
        data->prim1[ iEqu ] = ( * field->q1 )[ iEqu ][ ug.cId ];
        data->prim2[ iEqu ] = ( * field->q2 )[ iEqu ][ ug.cId ];
    }
    nscom.gama = ( * unsf.gama  )[ 0 ][ ug.cId ];
    gcom.cvol  = ( * ug.cvol  )[ ug.cId ];
    gcom.cvol1 = ( * ug.cvol1 )[ ug.cId ];
    gcom.cvol2 = ( * ug.cvol2 )[ ug.cId ];

    PrimToQ( data->prim , nscom.gama, data->q  );
    PrimToQ( data->prim1, nscom.gama, data->q1 );
    PrimToQ( data->prim2, nscom.gama, data->q2 );
}

void UNsUnstPrepareCriData( Unsteady * unst )
{
    UsdData * data = unst->data;
    UsdField * field = unst->field;
    for ( int iEqu = 0; iEqu < data->nEqu; ++ iEqu )
    {
        data->prim [ iEqu ] = ( * field->q  )[ iEqu ][ ug.cId ];
        data->prim1[ iEqu ] = ( * field->q1 )[ iEqu ][ ug.cId ];
        data->prim2[ iEqu ] = ( * field->q2 )[ iEqu ][ ug.cId ];
    }

    nscom.gama = ( * unsf.gama  )[ 0 ][ ug.cId ];

    for ( int iEqu = 0; iEqu < data->nEqu; ++ iEqu )
    {
        data->res [ iEqu ] = ( * field->res  )[ iEqu ][ ug.cId ];
    }

    PrimToQ( data->prim , nscom.gama, data->q  );
    PrimToQ( data->prim1, nscom.gama, data->q1 );
    PrimToQ( data->prim2, nscom.gama, data->q2 );
}

EndNameSpace