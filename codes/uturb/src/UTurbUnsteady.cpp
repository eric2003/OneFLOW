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

#include "UTurbUnsteady.h"
#include "SolverDef.h"
#include "UsdData.h"
#include "Com.h"
#include "UCom.h"
#include "UnsGrid.h"
#include "Zone.h"
#include "DataBase.h"
#include "TurbCom.h"
#include "UTurbCom.h"
#include "NsIdx.h"

BeginNameSpace( ONEFLOW )

UTurbUsdField::UTurbUsdField()
{
}

UTurbUsdField::~UTurbUsdField()
{
}

void UTurbUsdField::Init()
{
    this->InitBasic( TURB_SOLVER );
}

UTurbUnsteady::UTurbUnsteady()
{
    this->sTid = TURB_SOLVER;
    data = new TurbUsdData();
    field = new UTurbUsdField();
    data->Init();
    field->Init();

    this->srcFun = & UTurbUnstPrepareSrcData;
    this->criFun = & UTurbUnstPrepareCriData;

    ug.Init();
    uturbf.Init();
}

UTurbUnsteady::~UTurbUnsteady()
{
    delete data;
    delete field;
}

void UTurbUnstPrepareSrcData( Unsteady * unst )
{
    UsdData * data = unst->data;
    UsdField * field = unst->field;
    for ( int iEqu = 0; iEqu < data->nEqu; ++ iEqu )
    {
        data->prim [ iEqu ] = ( * field->q  )[ iEqu ][ ug.cId ];
        data->prim1[ iEqu ] = ( * field->q1 )[ iEqu ][ ug.cId ];
        data->prim2[ iEqu ] = ( * field->q2 )[ iEqu ][ ug.cId ];
    }

    gcom.cvol  = ( * ug.cvol  )[ ug.cId ];
    gcom.cvol1 = ( * ug.cvol1 )[ ug.cId ];
    gcom.cvol2 = ( * ug.cvol2 )[ ug.cId ];

    Real coef = 1.0;

    if ( data->nEqu >= 2 )
    {
        coef  = ( * uturbf.q_ns )[ IDX::IR ][ ug.cId ];
    }

    for ( int iEqu = 0; iEqu < data->nEqu; ++ iEqu )
    {
        data->q [ iEqu ] = coef * data->prim [ iEqu ];
        data->q1[ iEqu ] = coef * data->prim1[ iEqu ];
        data->q2[ iEqu ] = coef * data->prim2[ iEqu ];
    }
}

void UTurbUnstPrepareCriData( Unsteady * unst )
{
    UsdData * data = unst->data;
    UsdField * field = unst->field;
    for ( int iEqu = 0; iEqu < data->nEqu; ++ iEqu )
    {
        data->prim [ iEqu ] = ( * field->q  )[ iEqu ][ ug.cId ];
        data->prim1[ iEqu ] = ( * field->q1 )[ iEqu ][ ug.cId ];
        data->prim2[ iEqu ] = ( * field->q2 )[ iEqu ][ ug.cId ];
    }

    gcom.cvol  = ( * ug.cvol  )[ ug.cId ];
    gcom.cvol1 = ( * ug.cvol1 )[ ug.cId ];
    gcom.cvol2 = ( * ug.cvol2 )[ ug.cId ];

    Real coef = 1.0;

    if ( data->nEqu >= 2 )
    {
        coef  = ( * uturbf.q_ns )[ IDX::IR ][ ug.cId ];
    }

    for ( int iEqu = 0; iEqu < data->nEqu; ++ iEqu )
    {
        data->q [ iEqu ] = coef * data->prim [ iEqu ];
        data->q1[ iEqu ] = coef * data->prim1[ iEqu ];
        data->q2[ iEqu ] = coef * data->prim2[ iEqu ];
    }

    for ( int iEqu = 0; iEqu < data->nEqu; ++ iEqu )
    {
        data->res [ iEqu ] = ( * field->res  )[ iEqu ][ ug.cId ];
    }

}

EndNameSpace