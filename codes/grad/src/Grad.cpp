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

#include "Grad.h"
#include "UGrad.h"
#include "UCom.h"
#include "FieldImp.h"
#include "FaceTopo.h"
#include "BcRecord.h"
#include "UnsGrid.h"
#include "FaceMesh.h"
#include "CellMesh.h"
#include "CellTopo.h"
#include "Zone.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

Grad::Grad()
{
    istore = 0;
}

Grad::~Grad()
{
    ;
}

void Grad::CalcGrad()
{
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        //ONEFLOW::CalcGrad( ( * q )[ iEqu ], ( * dqdx )[ iEqu ], ( * dqdy )[ iEqu ], ( * dqdz )[ iEqu ] );
        ONEFLOW::CalcGradGGCellWeight( ( * q )[ iEqu ], ( * dqdx )[ iEqu ], ( * dqdy )[ iEqu ], ( * dqdz )[ iEqu ] );
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

void Grad::SwapBcGrad()
{
    UnsGrid * grid = Zone::GetUnsGrid();

    UploadInterfaceValue( grid, dqdx, namex, nEqu );
    UploadInterfaceValue( grid, dqdy, namey, nEqu );
    UploadInterfaceValue( grid, dqdz, namez, nEqu );

    DownloadInterfaceValue( grid, dqdx, namex, nEqu );
    DownloadInterfaceValue( grid, dqdy, namey, nEqu );
    DownloadInterfaceValue( grid, dqdz, namez, nEqu );
    if ( istore == 1 )
    {
        this->StoreBcGrad();
    }
}

void Grad::StoreBcGrad()
{
    UnsGrid * grid = Zone::GetUnsGrid();

    FaceTopo * faceTopo = grid->faceTopo;

    IntField & bcType = faceTopo->bcManager->bcRecord->bcType;

    int nBFace = bcType.size();

    for ( int iFace = 0; iFace < nBFace; ++ iFace )
    {
        int bc_type = bcType[ iFace ];

        int lc = faceTopo->lCell[ iFace ];
        int rc = faceTopo->rCell[ iFace ];

        for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
        {
            ( * bdqdx )[ iEqu ][ iFace ] = ( * dqdx )[ iEqu ][ lc ];
            ( * bdqdy )[ iEqu ][ iFace ] = ( * dqdy )[ iEqu ][ lc ];
            ( * bdqdz )[ iEqu ][ iFace ] = ( * dqdz )[ iEqu ][ lc ];
        }
    }

}

EndNameSpace