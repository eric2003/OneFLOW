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

#include "CellMesh.h"
#include "NodeMesh.h"
#include "FaceMesh.h"
#include "FaceTopo.h"
#include "BcRecord.h"
#include "CellTopo.h"
#include "UnsGrid.h"
#include "ActionState.h"
#include "DataBook.h"
#include "DataBaseIO.h"
#include "HXMath.h"


using namespace std;

BeginNameSpace( ONEFLOW )

CellMesh::CellMesh()
{
    cellTopo = new CellTopo();
}

CellMesh::~CellMesh()
{
    delete cellTopo;
}

UInt CellMesh::GetNumberOfCells()
{ 
    return cellTopo->GetNumberOfCells(); 
}

void CellMesh::AllocateMetrics( FaceMesh * faceMesh )
{
    UInt numberOfCells = this->GetNumberOfCells();
    UInt nBFace = faceMesh->GetNBFace();
    UInt nTCell = numberOfCells + nBFace;
    this->xcc.resize( nTCell );
    this->ycc.resize( nTCell );
    this->zcc.resize( nTCell );
    this->vol.resize( nTCell );

}

void CellMesh::AllocDist()
{
    UInt numberOfCells = this->GetNumberOfCells();
    dist.resize( numberOfCells );
}

void CellMesh::ReadDist()
{
    ActionState::dataBook->MoveToBegin();
    HXRead( ActionState::dataBook, dist );
}

void CellMesh::DumpDist()
{
    ActionState::dataBook->MoveToBegin();
    HXWrite( ActionState::dataBook, dist );
}

void CellMesh::CalcCellSpan( UnsGrid * grid )
{
    if ( this->span.size() ) return;
    int nCell = this->GetNumberOfCells();
    this->span.resize( nCell );
    CalcC2f( grid );

    FaceTopo * faceTopo = grid->faceTopo;
    LinkField & c2f = this->cellTopo->c2f;
    IntField & lcf = faceTopo->lCell;
    IntField & rcf = faceTopo->rCell;

    for ( int ic = 0; ic < nCell; ++ ic )
    {
        Real delta  = - LARGE;
        Real deltaj = 0.0;

        int fn = c2f[ ic ].size();

        for ( int iFace = 0; iFace < fn; ++ iFace )
        {
            int fId = c2f[ ic ][ iFace ];

            int lc = lcf[ fId ];
            int rc = rcf[ fId ];

            Real dx = xcc[ rc ] - xcc[ lc ];
            Real dy = ycc[ rc ] - ycc[ lc ];
            Real dz = zcc[ rc ] - zcc[ lc ];
            deltaj = DIST( dx, dy, dz );
            delta  = MAX( deltaj , delta );
        }

        span[ ic ] = delta;
    }
}

void CalcCellSpan( UnsGrid * grid )
{
    grid->cellMesh->CalcCellSpan( grid );
}

EndNameSpace