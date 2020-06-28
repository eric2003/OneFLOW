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

#include "UCom.h"
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

UGeom ug;

UGeom::UGeom()
{
    ;
}

UGeom::~UGeom()
{
    ;
}

void UGeom::Init()
{
    //UnsGrid * grid = Zone::GetUnsGrid();
    grid = Zone::GetUnsGrid();

    ug.nBFace = grid->nBFace;
    ug.nCell = grid->nCell;
    ug.nTCell = grid->nCell + grid->nBFace;
    ug.nFace = grid->nFace;

    this->SetStEd( F_TOTAL );
    this->CreateBcRegion();

    FaceTopo * faceTopo = grid->faceTopo;
    ug.lcf = & faceTopo->lCell;
    ug.rcf = & faceTopo->rCell;

    FaceMesh * faceMesh = grid->faceMesh;
    CellMesh * cellMesh = grid->cellMesh;
    CellTopo * cellTopo = grid->cellMesh->cellTopo;

    ug.xfn = & faceMesh->xfn;
    ug.yfn = & faceMesh->yfn;
    ug.zfn = & faceMesh->zfn;
    ug.vfn = & faceMesh->vfn;
    ug.farea = & faceMesh->area;

    ug.vfx = & faceMesh->vfx;
    ug.vfy = & faceMesh->vfy;
    ug.vfz = & faceMesh->vfz;

    ug.xfc = & faceMesh->xfc;
    ug.yfc = & faceMesh->yfc;
    ug.zfc = & faceMesh->zfc;

    ug.xcc = & cellMesh->xcc;
    ug.ycc = & cellMesh->ycc;
    ug.zcc = & cellMesh->zcc;

    ug.cvol  = & cellMesh->vol;
    ug.cvol1 = & cellMesh->vol;
    ug.cvol2 = & cellMesh->vol;

    ug.blankf = & cellTopo->blank;

    cellTopo->CalcC2f( faceTopo );

    ug.c2f = & cellTopo->c2f;

    //ug.ireconface = 0;
    ug.ireconface = 1;
}

void UGeom::CreateBcRegion()
{
    UnsGrid * grid = Zone::GetUnsGrid();
    BcRecord * bcRecord = grid->faceTopo->bcManager->bcRecord;
    bcRecord->CreateBcRegion();

    ug.bcRecord = bcRecord;
}

void UGeom::SetStEd( int flag )
{
    if ( flag == F_INNER )
    {
        this->ist = 0;
        this->ied = ug.nCell;
    }
    else if ( flag == F_GHOST )
    {
        this->ist = ug.nCell;
        this->ied = ug.nTCell;
    }
    else if ( flag == F_TOTAL )
    {
        this->ist = 0;
        this->ied = ug.nTCell;
    }
}

void UGeom::DumpCellFace( int cId )
{
    int nFace = ( * this->c2f )[ cId ].size();
    for ( int fId = 0; fId < nFace; ++ fId )
    {
        cout << ( * this->c2f )[ cId ][ fId ] << " ";
    }
    cout << endl;
}

void AddF2CField( MRField * cellField, MRField * faceField )
{
    int nEqu = cellField->GetNEqu();
    for ( int fId = 0; fId < ug.nBFace; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
        {
            ( * cellField )[ iEqu ][ ug.lc ] -= ( * faceField )[ iEqu ][ ug.fId ];
        }
    }

    for ( int fId = ug.nBFace; fId < ug.nFace; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
        {
            ( * cellField )[ iEqu ][ ug.lc ] -= ( * faceField )[ iEqu ][ ug.fId ];
            ( * cellField )[ iEqu ][ ug.rc ] += ( * faceField )[ iEqu ][ ug.fId ];
        }
    }
}

void AddF2CFieldDebug( MRField * cellField, MRField * faceField )
{
    int nEqu = cellField->GetNEqu();
    for ( int fId = 0; fId < ug.nBFace; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        if ( ug.lc == 9 || ug.rc == 9 )
        {
            cout << " fId = " << fId << "\n";
            int iEqu = 0;
            cout << ( * faceField )[ iEqu ][ ug.fId ] << "\n";
        }

        for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
        {
            ( * cellField )[ iEqu ][ ug.lc ] -= ( * faceField )[ iEqu ][ ug.fId ];
        }
        if ( ug.lc == 9 || ug.rc == 9 )
        {
            for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
            {
                int cc = 9;
                cout << ( * cellField )[ iEqu ][ cc ] << "\n";
            }
        }
    }

    for ( int fId = ug.nBFace; fId < ug.nFace; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];
        if ( ug.lc == 9 || ug.rc == 9 )
        {
            cout << " fId = " << fId << "\n";
            int iEqu = 0;
            cout << ( * faceField )[ iEqu ][ ug.fId ] << "\n";
        }

        for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
        {
            ( * cellField )[ iEqu ][ ug.lc ] -= ( * faceField )[ iEqu ][ ug.fId ];
            ( * cellField )[ iEqu ][ ug.rc ] += ( * faceField )[ iEqu ][ ug.fId ];
        }
        if ( ug.lc == 9 || ug.rc == 9 )
        {
            for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
            {
                int cc = 9;
                cout << ( * cellField )[ iEqu ][ cc ] << "\n";
            }
        }
    }

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        int cc = 9;
        cout << ( * cellField )[ iEqu ][ cc ] << "\n";
    }
    int kkk = 1;
}

EndNameSpace