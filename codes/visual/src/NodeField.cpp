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

#include "NodeField.h"
#include "UnsGrid.h"
#include "Zone.h"
#include "DataBase.h"
#include "FaceTopo.h"
#include "BcRecord.h"
#include "Boundary.h"

BeginNameSpace( ONEFLOW )

MRField * AllocNodeVar( int nEqu )
{
    UnsGrid * grid = Zone::GetUnsGrid();
    int nNode = grid->nNode;
    MRField * nf = new MRField( nEqu, nNode );
    return nf;
}

MRField * CreateNodeVar( const string & name )
{
    UnsGrid * grid = Zone::GetUnsGrid();
    MRField * cf = GetFieldPointer< MRField > ( grid, name );
    int nNode = grid->nNode;
    int nEqu = cf->GetNEqu();

    MRField * nf = AllocNodeVar( nEqu );
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        CalcNodeVar( ( * nf )[ iEqu ], ( * cf )[ iEqu ] );
    }
    return nf;
}

MRField * CreateNodeVar( RealField & qc )
{
    UnsGrid * grid = Zone::GetUnsGrid();
    MRField * fn = AllocNodeVar( 1 );
    CalcNodeVar( ( * fn )[ 0 ], qc );
    return fn;
}

void CalcNodeVar( RealField & qNodeField, RealField & qField )
{
    UnsGrid * grid = Zone::GetUnsGrid();
    FaceTopo * faceTopo = grid->faceTopo;
    LinkField & f2c = faceTopo->f2n;

    int nNode = grid->nNode;
    int nFace = grid->nFace;
    RealField nCount( nNode );
    nCount = 0.0;
    qNodeField = 0.0;

    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        int lc = faceTopo->lCell[ iFace ];
        int rc = faceTopo->rCell[ iFace ];

        int fnNode = f2c[ iFace ].size();
        for ( int iNode = 0; iNode < fnNode; ++ iNode )
        {
            int nodeId = f2c[ iFace ][ iNode ];

            qNodeField[ nodeId ] += qField[ lc ];
            nCount[ nodeId ] += 1;

            qNodeField[ nodeId ] += qField[ rc ];
            nCount[ nodeId ] += 1;
        }
    }  

    FixBcNodeVar( qNodeField, qField, nCount, BC::SYMMETRY     , true );
    FixBcNodeVar( qNodeField, qField, nCount, BC::SOLID_SURFACE, true );
    FixBcNodeVar( qNodeField, qField, nCount, BC::INTERFACE    , true );
    FixBcNodeVar( qNodeField, qField, nCount, BC::FARFIELD     , true );

    for ( int iNode = 0; iNode < nNode; ++ iNode )
    {
        qNodeField[ iNode ] /= ( nCount[ iNode ] + SMALL );
    }
}

void FixBcNodeVar( RealField & qNodeField, RealField & qField, RealField & nCount, int bcType, bool twoSide )
{
    UnsGrid * grid = Zone::GetUnsGrid();
    FaceTopo * faceTopo = grid->faceTopo;
    BcRecord * bcRecord = faceTopo->bcManager->bcRecord;
    LinkField & f2c = faceTopo->f2n;

    int nNode = grid->nNode;
    int nFace = grid->nFace;
    int nBFace = grid->nBFace;

    for ( int iFace = 0; iFace < nBFace; ++ iFace )
    {
        if ( bcRecord->bcType[ iFace ] != bcType ) continue;
        int lc = faceTopo->lCell[ iFace ];
        int rc = faceTopo->rCell[ iFace ];

        int fnNode = f2c[ iFace ].size();
        for ( int iNode = 0; iNode < fnNode; ++ iNode )
        {
            int nodeIndex = f2c[ iFace ][ iNode ];
            qNodeField [ nodeIndex ] = 0.0;
            nCount[ nodeIndex ] = 0;
        }
    }

    for ( int iFace = 0; iFace < nBFace; ++ iFace )
    {
        if ( bcRecord->bcType[ iFace ] != bcType ) continue;

        int lc = faceTopo->lCell[ iFace ];
        int rc = faceTopo->rCell[ iFace ];

        int fnNode = f2c[ iFace ].size();
        for ( int iNode = 0; iNode < fnNode; ++ iNode )
        {
            int nodeIndex = f2c[ iFace ][ iNode ];
            qNodeField [ nodeIndex ] += qField[ lc ];
            nCount[ nodeIndex ] += 1;
            if ( twoSide )
            {
                qNodeField[ nodeIndex ] += qField[ rc ];
                nCount[ nodeIndex ] += 1;
            }
        }
    }
}

EndNameSpace