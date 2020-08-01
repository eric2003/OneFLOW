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

#include "CellTopo.h"
#include "CellMesh.h"
#include "FaceTopo.h"
#include "UnsGrid.h"
#include "HXCgns.h"
#include "DataBase.h"

BeginNameSpace( ONEFLOW )


CellTopo::CellTopo()
{
    ;
}

CellTopo::~CellTopo()
{
    ;
}

void CellTopo::Alloc( int nCell )
{
    blank.resize( nCell );
    cellType.resize( nCell );
    cellToNode.resize( nCell );

    blank = 1;
}

void CellTopo::PushElement( int p1, int p2, int p3, int elementType )
{
    IntField element;
    element.push_back( p1 );
    element.push_back( p2 );
    element.push_back( p3 );

    cellToNode.push_back( element );
    cellType.push_back( elementType );
}

void CellTopo::PushElement( int p1, int p2, int p3, int p4, int elementType )
{
    if ( elementType == ONEFLOW::TRI_3 )
    {
        IntField element1;
        element1.push_back( p1 );
        element1.push_back( p2 );
        element1.push_back( p3 );

        IntField element2;
        element2.push_back( p3 );
        element2.push_back( p4 );
        element2.push_back( p1 );

        cellToNode.push_back( element1 );
        cellToNode.push_back( element2 );

        cellType.push_back( elementType );
        cellType.push_back( elementType );
    }
    else if ( elementType == ONEFLOW::QUAD_4 )
    {
        IntField element;
        element.push_back( p1 );
        element.push_back( p2 );
        element.push_back( p3 );
        element.push_back( p4 );

        cellToNode.push_back( element );
        cellType.push_back( elementType );
    }
}

void CellTopo::CalcC2f( FaceTopo * faceTopo )
{
    if ( c2f.size() != 0 ) return;

    int nCell = this->GetNumberOfCells();
    int nBFace = faceTopo->GetNBFace();
    int nFace = faceTopo->GetNFace();
	int startStrategy = ONEFLOW::GetDataValue< int >("startStrategy");

	if (startStrategy == 2|| startStrategy == 3)
	{
		c2f.resize(nCell + nBFace);  //ins时做了修改
	}
	else
	{
		c2f.resize(nCell);
	}

    for ( int iFace = 0; iFace < nBFace; ++ iFace )
    {
        int lc  = faceTopo->lCell[ iFace ];
        c2f[ lc  ].push_back( iFace );
    }

    for ( int iFace = nBFace; iFace < nFace; ++ iFace )
    {
        int lc  = faceTopo->lCell[ iFace ];
        int rc  = faceTopo->rCell[ iFace ];
        c2f[ lc ].push_back( iFace );
        c2f[ rc ].push_back( iFace );
    }
}

void CellTopo::CalcC2C( FaceTopo * faceTopo )
{
    faceTopo->CalcC2C( this->c2c );
}

void CalcC2f( UnsGrid * grid )
{
    FaceTopo * faceTopo = grid->faceTopo;
    CellTopo * cellTopo = grid->cellMesh->cellTopo;
    cellTopo->CalcC2f( faceTopo );
}

void CalcC2C( UnsGrid * grid )
{
    FaceTopo * faceTopo = grid->faceTopo;
    CellTopo * cellTopo = grid->cellMesh->cellTopo;
    cellTopo->CalcC2C( faceTopo );
}

EndNameSpace