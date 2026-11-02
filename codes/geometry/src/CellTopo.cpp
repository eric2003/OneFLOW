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

BeginNameSpace( ONEFLOW )


void CellTopo::Alloc( int nCells )
{
    blank.resize( nCells );
    eTypes.resize( nCells );
    elements.resize( nCells );

    blank = 1;
}

void CellTopo::PushElement( int p1, int p2, int p3, int elementType )
{
    IntField element;
    element.push_back( p1 );
    element.push_back( p2 );
    element.push_back( p3 );

    elements.push_back( element );
    eTypes.push_back( elementType );
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

        elements.push_back( element1 );
        elements.push_back( element2 );

        eTypes.push_back( elementType );
        eTypes.push_back( elementType );
    }
    else if ( elementType == ONEFLOW::QUAD_4 )
    {
        IntField element;
        element.push_back( p1 );
        element.push_back( p2 );
        element.push_back( p3 );
        element.push_back( p4 );

        elements.push_back( element );
        eTypes.push_back( elementType );
    }
}

void CellTopo::CalcC2f( FaceTopo & faceTopo )
{
    const HXSize_t nCells = this->GetCellCount();
    const HXSize_t nBFaces = faceTopo.GetBoundaryFaceCount();
    const HXSize_t nFaces = faceTopo.GetFaceCount();
    const IntField & leftCells = faceTopo.GetLeftCells();
    const IntField & rightCells = faceTopo.GetRightCells();

    if ( nBFaces > nFaces ||
         leftCells.size() != nFaces ||
         rightCells.size() != nFaces )
    {
        throw std::logic_error(
            "CellTopo::CalcC2f: inconsistent face-cell topology sizes" );
    }

    // Build derived adjacency separately so invalid input cannot leave a
    // partially populated cache behind.
    LinkField nextC2f;
    nextC2f.resize( nCells + nBFaces );

    for ( HXSize_t iFace = 0; iFace < nFaces; ++ iFace )
    {
        const int leftCell = leftCells[ iFace ];
        if ( leftCell < 0 ||
             static_cast< HXSize_t >( leftCell ) >= nCells )
        {
            throw std::out_of_range(
                "CellTopo::CalcC2f: left cell ID is out of range" );
        }

        nextC2f[ leftCell ].push_back( static_cast< int >( iFace ) );

        if ( iFace >= nBFaces )
        {
            const int rightCell = rightCells[ iFace ];
            if ( rightCell < 0 ||
                 static_cast< HXSize_t >( rightCell ) >= nCells )
            {
                throw std::out_of_range(
                    "CellTopo::CalcC2f: internal-face right cell ID is out of range" );
            }
            if ( rightCell == leftCell )
            {
                throw std::logic_error(
                    "CellTopo::CalcC2f: internal face must connect distinct cells" );
            }

            nextC2f[ rightCell ].push_back( static_cast< int >( iFace ) );
        }
    }

    c2f = nextC2f;
}

void CellTopo::CalcC2C( FaceTopo & faceTopo )
{
    faceTopo.CalcC2C( this->c2c );
}

void CalcC2f( UnsGrid & grid )
{
    FaceTopo & faceTopo = grid.GetFaceTopo();
    CellTopo & cellTopo = grid.GetCellMesh().GetCellTopo();
    cellTopo.CalcC2f( faceTopo );
}

void CalcC2C( UnsGrid & grid )
{
    FaceTopo & faceTopo = grid.GetFaceTopo();
    CellTopo & cellTopo = grid.GetCellMesh().GetCellTopo();
    cellTopo.CalcC2C( faceTopo );
}

EndNameSpace
