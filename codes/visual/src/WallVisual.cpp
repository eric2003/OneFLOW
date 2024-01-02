/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2024 He Xin and the OneFLOW contributors.
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

#include "WallVisual.h"
#include "Visualize.h"
#include "UnitElement.h"
#include "ElementHome.h"
#include "HXMath.h"
#include "HXCgns.h"
#include "Dimension.h"
#include <algorithm>
#include <fstream>
#include <iostream>


BeginNameSpace( ONEFLOW )

void PointDebug::AddPoint( Real x, Real y, Real z )
{
    this->xArray.push_back( x );
    this->yArray.push_back( y );
    this->zArray.push_back( z );
}

void PointDebug::GetPoint( int index, Real & x, Real & y, Real & z )
{
    x = this->xArray[ index ];
    y = this->yArray[ index ];
    z = this->zArray[ index ];
}

WallVisual::WallVisual()
{
    faceSet = new std::set< HXSort< IntField > >();
}

WallVisual::~WallVisual()
{
    delete faceSet;
}

void WallVisual::PushElement( int p1, int p2, int elementType )
{
    IntField element;
    element.push_back( p1 );
    element.push_back( p2 );

    this->eLink.push_back( element );
    this->elementType.push_back( elementType );
}

void WallVisual::PushElement( int p1, int p2, int p3, int elementType )
{
    IntField element;
    element.push_back( p1 );
    element.push_back( p2 );
    element.push_back( p3 );

    this->eLink.push_back( element );
    this->elementType.push_back( elementType );
}

void WallVisual::PushElement( int p1, int p2, int p3, int p4, int elementType )
{
    if ( elementType == TRI_3 )
    {
        IntField element1;
        element1.push_back( p1 );
        element1.push_back( p2 );
        element1.push_back( p3 );

        IntField element2;
        element2.push_back( p3 );
        element2.push_back( p4 );
        element2.push_back( p1 );

        this->eLink.push_back( element1 );
        this->eLink.push_back( element2 );

        this->elementType.push_back( elementType );
        this->elementType.push_back( elementType );
    }
    else if ( elementType == QUAD_4 )
    {
        IntField element;
        element.push_back( p1 );
        element.push_back( p2 );
        element.push_back( p3 );
        element.push_back( p4 );

        this->eLink.push_back( element );
        this->elementType.push_back( elementType );
    }
}

void WallVisual::BuildFaceTopo( IntField & faceNodeIndexArray, int loc_Face, int iCell, int face_type )
{
    IntField faceNodeIndexArraySort = faceNodeIndexArray;
    std::sort( faceNodeIndexArraySort.begin(), faceNodeIndexArraySort.end() );

    HXSort< IntField > faceForSorting;
    faceForSorting.value = faceNodeIndexArraySort;

    std::set< HXSort< IntField > >::iterator iter = faceSet->find( faceForSorting );
    if ( iter == faceSet->end() )
    {
        HXSize_t oldFaceNumber = faceSet->size();
        faceForSorting.index = oldFaceNumber;

        faceSet->insert( faceForSorting );

        HXSize_t fId = oldFaceNumber;
        HXSize_t newFaceNumber = fId + 1;
        lCell.resize( newFaceNumber );
        rCell.resize( newFaceNumber );
        lPos.resize( newFaceNumber );
        rPos.resize( newFaceNumber );
        faceType.resize( newFaceNumber );
        faceType[ fId ] = face_type;
        lCell[ fId ] = iCell;
        rCell[ fId ] = INVALID_INDEX;
        lPos[ fId ] = loc_Face;
        rPos[ fId ] = INVALID_INDEX;

        this->fLink.push_back( faceNodeIndexArray );
    }
    else
    {
        int fId = iter->index;
        rCell[ fId ] = iCell;
        rPos[ fId ] = loc_Face;
    }
}

void WallVisual::CalcOrderMap( int & nBFaces, IntField & orderMapping )
{
    int nFaces = rCell.size();
    orderMapping.resize( nFaces );
    int iCount = 0;
    nBFaces = 0;
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        int rc = rCell[ iFace ];
        if ( rc == INVALID_INDEX )
        {
            orderMapping[ iCount ++ ] = iFace;
            ++ nBFaces;
        }
    }

    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        int rc = rCell[ iFace ];
        if ( rc != INVALID_INDEX )
        {
            orderMapping[ iCount ++ ] = iFace;
        }
    }
}

void WallVisual::ConstructTopology2D()
{
}

void WallVisual::ConstructTopology3D()
{
    HXSize_t nCells = this->eLink.size();

    HXSort< IntField > faceForSorting;

    for ( HXSize_t iCell = 0; iCell < nCells; ++ iCell )
    {
        IntField & element = this->eLink[ iCell ];

        int eType = elementType[ iCell ];

        UnitElement * unitElement = ElementHome::GetUnitElement( eType );

        int elem_nFace = unitElement->GetElementFaceNumber();

        for ( int loc_Face = 0; loc_Face < elem_nFace; ++ loc_Face )
        {
            IntField & local_fn = unitElement->GetElementFace( loc_Face );
            int face_type = unitElement->GetFaceType( loc_Face );
            int nNodes = local_fn.size();
            IntField fn;
            for ( int iFacePoint = 0; iFacePoint < nNodes; ++ iFacePoint )
            {
                fn.push_back( element[ local_fn[ iFacePoint ] ] );
            }

            BuildFaceTopo( fn, loc_Face, iCell, face_type );
        }
    }

    IntField orderMapping;
    int nBFaces = 0;
    CalcOrderMap( nBFaces, orderMapping );

    ONEFLOW::Reorder( this->lCell, orderMapping );
    ONEFLOW::Reorder( this->rCell, orderMapping );
    ONEFLOW::Reorder( this->lPos, orderMapping );
    ONEFLOW::Reorder( this->rPos, orderMapping );
    ONEFLOW::Reorder( this->fLink, orderMapping );
    ONEFLOW::Reorder( this->faceType, orderMapping );

    for ( int iFace = 0; iFace < nBFaces; ++ iFace )
    {
        rCell[ iFace ] = iFace + nCells;
    }
}

void WallVisual::ConstructTopology()
{
    if ( Dim::dimension == THREE_D )
    {
        this->ConstructTopology3D();
    }
}

void WallVisual::Visual( std::fstream & file, StringField & titleOfTecplot, RealField2D & qNodeField )
{
    if ( Dim::dimension == TWO_D )
    {
        this->VisualLine( file, titleOfTecplot, qNodeField );
    }
    else
    {
        this->Visual3D( file, titleOfTecplot, qNodeField );
    }
}

void WallVisual::Visual3D( std::fstream & file, StringField & titleOfTecplot, RealField2D & qNodeField )
{
    std::ostringstream oss;
    for ( HXSize_t i = 0; i < titleOfTecplot.size(); ++ i )
    {
        oss << titleOfTecplot[ i ] << "\n";
    }

    int totalNumFaceNodes = GetTotalNumFaceNodes( this->fLink );
    int numberOfFaces = this->fLink.size();

    int numberOfNodes = xN.size();
    int numberOfCells = this->eLink.size();

    oss << "ZONE\n";
    oss << "ZoneType = FEPolygon\n";

    oss << "Nodes    = " << numberOfNodes << "\n";
    oss << "Faces    = " << numberOfFaces << "\n";
    oss << "Elements = " << numberOfCells << "\n";
    oss << "TotalNumFaceNodes = " << totalNumFaceNodes << "\n";
    oss << "NumConnectedBoundaryFaces = 0\n";
    oss << "TotalNumBoundaryConnections = 0\n";

    Plot::oss = & oss;
    Plot::DumpField( xN );
    Plot::DumpField( yN );
    Plot::DumpField( zN );

    HXSize_t nEqu = qNodeField.size();

    for ( HXSize_t iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        Plot::DumpField( qNodeField[ iEqu ] );
    }

    Plot::DumpFaceNodeLink( this->fLink );
    Plot::DumpFaceElementLink( lCell, numberOfCells );
    Plot::DumpFaceElementLink( rCell, numberOfCells );

    file << oss.str();
}

void WallVisual::VisualLine( std::fstream & file, StringField & titleOfTecplot, RealField2D & qNodeField )
{
    std::ostringstream oss;
    for ( HXSize_t i = 0; i < titleOfTecplot.size(); ++ i )
    {
        oss << titleOfTecplot[ i ] << "\n";
    }

    int numberOfFaces = this->fLink.size();
    int numberOfNodes = xN.size();
    int numberOfCells = this->eLink.size();

    oss << "ZONE I = " << numberOfNodes << " F = BLOCK " << "\n";

    Plot::oss = & oss;
    Plot::DumpField( xN );
    Plot::DumpField( yN );
    Plot::DumpField( zN );

    int nEqu = qNodeField.size();

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        Plot::DumpField( qNodeField[ iEqu ] );
    }
    file << oss.str();
}

EndNameSpace
