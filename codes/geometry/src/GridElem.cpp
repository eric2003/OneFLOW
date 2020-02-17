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

#include "GridElem.h"
#include "CgnsZone.h"
#include "HXCgns.h"
#include "UnsGrid.h"
#include "HXMath.h"
#include "ElemFeature.h"
#include "FaceTopo.h"
#include "FaceSolver.h"
#include "BcRecord.h"
#include "Boundary.h"
#include "PointFactory.h"
#include "NodeMesh.h"
#include <iostream>
#include <iomanip>
using namespace std;

BeginNameSpace( ONEFLOW )

GridElem::GridElem( HXVector< CgnsZone * > & cgnsZones )
{
    this->minLen = LARGE;
    this->maxLen = -LARGE;
    this->cgnsZones = cgnsZones;
    this->point_factory = new PointFactory();
    this->elem_feature = new ElemFeature();
    this->face_solver = new FaceSolver();
    this->elem_feature->face_solver = face_solver;
}

GridElem::~GridElem()
{
    cout << "delete this->point_factory;\n";
    delete this->point_factory;
    cout << "delete this->elem_feature;\n";
    delete this->elem_feature;
    cout << "delete this->face_solver;\n";
    delete this->face_solver;
    cout << "GridElem::~GridElem()\n";
}

void GridElem::PrepareUnsCompGrid()
{
    cout << " InitCgnsElements()\n";
    this->InitCgnsElements();
    cout << " ScanElements()\n";
    this->elem_feature->ScanElements();
    cout << " ScanBcFace()\n";
    this->ScanBcFace();

    //Continue to parse
    cout << " ScanElements()\n";
    this->elem_feature->ScanElements();
    this->GenerateCompElement();
}

void GridElem::ComputeMinMaxInfo()
{
    this->minLen =   LARGE;
    this->maxLen = - LARGE;
    for ( int iZone = 0; iZone < cgnsZones.size(); ++ iZone )
    {
        CgnsZone * cgnsZone = cgnsZones[ iZone ];

        this->minLen = MIN( this->minLen, cgnsZone->minLen );
        this->maxLen = MAX( this->maxLen, cgnsZone->maxLen );
    }
    cout << "   minLen = " << this->minLen << " maxLen = " << this->maxLen << endl;
}

void GridElem::InitCgnsElements()
{
    for ( int iZone = 0; iZone < cgnsZones.size(); ++ iZone )
    {
        CgnsZone * cgnsZone = cgnsZones[ iZone ];
        cgnsZone->InitElement( this );
    }
}

void GridElem::ScanBcFace()
{
    int nZone = this->cgnsZones.size();
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        this->cgnsZones[ iZone ]->ScanBcFace( this->elem_feature->face_solver );
    }
}

void GridElem::GenerateCompElement()
{
    int nElement =  this->elem_feature->eType->size();

    FaceTopo * faceTopo = this->face_solver->faceTopo;

    int nFace = this->face_solver->faceTopo->f2n.size();

    int nBFace = 0;

    //cout << " nFace = " << nFace << "\n";

    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        if ( iFace % 200000 == 0 ) 
        {
            cout << " iFace = " << iFace << " numberOfTotalFaces = " << nFace << endl;
        }

        int rc = ( faceTopo->rCell )[ iFace ];

        if ( rc == INVALID_INDEX )
        {
            faceTopo->bcManager->bcRecord->bcType.push_back( ( * this->face_solver->faceBcType )[ iFace ] );
            faceTopo->bcManager->bcRecord->bcRegion.push_back( ( * this->face_solver->faceBcKey )[ iFace ] );
            ++ nBFace;
        }
    }

    this->point_factory->InitC2g();

}

void GridElem::GenerateCompGrid( Grid * gridIn )
{
    UnsGrid * grid = UnsGridCast ( gridIn );

    grid->nCell = this->elem_feature->eType->size();
    cout << "   nCell = " << grid->nCell << endl;

    int nNode = this->point_factory->c2g.size();
    grid->nodeMesh->CreateNodes( nNode );
    grid->nNode = nNode;

    for ( int iNode = 0; iNode < nNode; ++ iNode )
    {
        int nodeIndex = this->point_factory->c2g[ iNode ];

        grid->nodeMesh->xN[ iNode ] = this->point_factory->pointList[ nodeIndex ].x;
        grid->nodeMesh->yN[ iNode ] = this->point_factory->pointList[ nodeIndex ].y;
        grid->nodeMesh->zN[ iNode ] = this->point_factory->pointList[ nodeIndex ].z;
    }

    this->ComputeBoundaryType( grid );
    this->ReorderLink( grid );

    cout << "\n-->All the computing information is ready\n";
}

void GridElem::ComputeBoundaryType( UnsGrid * grid )
{
    cout << "\n-->Set boundary condition......\n";
    delete grid->faceTopo;
    grid->faceTopo = this->face_solver->faceTopo;
    grid->faceTopo->grid = grid;
    this->face_solver->faceTopo = 0;
    int nFace = grid->faceTopo->f2n.size();
    cout << " nFace = " << nFace << "\n";
     
    BcRecord * bcRecord = grid->faceTopo->bcManager->bcRecord;
    int nBFace = bcRecord->bcType.size();

    grid->nBFace = nBFace;

    cout << " nBFace = " << nBFace << "\n";

    BcTypeMap * bcTypeMap = new BcTypeMap();
    bcTypeMap->Init();

    IntField cgnsBcArray = bcRecord->bcType;

    IntSet originalBcSet, finalBcSet;
    int iCount = 0;
    for ( int iFace = 0; iFace < nBFace; ++ iFace )
    {
        int cgnsBcType = bcRecord->bcType[ iFace ];
        int bcNameId = bcRecord->bcRegion[ iFace ];
        int bcType = bcTypeMap->Cgns2OneFlow( cgnsBcType );

        bcRecord->bcType[ iCount ] = bcType;

        originalBcSet.insert( cgnsBcType );
        finalBcSet.insert( bcType );
        ++ iCount;
    }

    delete bcTypeMap;

    IntField nBFaceSub;

    for ( IntSet::iterator iter = originalBcSet.begin(); iter != originalBcSet.end(); ++ iter )
    {
        int iCount = 0;
        for ( int iFace = 0; iFace < nBFace; ++ iFace )
        {
            int cgnsBcType = cgnsBcArray[ iFace ];
            if ( cgnsBcType == * iter )
            {
                ++ iCount;
            }
        }
        nBFaceSub.push_back( iCount );
    }

    cout << " Original Boundary Condition Number is " << originalBcSet.size() << endl;
    iCount = 0;
    for ( IntSet::iterator iter = originalBcSet.begin(); iter != originalBcSet.end(); ++ iter )
    {
        int oriBcType = * iter;
        cout << " Boundary Type = " << setw( 3 ) << oriBcType << "  Name = " << setiosflags(ios::left) << setw( 23 ) << GetCgnsBcName( oriBcType );
        cout << " Face = " << nBFaceSub[ iCount ] << endl;
        ++ iCount;
    }
    cout << endl;
    cout << " Final Boundary Condition Number is " << finalBcSet.size() << endl;
    cout << " Boundary Type : ";
    for ( IntSet::iterator iter = finalBcSet.begin(); iter != finalBcSet.end(); ++ iter )
    {
        cout << * iter << " ";
    }
    cout << endl;
}

void GridElem::ReorderLink( UnsGrid * grid )
{
    FaceTopo * faceTopo = grid->faceTopo;
    int nFace = faceTopo->faceType.size();
    grid->nFace = nFace;

    IntField f1map( nFace ), f2map( nFace );
    int iCount = 0;
    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        int rc = faceTopo->rCell[ iFace ];
        if ( rc == INVALID_INDEX )
        {
            f1map[ iFace ] = iCount;
            f2map[ iCount ] = iFace;
            ++ iCount;
        }
    }

    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        int rc = faceTopo->rCell[ iFace ];
        if ( rc != INVALID_INDEX )
        {
            f1map[ iFace ] = iCount;
            f2map[ iCount ] = iFace;
            ++ iCount;
        }
    }
    faceTopo->faceToNodeNew.resize( nFace );
    faceTopo->lCellNew.resize( nFace );
    faceTopo->rCellNew.resize( nFace );
    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        int jFace = f2map[ iFace ];
        faceTopo->faceToNodeNew[ iFace ] = faceTopo->f2n[ jFace ];
        faceTopo->lCellNew[ iFace ] = faceTopo->lCell[ jFace ];
        faceTopo->rCellNew[ iFace ] = faceTopo->rCell[ jFace ];
    }
    faceTopo->f2n = faceTopo->faceToNodeNew;
    faceTopo->lCell = faceTopo->lCellNew;
    faceTopo->rCell = faceTopo->rCellNew;
}

GridElemS::GridElemS()
{
    ;
}

GridElemS::~GridElemS()
{
    for ( int i = 0; i < this->data.size(); ++ i )
    {
        delete this->data[ i ];
    }
}

void GridElemS::AddGridElem( GridElem * gridElem )
{
    this->data.push_back( gridElem );
}

void GridElemS::AddGridElem( HXVector< CgnsZone * > cgnsZones, int iZone )
{
    GridElem * gridElem = new GridElem( cgnsZones );
    this->AddGridElem( gridElem );
}

GridElem * GridElemS::GetGridElem( int iGridElem )
{
    return this->data[ iGridElem ];
}

EndNameSpace