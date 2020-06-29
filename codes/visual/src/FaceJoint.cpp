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

#include "FaceJoint.h"
#include "PointSearch.h"
#include "WallVisual.h"
#include "HXCgns.h"
#include "HXMath.h"
#include "Tolerence.h"

BeginNameSpace( ONEFLOW )

FaceJointManager::FaceJointManager()
{
    global = new FaceJoint();
}

FaceJointManager::~FaceJointManager()
{
    delete global;
    int nZone = this->patch.size();
    for ( int zId = 0; zId < nZone; ++ zId )
    {
        delete this->patch[ zId ];
    }
}

void FaceJointManager::ConstructPointIndex()
{
    this->global->ConstructPointIndex();

    int nLocal = this->patch.size();
    for ( int iLocal = 0; iLocal < nLocal; ++ iLocal )
    {
        FaceJoint * local = this->patch[ iLocal ];
        local->ConstructPointIndex();
        local->ConstructPointIndexMap( this->global );
    }
}

void FaceJointManager::CalcNodeValue()
{
    this->global->CalcNodeValue();

    int nLocal = this->patch.size();
    for ( int iLocal = 0; iLocal < nLocal; ++ iLocal )
    {
        FaceJoint * local = this->patch[ iLocal ];
        local->RemapNodeValue( this->global );
    }
}

FaceJoint::FaceJoint()
{
    pmin.resize( 3 );
    pmax.resize( 3 );
    ps = new PointSearch();
    wallVisual = new WallVisual();
    isValid = false;
}

FaceJoint::~FaceJoint()
{
    delete ps;
    delete wallVisual;
}

void FaceJoint::CalcBoundBox()
{
    pmin[ 0 ] = LARGE;
    pmin[ 1 ] = LARGE;
    pmin[ 2 ] = LARGE;

    pmax[ 0 ] = - LARGE;
    pmax[ 1 ] = - LARGE;
    pmax[ 2 ] = - LARGE;

    FaceJoint::PointLink & fvp = this->fvp;
    int numberOfWallFaces = this->GetSize();

    for ( int iFace = 0; iFace < numberOfWallFaces; ++ iFace )
    {
        FaceJoint::PointField & ptList = fvp[ iFace ];
        int numberOfNodes = ptList.size();
        for ( int iNode = 0; iNode < numberOfNodes; ++ iNode )
        {
            FaceJoint::PointType & point = ptList[ iNode ];
            pmin[ 0 ] = MIN( pmin[ 0 ], point.x );
            pmin[ 1 ] = MIN( pmin[ 1 ], point.y );
            pmin[ 2 ] = MIN( pmin[ 2 ], point.z );

            pmax[ 0 ] = MAX( pmax[ 0 ], point.x );
            pmax[ 1 ] = MAX( pmax[ 1 ], point.y );
            pmax[ 2 ] = MAX( pmax[ 2 ], point.z );
        }
    }

    dismin = LARGE;
    dismax = - LARGE;

    Real ptTol = Tolerence::GetTol();

    for ( int iFace = 0; iFace < numberOfWallFaces; ++ iFace )
    {
        FaceJoint::PointField & pointArray = fvp[ iFace ];
        int numberOfNodes = pointArray.size();
        for ( int iNode = 0; iNode < numberOfNodes - 1; ++ iNode )
        {
            int p1 = iNode;
            FaceJoint::PointType & point1 = pointArray[ p1 ];
            for ( int jNode = iNode + 1; jNode < numberOfNodes; ++ jNode )
            {
                int p2 = jNode;
                FaceJoint::PointType & point2 = pointArray[ p2 ];
                Real dx = point2.x - point1.x;
                Real dy = point2.y - point1.y;
                Real dz = point2.z - point1.z;
                Real ds = DIST( dx, dy, dz );

                if ( ds <= ptTol ) continue;
                dismin = MIN( dismin, ds );
                dismax = MAX( dismax, ds );
            }
        }
    }
}

void FaceJoint::ConstructPointIndex()
{
    if ( ! this->isValid ) return;
    this->CalcBoundBox();

    FaceJoint::PointLink & fvp = this->fvp;

    Real tolerance = this->dismin / 4;

    this->ps->Initialize( this->pmin, this->pmax, tolerance );

    RealField coor( 3 );
    LinkField & fLink = this->fLink;

    int nWallFace = this->GetSize();

    for ( int iFace = 0; iFace < nWallFace; ++ iFace )
    {
        FaceJoint::PointField & ptList = fvp[ iFace ];
        int numberOfNodes = ptList.size();
        IntField face;
        for ( int iNode = 0; iNode < numberOfNodes; ++ iNode )
        {
            FaceJoint::PointType & point = ptList[ iNode ];
            int index = this->ps->AddPoint( point.x, point.y, point.z );
            face.push_back( index );
        }
        fLink.push_back( face );
    }

    int numberOfPoints = this->ps->GetNPoint();

    IntField & weightId = this->weightId;
    weightId.resize( numberOfPoints, 0 );

    RealField & fnv = this->fnv;
    fnv.resize( numberOfPoints, 0 );
}

void FaceJoint::ConstructPointIndexMap( FaceJoint * globalBasicWall )
{
    if ( ! this->isValid ) return;

    int numberOfPoints = this->ps->GetNPoint();
    this->l2g.resize( numberOfPoints );
    for ( int iPoint = 0; iPoint < numberOfPoints; ++ iPoint )
    {
        Real xm, ym, zm;
        this->ps->GetPoint( iPoint, xm, ym, zm );
        int index = globalBasicWall->ps->AddPoint( xm, ym, zm );
        this->l2g[ iPoint ] = index;
    }
}

void FaceJoint::CalcNodeValue()
{
    int numberOfWallFaces = this->GetSize();

    LinkField & fLink = this->fLink;
    FaceJoint::PointLink & faceVertexArray = this->fvp;

    IntField & weightId = this->weightId;
    RealField & fnv = this->fnv;
    RealField & fcv = this->fcv;

    RealField coor( 3 );

    for ( int iFace = 0; iFace < numberOfWallFaces; ++ iFace )
    {
        FaceJoint::PointField & pointArray = faceVertexArray[ iFace ];
        int numberOfNodes = pointArray.size();
        for ( int iNode = 0; iNode < numberOfNodes; ++ iNode )
        {
            FaceJoint::PointType & point = pointArray[ iNode ];
            int index = this->ps->AddPoint( point.x, point.y, point.z );
            fnv[ index ] += fcv[ iFace ];
            weightId[ index ] ++;
        }
    }

    int nNode = fnv.size();
    for ( int iNode = 0; iNode < nNode; ++ iNode )
    {
        fnv[ iNode ] /= weightId[ iNode ];
    }
}

void FaceJoint::RemapNodeValue( FaceJoint * globalBasicWall )
{
    if ( ! this->isValid ) return;

    int numberOfPoints = this->ps->GetNPoint();
    for ( int iPoint = 0; iPoint < numberOfPoints; ++ iPoint )
    {
        int index = this->l2g[ iPoint ];
        this->fnv[ iPoint ] = globalBasicWall->fnv[ index ];
    }
}

void FaceJoint::AddFacePoint( int nSolidCell, FaceJoint::PointLink & ptLink )
{
    FaceJoint::PointLink & fvp = this->fvp;

    for ( int iElement = 0; iElement < nSolidCell; ++ iElement )
    {
        fvp.push_back( ptLink[ iElement ] );
    }
}

void FaceJoint::AddFaceCenterValue( int nSolidCell, RealField & fcvIn )
{
    RealField & fcv = this->fcv;

    for ( int iElement = 0; iElement < nSolidCell; ++ iElement )
    {
        fcv.push_back( fcvIn[ iElement ] );
    }
}

void FaceJoint::Visual( fstream & file )
{
    if ( ! this->isValid ) return;
    LinkField & fLink = this->fLink;
    RealField & fnv = this->fnv;

    int count2 = 0;
    int count3 = 0;
    int count4 = 0;

    int numberOfWallFaces = this->GetSize();
    for ( int iFace = 0; iFace < numberOfWallFaces; ++ iFace )
    {
        IntField & facePointIndex = fLink[ iFace ];
        int numberOfNodes = facePointIndex.size();
        IntField diffId;
        int p0 = facePointIndex[ 0 ];
        diffId.push_back( facePointIndex[ 0 ] );
        for ( int iNode = 1; iNode < numberOfNodes; ++ iNode )
        {
            int p1 = facePointIndex[ iNode - 1 ];
            int p2 = facePointIndex[ iNode ];
            if ( p2 == p1 ) continue;
            if ( iNode == numberOfNodes - 1 )
            {
                if ( p2 == p0 )
                {
                    continue;
                }
            }
            diffId.push_back( p2 );
        }
        if ( diffId.size() == 2 )
        {
            this->wallVisual->PushElement( diffId[ 0 ], diffId[ 1 ],  BAR_2 );
            ++ count2;
        }
        else if ( diffId.size() == 3 )
        {
            this->wallVisual->PushElement( diffId[ 0 ], diffId[ 1 ], diffId[ 2 ], TRI_3 );
            ++ count3;
        }
        else if ( diffId.size() == 4 )
        {
            this->wallVisual->PushElement( diffId[ 0 ], diffId[ 1 ], diffId[ 2 ], diffId[ 3 ], QUAD_4 );
            ++ count4;
        }
    }
    //cout << " number of TRI_3 is " << count3 << "\n";
    //cout << " number of QUAD_4 is " << count4 << "\n";
    //cout << " number of total faces are " << numberOfWallFaces << "\n";

    int numberOfPoints = this->ps->GetNPoint();

    for ( int iNode = 0; iNode < numberOfPoints; ++ iNode )
    {
        Real xm, ym, zm;

        this->ps->GetPoint( iNode, xm, ym, zm );

        this->wallVisual->xN.push_back( xm );
        this->wallVisual->yN.push_back( ym );
        this->wallVisual->zN.push_back( zm );
    }

    this->wallVisual->ConstructTopology();

    StringField titleOfTecplot;
    titleOfTecplot.push_back( "title=\"THE FLOW FIELD OF ONEFLOW\"" );
    titleOfTecplot.push_back( "variables=" );
    titleOfTecplot.push_back( "\"x\"" );
    titleOfTecplot.push_back( "\"y\"" );
    titleOfTecplot.push_back( "\"z\"" );
    titleOfTecplot.push_back( "\"heatflux\"" );

    int numberOfEquations = 1;
    RealField2D qNodeField;
    qNodeField.push_back( fnv );
    this->wallVisual->Visual( file, titleOfTecplot, qNodeField );
}


EndNameSpace