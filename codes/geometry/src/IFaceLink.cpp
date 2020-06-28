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

#include "IFaceLink.h"
#include "InterFace.h"
#include "Grid.h"
#include "PointSearch.h"
#include "FaceSearch.h"
#include "CgnsPeriod.h"
#include "NodeMesh.h"
#include <algorithm>
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )
IFaceLink::IFaceLink( Grids & grids )
{
    this->grids = grids;

    int nZone = grids.size();
    this->l2g.resize( nZone );

    this->face_search = new FaceSearch();
    this->point_search = new PointSearch();
    this->point_search->Initialize( grids );
}

IFaceLink::~IFaceLink()
{
    ;
}

void IFaceLink::Init( Grid * grid )
{
    int zid = grid->id;
    int nIFace = grid->interFace->nIFace;

    this->l2g[ zid ].resize( nIFace );
}

void IFaceLink::AddFace( const IntField & facePointIndexes )
{
    this->face_search->AddFace( facePointIndexes );
}

void IFaceLink::CreateLink( IntField & faceNode, int zid, int lCount )
{
    int nTIFace = this->gI2Zid.size();

    this->AddFace( faceNode );

    std::sort( faceNode.begin(), faceNode.end() );
    HXSort< IntField > face( faceNode, nTIFace );
    set < HXSort< IntField > >::iterator iter = this->inFaceList.find( face );

    if ( iter == this->inFaceList.end() )
    {
        this->inFaceList.insert( face );
        this->l2g[ zid ][ lCount ] = nTIFace;
        IntField zids;
        IntField lIid;
        zids.push_back( zid );
        lIid.push_back( lCount );

        this->gI2Zid.push_back( zids );
        this->g2l.push_back( lIid );
    }
    else
    {
        int gIid = iter->index;
        this->l2g[ zid ][ lCount ] = gIid;
        this->gI2Zid[ gIid ].push_back( zid );
        this->g2l[ gIid ].push_back( lCount );
        //this->inFaceList.erase( iter ); //??
    }
}

void IFaceLink::ReconstructInterFace()
{
    this->face_search->CalcNewFaceId( this );
}

void IFaceLink::InitNewLgMapping()
{
    this->gI2ZidNew = this->gI2Zid;
    this->g2lNew = this->g2l;
}

void IFaceLink::UpdateLgMapping()
{
    this->gI2Zid = this->gI2ZidNew;
    this->g2l = this->g2lNew;
}

void IFaceLink::MatchInterfaceTopology( Grid * grid )
{
    InterFace * interFace = grid->interFace;
    if ( ! interFace ) return;

    int nPeoridic = 0;

    int nIFace = this->l2g[ grid->id ].size();

    for ( int iIFace = 0; iIFace < nIFace; ++ iIFace )
    {
        int gIFace = this->l2g[ grid->id ][ iIFace ];
        bool flag = false;
        int nIZone = this->gI2Zid[ gIFace ].size();

        if ( nIZone != 2 )
        {
            if ( nIZone > 2 )
            {
                //cout << " More than two faces coincide\n";
            }
            else
            {
                ++nPeoridic;
                //cout << " Less than two faces coincide\n";
            }
            //cout << " Current ZoneIndex  = " << grid->id << endl;
            //cout << " nIZone = " << nIZone << endl;
            //cout << " LocalInterface Index = " << iIFace << " nIFace = " << nIFace << endl;
        }

        for ( int iIZone = 0; iIZone < nIZone; ++ iIZone )
        {
            int nZid = this->gI2Zid [ gIFace ][ iIZone ];
            int lId  = this->g2l[ gIFace ][ iIZone ];
            if ( ( nZid != grid->id ) ||
                 ( lId  != iIFace   ) )
            {
                interFace->zoneId[ iIFace ] = nZid;
                interFace->localInterfaceId[ iIFace ] = lId;
                flag = true;
                break;
            }
        }

        //if ( ! flag )
        //{
        //    cout << "LocalInterface Index = " << iIFace << " There is a problem in the input grid. Please check it carefully!\n";
        //}
    }
    cout << " Total peoridic boundary faces = " << nPeoridic << "\n";
    if ( nPeoridic != 0 )
    {
        //this->MatchPeoridicInterface( grid );
    }
    int kkk = 1;
}

void IFaceLink::MatchPeoridicInterface( Grid * grid )
{
    InterFace * interFace = grid->interFace;
    if ( ! interFace ) return;

    int nPeoridic = 0;

    int nIFace = this->l2g[ grid->id ].size();

    for ( int iIFace = 0; iIFace < nIFace; ++ iIFace )
    {
        int gIFace = this->l2g[ grid->id ][ iIFace ];
        bool flag = false;
        int nIZone = this->gI2Zid[ gIFace ].size();

        if (nIZone == 2) continue;

        int iIZone = 0;

        int nZid = this->gI2Zid [ gIFace ][ iIZone ];
        int lId  = this->g2l[ gIFace ][ iIZone ];

        FaceSort * faceSort = this->face_search->faceArray[ gIFace ];
        IntField & nodeId = faceSort->nodeId;

        RealField xList, yList, zList;
        this->point_search->GetFaceCoorList( nodeId, xList, yList, zList );

        RealField xxList, yyList, zzList;
        f2fmap.FindFace( xList, yList, zList, xxList, yyList, zzList );

        int nNode = xxList.size( );
        IntField faceNode_period;

        for ( int i = 0; i < nNode; ++ i )
        {
            Real xm = xxList[ i ];
            Real ym = yyList[ i ];
            Real zm = zzList[ i ];
            int id = this->point_search->FindPoint( xm, ym, zm );
            faceNode_period.push_back( id );
        }

        int faceId_period = this->face_search->FindFace( faceNode_period );

        int nZid_period = this->gI2Zid [ faceId_period ][ iIZone ];
        int lId_period  = this->g2l[ faceId_period ][ iIZone ];

        interFace->zoneId[ iIFace ] = nZid_period;
        interFace->localInterfaceId[ iIFace ] = lId_period;
        int kkk = 1;
    }
    int kkk = 1;
}

void GetFaceCoorList( IntField & faceNode, RealField & xList, RealField & yList, RealField & zList, NodeMesh * nodeMesh )
{
    int nPoint = faceNode.size();
    for ( int iNode = 0; iNode < nPoint; ++ iNode )
    {
        int gN = faceNode[ iNode ];
        xList[ iNode ] = nodeMesh->xN[ gN ];
        yList[ iNode ] = nodeMesh->yN[ gN ];
        zList[ iNode ] = nodeMesh->zN[ gN ];
    }
}

void GetCoorIdList( IFaceLink * iFaceLink, RealField & xList, RealField & yList, RealField & zList, int nPoint, IntField & pointId )
{
    for ( int iNode = 0; iNode < nPoint; ++ iNode )
    {
        Real xm = xList[ iNode ];
        Real ym = yList[ iNode ];
        Real zm = zList[ iNode ];

        pointId[ iNode ] = iFaceLink->point_search->AddPoint( xm, ym, zm );
    }
}

EndNameSpace