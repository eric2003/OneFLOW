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

#include "IFaceLink.h"
#include "Constant.h"
#include "InterFace.h"
#include "Grid.h"
#include "PointLocator.h"
#include "FaceSearch.h"
#include "CgnsPeriod.h"
#include "NodeMesh.h"
#include <algorithm>
#include <iostream>

BeginNameSpace( ONEFLOW )
IFaceLink::IFaceLink( Grids & grids )
{
    this->grids = grids;

    int nZone = grids.size();
    this->l2g.resize( nZone );

    this->face_search = std::make_unique< FaceSearch >();
    this->point_search = std::make_unique< PointLocator >();
    this->point_search->Initialize( grids );
}

IFaceLink::~IFaceLink() = default;

void IFaceLink::Init( Grid * grid )
{
    int zid = grid->id;
    int nIFaces = grid->interFace->nIFaces;

    this->l2g[ zid ].resize( nIFaces );
}

void IFaceLink::AddFace( const IntField & facePointIndexes )
{
    this->face_search->AddFace( facePointIndexes );
}

void IFaceLink::CreateLink( IntField & faceNode, int zid, int lCount )
{
    // Add face to the face list
    this->AddFace(faceNode);

    // Find or add the face (HXLookup automatically sorts the nodes)
    auto [gIid, isNew] = this->faceLookup.FindOrAdd(faceNode);

    if ( isNew )
    {
        // New face: update local-to-global mapping
        this->l2g[zid][lCount] = gIid;

        // Initialize face connectivity data
        IntField zids;
        IntField lIid;
        zids.push_back(zid);
        lIid.push_back(lCount);
        this->gI2Zid.push_back(std::move(zids));
        this->g2l.push_back(std::move(lIid));
    }
    else
    {
        // Existing face: use the existing global face ID
        this->l2g[zid][lCount] = gIid;
        this->gI2Zid[gIid].push_back(zid);
        this->g2l[gIid].push_back(lCount);
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
    InterFace * interFace = grid->interFace.get();
    if ( ! interFace ) return;

    int missingPeriodicPartnerCount = 0;

    int nIFaces = this->l2g[ grid->id ].size();

    for ( int iIFace = 0; iIFace < nIFaces; ++ iIFace )
    {
        int gIFace = this->l2g[ grid->id ][ iIFace ];
        int nIZone = this->gI2Zid[ gIFace ].size();

        if ( nIZone != 2 )
        {
            if ( nIZone > 2 )
            {
                //std::cout << " More than two faces coincide\n";
            }
            else
            {
                ++missingPeriodicPartnerCount;
                //std::cout << " Less than two faces coincide\n";
            }
            //std::cout << " Current ZoneIndex  = " << grid->id << std::endl;
            //std::cout << " nIZone = " << nIZone << std::endl;
            //std::cout << " LocalInterface Index = " << iIFace << " nIFaces = " << nIFaces << std::endl;
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
                break;
            }
        }

        //if ( ! flag )
        //{
        //    std::cout << "LocalInterface Index = " << iIFace << " There is a problem in the input grid. Please check it carefully!\n";
        //}
    }
    std::cout << " Periodic boundary faces missing a partner = "
              << missingPeriodicPartnerCount << "\n";

}

void IFaceLink::MatchPeriodicInterface( Grid * grid )
{
    InterFace * interFace = grid->interFace.get();
    if ( ! interFace ) return;

    int nIFaces = this->l2g[ grid->id ].size();

    for ( int iIFace = 0; iIFace < nIFaces; ++ iIFace )
    {
        int gIFace = this->l2g[ grid->id ][ iIFace ];
        int nIZone = this->gI2Zid[ gIFace ].size();

        if (nIZone == 2) continue;

        // faceArray now stores IntField directly
        const IntField & nodeId = this->face_search->faceArray[ gIFace ];

        RealField xList, yList, zList;
        this->point_search->GetFaceCoorList( nodeId, xList, yList, zList );

        RealField xxList, yyList, zzList;
        f2fmap.FindFace( xList, yList, zList, xxList, yyList, zzList );

        int nNodes = xxList.size( );
        IntField faceNode_period;

        for ( int i = 0; i < nNodes; ++ i )
        {
            Real xm = xxList[ i ];
            Real ym = yyList[ i ];
            Real zm = zzList[ i ];
            int id = this->point_search->FindPoint( xm, ym, zm );
            if ( id == INVALID_INDEX )
            {
                faceNode_period.clear();
                break;
            }
            faceNode_period.push_back( id );
        }

        if ( faceNode_period.size() != nNodes ) continue;

        int faceId_period = this->face_search->FindFace( faceNode_period );
        if ( faceId_period == INVALID_INDEX ) continue;

        const IntField & periodicZones = this->gI2Zid[ faceId_period ];
        const IntField & periodicLocalIds = this->g2l[ faceId_period ];
        if ( periodicZones.empty() || periodicLocalIds.empty() ) continue;

        int nZid_period = periodicZones[ 0 ];
        int lId_period  = periodicLocalIds[ 0 ];

        interFace->zoneId[ iIFace ] = nZid_period;
        interFace->localInterfaceId[ iIFace ] = lId_period;

    }

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
