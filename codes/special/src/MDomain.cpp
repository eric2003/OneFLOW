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

#include "BlockFaceSolver.h"
#include "MLine.h"
#include "SimpleDomain.h"
#include "LineMachine.h"
#include "DomainMachine.h"
#include "MDomain.h"
#include "SDomain.h"
#include "HXPointer.h"
#include "CurveInfo.h"
#include "SegmentCtrl.h"
#include "BlockElem.h"
#include "HXCgns.h"
#include "Dimension.h"
#include <algorithm>
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

MDomain::MDomain()
{
    ;
}

MDomain::~MDomain()
{
    int nSDomain = sDomainList.size();
    for ( int iSDomain = 0; iSDomain < nSDomain; ++ iSDomain )
    {
        delete sDomainList[ iSDomain ];
    }
}

void MDomain::CalcSubDomainCtrlCoor()
{
    while ( true )
    {
        bool local_flag = true;
        int nDomain = sDomainList.size();
        for ( int iDomain = 0; iDomain < nDomain; ++ iDomain )
        {
            SDomain * sDomain = sDomainList[ iDomain ];
            bool flag = sDomain->CalcSingleDomainCoor();
            local_flag = ( flag && local_flag );
        }
        if ( local_flag ) break;
    };
}

void MDomain::AddSubDomain( int fid, IntField & lineList, IntField & posList )
{
    SDomain * sdomain = new SDomain();
    sdomain->coorMap = this->coorMap;
    sdomain->SetDomain( fid, lineList, posList );
    sDomainList.push_back( sdomain );
}

SDomain * MDomain::FindSDomain( int fid )
{
    int nSDomain = sDomainList.size();
    for ( int iSDomain = 0; iSDomain < nSDomain; ++ iSDomain )
    {
        SDomain * sDomain = sDomainList[ iSDomain ];
        if ( sDomain->domain_id == fid )
        {
            return sDomain;
        }
    }
    return 0;
}

void MDomain::CalcCoor()
{
    int closedLine = 1;
    this->CalcBcCoor( this->coorMap, closedLine );
    this->CalcSubDomainCtrlCoor();
}

int MDomain::GetNsubDomain()
{
    return this->sDomainList.size();
}

void MDomain::ConstructMultiDomainTopo()
{
    int nSDomain = sDomainList.size();
    for ( int iSDomain = 0; iSDomain < nSDomain; ++ iSDomain )
    {
        SDomain * sDomain = sDomainList[ iSDomain ];
        sDomain->ConstructSDomainCtrlPoint();
        sDomain->ConstructDomainTopo();
    }

    this->ConstructMultiLineToDomainMap();
    this->ConstructMultiPointToDomainMap();
    this->ConstructMultiPointToPointMap();
    this->ConstructBcPoint();
    this->ConstructCtrlPoint();
}

void MDomain::ConstructMultiLineToDomainMap()
{
    for ( int iDomain = 0; iDomain < this->sDomainList.size(); ++ iDomain )
    {
        SDomain * sdomain = this->sDomainList[ iDomain ];
        sdomain->ConstructLineToDomainMap( this->lineToDomainMap );
    }

    int kkk = 1;
}

void MDomain::ConstructMultiPointToDomainMap()
{
    for ( int iDomain = 0; iDomain < this->sDomainList.size(); ++ iDomain )
    {
        SDomain * sdomain = this->sDomainList[ iDomain ];
        sdomain->ConstructPointToDomainMap( this->pointToDomainMap );
    }
}

void MDomain::ConstructMultiPointToPointMap()
{
    for ( int iDomain = 0; iDomain < this->sDomainList.size(); ++ iDomain )
    {
        SDomain * sdomain = this->sDomainList[ iDomain ];
        sdomain->ConstructPointToPointMap( this->pointToPointMap );
    }
}

void MDomain::CreateInpFaceList( HXVector< Face2D * > &facelist )
{
    for ( int iDomain = 0; iDomain < this->sDomainList.size(); ++ iDomain )
    {
        SDomain * sDomain = this->sDomainList[ iDomain ];
        Face2D * face2d = new Face2D();
        face2d->face_id = sDomain->domain_id;
        face2d->ctrlpoints = sDomain->ctrlpoints;
        BlkF2C & face_struct = blkFaceSolver.myFaceSolver.face2Block[ face2d->face_id ];
        face2d->bcType = face_struct.bctype;
        face2d->CalcStEd( coorMap );
        facelist.push_back( face2d );
    }
}

void MDomain::SetBlkBcMesh( Block3D * blk3d )
{
    int nSDomain = sDomainList.size();
    for ( int iSDomain = 0; iSDomain < nSDomain; ++ iSDomain )
    {
        SDomain * sDomain = sDomainList[ iSDomain ];
        sDomain->SetBlkBcMesh( blk3d );
    }
}

EndNameSpace