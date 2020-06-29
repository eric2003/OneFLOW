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

#include "BcRecord.h"
#include "Boundary.h"
#include "Dimension.h"
#include "InterFace.h"
#include "IFaceLink.h"
#include "FaceSearch.h"
#include "HXMath.h"
#include "HXStd.h"

BeginNameSpace( ONEFLOW )

BcInfo::BcInfo()
{
    ;
}

BcInfo::~BcInfo()
{
    ;
}

BcRecord::BcRecord()
{
    bcInfo = 0;
}

BcRecord::~BcRecord()
{
    delete bcInfo;
}

void BcRecord::CreateBcRegion()
{
    if ( bcInfo ) return;
    this->bcInfo = new BcInfo();

    IntSet bcTypeSet;
    IntSet bcUserTypeSet;

    int nBFace = this->GetNBFace();

    for ( int iFace = 0; iFace < nBFace; ++ iFace )
    {
        int bcType = this->bcType[ iFace ];
        int bcRegion = this->bcRegion[ iFace ];
        bcTypeSet.insert( bcType );
        if ( bcType == BC::GENERIC_2 )
        {
            bcUserTypeSet.insert( bcRegion );
        }
    }

    ONEFLOW::Set2Array( bcTypeSet, bcInfo->bcType );

    int nRegion = bcInfo->bcType.size();
    bcInfo->bcFace.resize( nRegion );
    bcInfo->bcRegion.resize( nRegion );
    bcInfo->bcdtkey.resize( nRegion );

    for ( int ir = 0; ir < nRegion; ++ ir )
    {
        int targetBcType = bcInfo->bcType[ ir ];
        for ( int iFace = 0; iFace < nBFace; ++ iFace )
        {
            int bcType = this->bcType[ iFace ];
            int bcdtkey = this->bcdtkey[ iFace ];
            if ( bcType == targetBcType )
            {
                int bcRegion = this->bcRegion[ iFace ];
                bcInfo->bcFace[ ir ].push_back( iFace );
                bcInfo->bcRegion[ ir ].push_back( bcRegion );
                bcInfo->bcdtkey[ ir ].push_back( bcdtkey );
            }
        }
    }
    int kkk = 1;
}

int BcRecord::GetNBFace()
{
    return bcType.size();
}

void BcRecord::Init( UInt nBFace )
{
    this->bcType.resize( nBFace );
    this->bcdtkey.resize( nBFace );
    this->bcRegion.resize( nBFace );
}

int BcRecord::CalcNIFace()
{
    int nBFace = this->GetNBFace();

    int nIFace = 0;
    for ( int iFace = 0; iFace < nBFace; ++ iFace )
    {
        if ( BC::IsInterfaceBc( this->bcType[ iFace ] ) )
        {
            ++ nIFace;
        }
    }

    return nIFace;
}

int BcRecord::CalcNumWallFace()
{
    int nBFace = this->GetNBFace();

    int nWFace = 0;
    for ( int iFace = 0; iFace < nBFace; ++ iFace )
    {
        if ( BC::IsWallBc( this->bcType[ iFace ] ) )
        {
            ++ nWFace;
        }
    }

    return nWFace;
}

void BcRecord::CreateI2B( InterFace * interFace )
{
    if ( ! interFace ) return;

    int nBFace = this->GetNBFace();

    int iFace  = 0;
    for ( int iBFace = 0; iBFace < nBFace; ++ iBFace )
    {
        if ( ! BC::IsInterfaceBc( this->bcType[ iBFace ] ) )
        {
            continue;
        }

        interFace->i2b[ iFace ] = iBFace;
        ++ iFace;
    }
}

BcManager::BcManager()
{
    bcRecord = new BcRecord();
    bcRecordNew = new BcRecord();
}

BcManager::~BcManager()
{
    delete bcRecord;
    delete bcRecordNew;
}

void BcManager::PreProcess()
{
    int nBFace = this->bcRecord->GetNBFace();
    bcFlag.resize( nBFace, 1 );
}

bool BcManager::ExistInterface()
{
    int nBFace = this->bcRecord->GetNBFace();

    for ( int iFace = 0; iFace < nBFace; ++ iFace )
    {
        int bcType = this->bcRecord->bcType[ iFace ];
        if ( BC::IsInterfaceBc( bcType ) )
        {
            return true;
        }
    }
    return false;
}

void BcManager::Update()
{
    * this->bcRecord = * this->bcRecordNew;
    int nBFace = this->bcRecord->bcType.size();
    this->bcFlag.resize( nBFace );
    this->bcFlag = 1;
}

void BcManager::CalcBcType( IntField & bcTypeList )
{
    IntSet bcTypeSet;

    int nBFace = this->bcRecord->GetNBFace();

    for ( int iFace = 0; iFace < nBFace; ++ iFace )
    {
        int bcType = this->bcRecord->bcType[ iFace ];
        bcTypeSet.insert( bcType );
    }

    for ( IntSet::iterator iter = bcTypeSet.begin(); iter != bcTypeSet.end(); ++ iter )
    {
        bcTypeList.push_back( * iter );
    }
}

void BasicRegion::SetRegion( int ist, int ied, int jst, int jed )
{
    this->start[ 0 ] = ist;
    this->end  [ 0 ] = ied;
    this->start[ 1 ] = jst;
    this->end  [ 1 ] = jed;
    this->start[ 2 ] = 1;
    this->end  [ 2 ] = 1;
}

void BasicRegion::SetRegion( int ist, int ied, int jst, int jed, int kst, int ked )
{
    this->start[ 0 ] = ist;
    this->end  [ 0 ] = ied;
    this->start[ 1 ] = jst;
    this->end  [ 1 ] = jed;
    this->start[ 2 ] = kst;
    this->end  [ 2 ] = ked;
}

TestRegion::TestRegion()
{
    ;
}

TestRegion::~TestRegion()
{
    ;
}

void TestRegion::Run( BasicRegion * r, int dimension ) 
{
    int n = dimension;
    for ( int i = 0; i < n; ++ i )
    {
        p1[ i ] = r->start[ i ];
        p2[ i ] = r->end[ i ];
    }

    for ( int i = 0; i < n; ++ i )
    {
        if ( p1[ i ] == p2[ i ] )
        {
            a[ i ] = 0;
        }
        else
        {
            if ( dimension == THREE_D )
            {
                if ( p1[ i ] > 0 )
                {
                    a[ i ] = 1;
                }
                else
                {
                    a[ i ] = 2;
                }
            }
            else
            {
                a[ i ] = 1;
            }
        }
    }

    for ( int i = 0; i < n; ++ i )
    {
        sign[ i ] = 1;
        if ( ABS( p1[ i ] ) > ABS( p2[ i ] ) )
        {
            sign[ i ] = - 1;
        }
    }

    for ( int m = 0; m < n; ++ m )
    {
        int aa = MIN( ABS( p1[ m ] ), ABS( p2[ m ] ) );
        int bb = MAX( ABS( p1[ m ] ), ABS( p2[ m ] ) );
        p1[ m ] = aa;
        p2[ m ] = bb;
    }
}

TestRegionM::TestRegionM()
{
    ;
}

TestRegionM::~TestRegionM()
{
    ;
}

void TestRegionM::Run( BcRegion * bcRegion, int dimension )
{
    s.Run( bcRegion->s, dimension );
    t.Run( bcRegion->t, dimension );

    for ( int i = 0; i < 3; ++ i )
    {
        itransform[ i ] = ( i + 1 );
    }

    for ( int i = 0; i < dimension; ++ i )
    {
        int v = s.a[ i ];
        for ( int j = 0; j < dimension; ++ j )
        {
            if ( s.a[ i ] == t.a[ j ] )
            {
                int sign = s.sign[ i ] * t.sign[ j ];
                itransform[ i ] = ( j + 1 ) * sign;
            }
        }
    }
}

BcRegion::BcRegion( int zid, int rid )
{
    s = new BasicRegion();
    t = new BasicRegion();

    this->rid = rid;
    s->zid = zid;
}

BcRegion::~BcRegion()
{
    delete s;
    delete t;
}

void BcRegion::GetNormalizeIJKRegion( int & ist, int & ied, int & jst, int & jed, int & kst, int & ked )
{
    ist = MIN( ABS( this->s->start[ 0 ] ), ABS( this->s->end[ 0 ] ) );
    ied = MAX( ABS( this->s->start[ 0 ] ), ABS( this->s->end[ 0 ] ) );
    jst = MIN( ABS( this->s->start[ 1 ] ), ABS( this->s->end[ 1 ] ) );
    jed = MAX( ABS( this->s->start[ 1 ] ), ABS( this->s->end[ 1 ] ) );
    kst = MIN( ABS( this->s->start[ 2 ] ), ABS( this->s->end[ 2 ] ) );
    ked = MAX( ABS( this->s->start[ 2 ] ), ABS( this->s->end[ 2 ] ) );
}

int BcRegion::CalcRegionCells()
{
    int imin, imax, jmin, jmax, kmin, kmax;
    this->GetNormalizeIJKRegion( imin, imax, jmin, jmax, kmin, kmax );

    int di = ONEFLOW::MAX( imax - imin, 1 );
    int dj = ONEFLOW::MAX( jmax - jmin, 1 );
    int dk = ONEFLOW::MAX( kmax - kmin, 1 );

    int nRegionCells = di * dj * dk;

    return nRegionCells;
}

BcRegionGroup::BcRegionGroup()
{
    regions = 0;
}

BcRegionGroup::~BcRegionGroup()
{
    if ( regions )
    {
        int nBcRegions = regions->size();
        for ( int ir = 0; ir < nBcRegions; ++ ir )
        {
            delete ( * regions )[ ir ];
        }
        delete regions;
    }
}

void BcRegionGroup::Create( int nBcRegions )
{
    regions = new HXVector< BcRegion * >( nBcRegions );
}

void BcRegionGroup::SetBcRegion( int ir, BcRegion * bcRegion )
{
    ( * regions )[ ir ] = bcRegion;
}

BcRegion *  BcRegionGroup::GetBcRegion( int ir )
{
    return ( * regions )[ ir ];
}

EndNameSpace