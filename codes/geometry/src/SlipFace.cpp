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

#include "SlipFace.h"
#include "Grid.h"
#include "UnsGrid.h"
#include "FaceMesh.h"
#include "Zone.h"
#include "ZoneState.h"
#include "FaceTopo.h"
#include "BcRecord.h"
#include "Boundary.h"
#include "DataStorage.h"
#include "Parallel.h"
#include "DataBook.h"
#include "SolverDef.h"
#include "LogFile.h"
#include "Parallel.h"
#include "HXMath.h"
#include "DataStorage.h"
#include <algorithm>
using namespace std;

BeginNameSpace( ONEFLOW )

SlipFace::SlipFace()
{
    this->nSlipFace = 0;
}

SlipFace::~SlipFace()
{
}

void SlipFace::Set( int nSlipFace, Grid * parent )
{
    this->Resize( nSlipFace );
    this->parent = parent;
    if ( this->parent )
    {
        this->zoneid = parent->id;
    }
}

void SlipFace::Resize( int nSlipFace )
{
    if ( nSlipFace <= 0 ) nSlipFace = 0;
    this->nSlipFace = nSlipFace;
    s2b.resize( nSlipFace );
    bcIdList.resize( nSlipFace );
    xfcList.resize( nSlipFace );
    yfcList.resize( nSlipFace );
    zfcList.resize( nSlipFace );

    zidList.resize( nSlipFace );
    tslipList.resize( nSlipFace );
    distList.resize( nSlipFace, LARGE );
}

void SlipFace::InitDist()
{
    int nSlipFace = distList.size();
    zidList = -1;
    tslipList = -1;
    distList = LARGE;
}

void SlipFace::Init()
{
    Grid * grid = Zone::GetGrid();
    IntField bcIdList;
}

SlipfacePair * SlipFace::GetSlipfacePair( int iNei )
{
    return this->slipfacePairs[ iNei ];
}

void SlipFace::InitNeighborZoneInfo()
{
    int nZone = ZoneState::nZones;

    IntField flags;
    this->InitNeighborFlag( flags );

    //The adjacent blocks of this block are calculated
    this->AllocateNeighbor();

    int iNei = 0;
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        if ( flags[ iZone ] )
        {
            this->InitNeighborZoneInfo( iNei, iZone );

            ++ iNei;
        }
    }
}

void SlipFace::InitNeighborZoneInfo( int iNei, int iZone )
{
    SlipfacePair * slipfacePair = slipfacePairs[ iNei ];
    slipfacePair->zid = this->zoneid;
    slipfacePair->nzid = iZone;

    this->z2n.insert( pair< int, int >( iZone, iNei ) );

    this->FillRecvId( iNei );
}

void SlipFace::InitNeighborFlag( IntField & flags )
{
    flags.resize( ZoneState::nZones, 0 );

    int nSlipFace = zidList.size();

    for ( int iFace = 0; iFace < this->nSlipFace; ++ iFace )
    {
        flags[ this->zidList[ iFace ] ] = 1;
    }

    nNeighbor = 0;
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        nNeighbor += flags[ iZone ];
    }
    int kkk = 1;
}

void SlipFace::AllocateNeighbor()
{
    this->slipfacePairs.resize( nNeighbor );

    for ( int iNei = 0; iNei < nNeighbor; ++ iNei )
    {
        this->slipfacePairs[ iNei ] = new SlipfacePair();
    }
}

void SlipFace::FillRecvId( int iNei )
{
    SlipfacePair * slipfacePair = slipfacePairs[ iNei ];
    slipfacePair->idrecv.resize( 0 );

    for ( int iFace = 0; iFace < this->nSlipFace; ++ iFace )
    {
        if ( this->zidList[ iFace ] == slipfacePair->nzid )
        {
            slipfacePair->idrecv.push_back( iFace );
        }
    }
}

void SlipFace::CalcSendId( int iNei, IntField & idsend )
{
    SlipfacePair * slipfacePair = slipfacePairs[ iNei ];
    idsend.resize( 0 );

    for ( int iFace = 0; iFace < this->nSlipFace; ++ iFace )
    {
        if ( this->zidList[ iFace ] == slipfacePair->nzid )
        {
            idsend.push_back( this->tslipList[ iFace ] );
        }
    }
}

void SlipFace::SetSendId( int zid, IntField & idsend )
{
    int iNei = this->z2n[ zid ];
    SlipfacePair * slipfacePair = slipfacePairs[ iNei ];
    slipfacePair->idsend = idsend;
}

LocalSlipFace::LocalSlipFace()
{
    ;
}

LocalSlipFace::~LocalSlipFace()
{
    ;
}

void LocalSlipFace::AddSlipFace( SlipFace * slipFace )
{
    this->data.push_back( slipFace );
}

void LocalSlipFace::PackData()
{
    ;
}

void LocalSlipFace::InitDist()
{
    for ( int i = 0; i < data.size(); ++ i )
    {
        SlipFace * slipface = data[ i ];
        slipface->InitDist();
    }
}

GlobalSlipFace::GlobalSlipFace()
{
    ;
}

GlobalSlipFace::~GlobalSlipFace()
{
    for ( int i = 0; i < this->data.size(); ++ i )
    {
        delete this->data[ i ];
    }
}

void GlobalSlipFace::AddSlipFace( SlipFace * slipFace )
{
    this->data.push_back( slipFace );
}

void GlobalSlipFace::Swap()
{
    for ( int proc = 0; proc < Parallel::nProc; ++ proc )
    {
        DataBook * dataBook = new DataBook();
        if ( proc == Parallel::pid )
        {
            Init( dataBook );
        }

        HXBcast( dataBook, proc );
        Trans( dataBook );
        delete dataBook;
    }
}

void GlobalSlipFace::Init( DataBook * dataBook )
{
    int nSize = localSlipFace->data.size();
    HXWrite( dataBook, nSize );
    for ( int i = 0; i < nSize; ++ i )
    {
        SlipFace * slipFace = localSlipFace->data[ i ];
        int nSlip = slipFace->xfcList.size();
        HXWrite( dataBook, nSlip );
        HXWrite( dataBook, slipFace->zoneid );
        HXWrite( dataBook, slipFace->xfcList );
        HXWrite( dataBook, slipFace->yfcList );
        HXWrite( dataBook, slipFace->zfcList );
    }
}

void GlobalSlipFace::Trans( DataBook * dataBook )
{
    dataBook->MoveToBegin();
    int nSize = -1;
    HXRead( dataBook, nSize );
    for ( int i = 0; i < nSize; ++ i )
    {
        SlipFace * slipFace = new SlipFace();
        this->AddSlipFace( slipFace );
        int nSlip = -1;
        HXRead( dataBook, nSlip );
        slipFace->Set( nSlip, 0 );
        HXRead( dataBook, slipFace->zoneid );
        HXRead( dataBook, slipFace->xfcList );
        HXRead( dataBook, slipFace->yfcList );
        HXRead( dataBook, slipFace->zfcList );
    }
}

void GlobalSlipFace::CalcDist()
{
    localSlipFace->InitDist();
    for ( int i = 0; i < localSlipFace->data.size(); ++ i )
    {
        SlipFace * slipface = localSlipFace->data[ i ];
        this->CalcDist( slipface );
    }
}

void GlobalSlipFace::CalcDist( SlipFace * slipface )
{
    for ( int iSlip = 0; iSlip < slipface->nSlipFace; ++ iSlip )
    {
        Real xfc = slipface->xfcList[ iSlip ];
        Real yfc = slipface->yfcList[ iSlip ];
        Real zfc = slipface->zfcList[ iSlip ];
        Real dst = slipface->distList[ iSlip ];
        int zid = slipface->zidList[ iSlip ];
        int isbc = slipface->tslipList[ iSlip ];

        for ( int i = 0; i < data.size(); ++ i )
        {
            SlipFace * slipface1 = this->data[ i ];
            if ( slipface1->zoneid == slipface->zoneid ) continue;
            this->Calc( xfc, yfc, zfc, dst, zid, isbc, slipface1 );
        }

        slipface->distList[ iSlip ] = dst;
        slipface->tslipList[ iSlip ] = isbc;
        slipface->zidList[ iSlip ] = zid;
    }
}

void GlobalSlipFace::Calc( Real xfc, Real yfc, Real zfc, Real & dst, int & zid, int & isbc, SlipFace * slipface )
{
    for ( int iSlip = 0; iSlip < slipface->nSlipFace; ++ iSlip )
    {
        Real xfc1 = slipface->xfcList[ iSlip ];
        Real yfc1 = slipface->yfcList[ iSlip ];
        Real zfc1 = slipface->zfcList[ iSlip ];
        Real dx = xfc1 - xfc;
        Real dy = yfc1 - yfc;
        Real dz = zfc1 - zfc;
        Real ds = SQR( dx, dy, dz );
        if ( ds < dst )
        {
            dst = ds;
            zid = slipface->zoneid;
            isbc = iSlip;
        }
    }
}

LocalSlipFace * localSlipFace;
GlobalSlipFace * globalSlipFace;
SlipFaceTopo * slipFaceTopo;
void CreateSlip();
void FreeSlip();

void CreateSlip()
{
    localSlipFace = new LocalSlipFace();
    globalSlipFace = new GlobalSlipFace();
    slipFaceTopo = new SlipFaceTopo();
}

void FreeSlip()
{
    delete localSlipFace;
    delete globalSlipFace;
    delete slipFaceTopo;
}

void InitSlipFaceTopo()
{
    CreateSlip();

    int nZone = ZoneState::nZones;

    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        if ( ! ZoneState::IsValidZone( iZone ) ) continue;

        ZoneState::zid = iZone;

        UnsGrid * grid = Zone::GetUnsGrid();

        BcRecord * bcRecord = grid->faceTopo->bcManager->bcRecord;
        int nBFace = bcRecord->GetNBFace();

        int nSlipFace = 0;
        for ( int iBFace = 0; iBFace < nBFace; ++ iBFace )
        {
            int bcType = bcRecord->bcType[ iBFace ];
            if ( BC::IsSlipfaceBc( bcType ) )
            {
                ++ nSlipFace;
            }
        }

        SlipFace * slipFace = grid->slipFace;
        slipFace->Set( nSlipFace, grid );

        localSlipFace->AddSlipFace( slipFace );

        int iFace = 0;
        for ( int iBFace = 0; iBFace < nBFace; ++ iBFace )
        {
            int bcType = bcRecord->bcType[ iBFace ];
            if ( ! BC::IsSlipfaceBc( bcType ) )
            {
                continue;
            }

            slipFace->s2b[ iFace ++ ] = iBFace;
        }

        RealField & xfc = grid->faceMesh->xfc;
        RealField & yfc = grid->faceMesh->yfc;
        RealField & zfc = grid->faceMesh->zfc;

        for ( int iSlip = 0; iSlip < nSlipFace; ++ iSlip )
        {
            int iBFace = slipFace->s2b[ iSlip ];
            slipFace->xfcList[ iSlip ] = xfc[ iBFace ];
            slipFace->yfcList[ iSlip ] = yfc[ iBFace ];
            slipFace->zfcList[ iSlip ] = zfc[ iBFace ];

        }
        int kkk = 1;
    }

    globalSlipFace->Swap();
    globalSlipFace->CalcDist();

    slipFaceTopo->InitZoneNeighborsInfo();
    slipFaceTopo->SwapNeighborZoneInfo();

    FreeSlip();

    int kkk = 1;
}

SlipfacePair::SlipfacePair()
{
    ;
}

SlipfacePair::~SlipfacePair()
{
    ;
}

SlipFaceTopo::SlipFaceTopo()
{
    ;
}

SlipFaceTopo::~SlipFaceTopo()
{
}

void SlipFaceTopo::InitZoneNeighborsInfo()
{
    int nZone = ZoneState::nZones;

    this->data.resize( nZone );

    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        if ( ! ZoneState::IsValidZone( iZone ) ) continue;

        Grid * grid = Zone::GetGrid( iZone );

        grid->slipFace->InitNeighborZoneInfo();
    }

    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        if ( ! ZoneState::IsValidZone( iZone ) ) continue;

        Grid * grid = Zone::GetGrid( iZone );

        IntField & t = this->data[ iZone ];

        for ( int iNei = 0; iNei < grid->slipFace->nNeighbor; ++ iNei )
        {
            SlipfacePair * slipfacePair = grid->slipFace->GetSlipfacePair( iNei );

            t.push_back( slipfacePair->nzid );
        }
        std::sort( t.begin(), t.end() );
    }
}

void SlipFaceTopo::SwapNeighborZoneInfo()
{
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        int pid = ZoneState::pid[ iZone ];

        IntField & t = this->data[ iZone ];
        int nNeighbor = t.size();

        ONEFLOW::HXBcast( & nNeighbor, 1, pid );

        t.resize( nNeighbor );

        if ( t.size() == 0 ) continue;

        ONEFLOW::HXBcast( & t[ 0 ], nNeighbor, pid );
    }
    this->SwapNeighborsSendContent();
}

void SlipFaceTopo::SwapNeighborsSendContent()
{
    int gl = 0;
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        int spid = ZoneState::pid[ iZone ];

        IntField & t = this->data[ iZone ];
        int nNeighbor = t.size();

        for ( int iNei = 0; iNei < nNeighbor; ++ iNei )
        {
            int nZid = t[ iNei ];
            int rpid = ZoneState::pid[ nZid ];
            int nSend = 0;
            IntField idsend;

            if ( Parallel::pid == spid )
            {
                Grid * grid = Zone::GetGrid( iZone );
                SlipfacePair * slipfacePair = grid->slipFace->GetSlipfacePair( iNei );

                grid->slipFace->CalcSendId( iNei, idsend );

                nSend = idsend.size();
            }

            ONEFLOW::HXSwapData( & nSend, 1, spid, rpid, iZone + gl * ZoneState::nZones );

            if ( Parallel::pid  == rpid && 
                          spid  != rpid )
            {
                idsend.resize( nSend );
            }

            ONEFLOW::HXSwapData( & idsend[ 0 ], nSend, spid, rpid, iZone + gl * ZoneState::nZones );

            if ( Parallel::pid == rpid )
            {
                Grid * gridN = Zone::GetGrid( nZid );
                gridN->slipFace->SetSendId( iZone, idsend );
            }
        }
    }
}


EndNameSpace