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


#pragma once
#include "HXDefine.h"
#include <vector>
#include <string>
#include <map>
using namespace std;

BeginNameSpace( ONEFLOW )

class Grid;
class DataStorage;
class InterfaceInfo;

const int MAX_GHOST_LEVELS = 2;

class InterfacePair;
class InterFace;

bool IsValid( InterFace * interFace );

class InterFace
{
public:
    InterFace();
    InterFace( int nIFace, Grid * parent = 0 );
    ~InterFace();
public:
    void Set( int nIFace, Grid * parent = 0 );
public:
    int   nIFace;
    int   nNeighbor;     //no of neighbors
    IntField zoneId;             //the ID of another zone corresponding to the interface of this block.
    IntField localInterfaceId;   //the ID of another face corresponding to the interface of this block.
    IntField localCellId;        //the ID of another cell corresponding to the interface of this block.
    IntField i2b; //the serial number correspondence in all boundaries of the interface of this block.
    IntField idir;
    Grid * parent;
    map< int, int > z2n;
public:
    HXVector< DataStorage * > dataSend;
    HXVector< DataStorage * > dataRecv;
    HXVector< InterfacePair * > interFacePairs;
public:
    void AllocSendRecv();
    void DeAllocSendRecv();
    void Resize( int nIFace );
    void InitNeighborFlag( IntField & flags );
    void InitNeighborZoneInfo();
    void InitNeighborZoneInfo( int iNei, int iZone );
    void CalcSendId( int iNei, IntField & idsend );
    void FillRecvId( int iNei );
    void SetSendId( int zid, IntField & idsend );
    void AllocateNeighbor();
    int CalcNIFace( int iNei );
    IntField & GetInterfaceId( int neiId, int iSr );
};

class InterfacePair
{
public:
    InterfacePair();
    ~InterfacePair();
public:
    //zid, nzid is the block number of the neighboring block (neighbor)
    int zid, nzid;
    //nIFaces is the number of interfaces corresponding to this block of adjacent blocks (not necessarily all)
    int nIFace;
    IntField idrecv; // using in receiving, interface number in the current zone
    IntField idsend; // using in sending,   interface number in the tagret  zone
protected:
    DataStorage * dataSend;
    DataStorage * dataRecv;
};

class InterFaceTopo
{
public:
    InterFaceTopo ();
    ~InterFaceTopo();
public:
    LinkField data;
public:
    void InitZoneNeighborsInfo();
    void SwapNeighborsSendContent();
    void SwapNeighborZoneInfo();
};

void InitInterfaceTopo();

class InterFaceState
{
public:
    InterFaceState();
    ~InterFaceState();
public:
    static InterFace * interFace;
};

extern InterFaceTopo interFaceTopo;

EndNameSpace