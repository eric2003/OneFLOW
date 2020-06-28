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
class SlipfacePair;

class SlipFace
{
public:
    SlipFace();
    ~SlipFace();
public:
    Grid * parent;
    int zoneid;
    int nSlipFace;
    IntField s2b;
    IntField bcIdList;
    RealField xfcList;
    RealField yfcList;
    RealField zfcList;
    IntField zidList;
    IntField tslipList;
    RealField distList;
public:
    void Set( int nSlipFace, Grid * parent );
    void Resize( int nSlipFace );
    void InitDist();
    void Init();
public:
    map< int, int > z2n;
    int   nNeighbor;     //no of neighbors
    HXVector< SlipfacePair * > slipfacePairs;
    SlipfacePair * GetSlipfacePair( int iNei );
    void InitNeighborZoneInfo();
    void InitNeighborZoneInfo( int iNei, int iZone );
    void InitNeighborFlag( IntField & flags );
    void AllocateNeighbor();
    void FillRecvId( int iNei );
    void CalcSendId( int iNei, IntField & idsend );
    void SetSendId( int zid, IntField & idsend );
};

class LocalSlipFace
{
public:
    LocalSlipFace();
    ~LocalSlipFace();
public:
    HXVector< SlipFace * > data;
    void AddSlipFace( SlipFace * slipFace );
    void PackData();
    void InitDist();
};

class DataBook;

class GlobalSlipFace
{
public:
    GlobalSlipFace();
    ~GlobalSlipFace();
public:
    HXVector< SlipFace * > data;
public:
    void AddSlipFace( SlipFace * slipFace );
    void Swap();
    void Init( DataBook * dataBook );
    void Trans( DataBook * dataBook );
    void CalcDist();
    void CalcDist( SlipFace * slipface );
    void Calc( Real xfc, Real yfc, Real zfc, Real & dst, int & zid, int & isbc, SlipFace * slipface );
};


extern LocalSlipFace  * localSlipFace;
extern GlobalSlipFace * globalSlipFace;

void InitSlipFaceTopo();

class DataStorage;

class SlipfacePair
{
public:
    SlipfacePair();
    ~SlipfacePair();
public:
    //zid, nzid is the block number of the neighboring block (neighbor)
    int zid, nzid;
    IntField idrecv; // using in receiving, interface number in the current zone
    IntField idsend; // using in sending,   interface number in the tagret  zone
    int GetNSend() { return static_cast<int> (idsend.size()); };
    int GetNRecv() { return static_cast<int> (idrecv.size()); };
protected:
    DataStorage * dataSend;
    DataStorage * dataRecv;
};


class SlipFaceTopo
{
public:
    SlipFaceTopo ();
    ~SlipFaceTopo();
public:
    LinkField data;
public:
    void InitZoneNeighborsInfo();
    void SwapNeighborsSendContent();
    void SwapNeighborZoneInfo();
};

EndNameSpace