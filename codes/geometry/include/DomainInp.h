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
#include <fstream>
using namespace std;

BeginNameSpace( ONEFLOW )

class GridMediator;

class PBlk
{
public:
    PBlk();
    ~PBlk();
public:
    int blk, fid;
    int i, j, k;
    bool operator < ( const PBlk * rhs ) const;
};

class ComparePBlk
{
public:
    bool operator()( const PBlk * lhs, const PBlk * rhs ) const 
    {
        return lhs->operator< ( rhs );
    }
};

class PatchBox
{
public:
    PatchBox();
    ~PatchBox();
public:
    int imin, imax;
    int jmin, jmax;
    int kmin, kmax;
    IntField ilist, jlist, klist;
    IntField idlist;
    IntField nii, njj;
    IntField idmap;
    IntField smin, smax;
    int a, va, b, vb, c;
public:
    void Init();
    void Set( int imin, int imax, int jmin, int jmax, int kmin, int kmax );
    void GenList();
    void CalcNormal();
    void CalcIdMap( PatchBox * box );
};

class MultiDomain
{
public:
    MultiDomain();
    ~MultiDomain();
public:
    HXVector< PatchBox * > boxlist1;
    HXVector< PatchBox * > boxlist2;
    IntField zoneid1, zoneid2;
    IntField fid1, fid2;
public:
    void Add( int zid1, int fid1, int zid2, int fid2, PatchBox * box1, PatchBox * box2 );

};

class PBlkSet
{
public:
    PBlkSet();
    ~PBlkSet();
public:
    IntField id;
    HXVector< set< PBlk *, ComparePBlk > * > pinfo;
public:
    void ReSize( int nSize );
    void Add( int idx, PBlk * pblk );
    void Analysys();
    void CalcDomainPatch( int nZone, MultiDomain & md );
    void CalcDomainPatch( int iZone, int jZone, MultiDomain & md );
    bool CrossDomain( int iZone, int idomain, int jZone, int jdomain, PatchBox & box_i, PatchBox & box_j );
    bool BlkDomainInSet( int zoneid, int idomain, set< PBlk *, ComparePBlk > * pset, PBlk *& pblk );
    bool CheckPBlkList( HXVector< PBlk * > & pblk_list, PatchBox & box );

};

class StrGrid;
class IjkBox
{
public:
    IjkBox(){};
    ~IjkBox(){};
public:
    IntField imin, imax, jmin, jmax, kmin, kmax;
    IntField idir;
    IntField valid;
public:
    void CreateBox( StrGrid * grid );
    void CreateBox( int ni, int nj, int nk );
    void Add( int iminIn, int imaxIn, int jminIn, int jmaxIn, int kminIn, int kmaxIn, int idirIn );
};

class IjkSlice
{
public:
    IjkSlice();
    ~IjkSlice();
public:
    int imin, imax, jmin, jmax, kmin, kmax;
    int idir;
public:
    void CalcIdir();
    void Run();
    void RunI();
};

class BcRegionGroup;

class StrGrid;
class PointSearch;

class DomainInp
{
public:
    DomainInp();
    ~DomainInp();
public:
    void Run();
public:
    void GeneInp();
    void OutputInp( GridMediator * gridMediator );
    void CalcFacePoint( StrGrid * grid, PointSearch * pointSearch, IjkBox * ijkBox, int zId, PBlkSet * pblkSet );
    void CalcDomainPatch( int nZone, GridMediator * gridMediator );
    void CalcDomainPatch( int iZone, int jZone, GridMediator * gridMediator );
    void GetId( int zid, int i, int j, int k, int & id, GridMediator * gridMediator, PointSearch * pointSearch );
    void DumpCoor( int zid, int i, int j, int k, GridMediator * gridMediator, fstream & file );
    void Dump( MultiDomain * md, GridMediator * gridMediator, PointSearch * pointSearch );
    void FindPhysicalPatch( StrGrid * grid, MultiDomain * md, int zid, IjkBox * ijkBox );
};


EndNameSpace