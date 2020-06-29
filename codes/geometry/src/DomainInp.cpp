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

#include "DomainInp.h"
#include "GridPara.h"
#include "CgnsFactory.h"
#include "GridMediator.h"
#include "Su2Grid.h"
#include "DataBase.h"
#include "ClassicGrid.h"
#include "StrGrid.h"
#include "PointSearch.h"
#include "BcRecord.h"
#include "HXMath.h"
#include "FileUtil.h"
#include "Prj.h"
#include "Partition.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

PBlk::PBlk()
{
    ;
}

PBlk::~PBlk()
{
    ;
}

bool PBlk::operator < ( const PBlk * rhs ) const
{
    if ( this->blk != rhs->blk )
    {
        return this->blk < rhs->blk;
    }

    if ( this->fid != rhs->fid )
    {
        return this->fid < rhs->fid;
    }

    if ( this->i != rhs->i )
    {
        return this->i < rhs->i;
    }

    if ( this->j != rhs->j )
    {
        return this->j < rhs->j;
    }

    return this->k < rhs->k;
}

PatchBox::PatchBox()
{
    Init();
}

PatchBox::~PatchBox()
{
    ;
}

void PatchBox::Init()
{
    idmap.resize( 4 );
    for ( int i = 0; i < idmap.size(); ++ i )
    {
        idmap[ i ] = i;
    }
    smin.resize( 3 );
    smax.resize( 3 );
}

void PatchBox::Set( int imin, int imax, int jmin, int jmax, int kmin, int kmax )
{
    this->imin = imin;
    this->imax = imax;

    this->jmin = jmin;
    this->jmax = jmax;

    this->kmin = kmin;
    this->kmax = kmax;
}

void PatchBox::GenList()
{
    int i1, j1, k1;
    int i2, j2, k2;
    int i3, j3, k3;
    int i4, j4, k4;
    if ( imin == imax )
    {
        i1 = imin;
        j1 = jmin;
        k1 = kmin;

        i2 = imin;
        j2 = jmax;
        k2 = kmin;

        i3 = imin;
        j3 = jmax;
        k3 = kmax;

        i4 = imin;
        j4 = jmin;
        k4 = kmax;
    }
    else if ( jmin == jmax )
    {
        i1 = imin;
        j1 = jmin;
        k1 = kmin;

        i2 = imin;
        j2 = jmin;
        k2 = kmax;

        i3 = imax;
        j3 = jmin;
        k3 = kmax;

        i4 = imax;
        j4 = jmin;
        k4 = kmin;
    }
    else if ( kmin == kmax )
    {
        i1 = imin;
        j1 = jmin;
        k1 = kmin;

        i2 = imax;
        j2 = jmin;
        k2 = kmin;

        i3 = imax;
        j3 = jmax;
        k3 = kmin;

        i4 = imin;
        j4 = jmax;
        k4 = kmin;
    }

    ilist.resize( 0 );
    jlist.resize( 0 );
    klist.resize( 0 );

    ilist.push_back( i1 );
    jlist.push_back( j1 );
    klist.push_back( k1 );

    ilist.push_back( i2 );
    jlist.push_back( j2 );
    klist.push_back( k2 );

    ilist.push_back( i3 );
    jlist.push_back( j3 );
    klist.push_back( k3 );

    ilist.push_back( i4 );
    jlist.push_back( j4 );
    klist.push_back( k4 );
}

void PatchBox::CalcNormal()
{
    int p1 = idmap[ 0 ];
    int p2 = idmap[ 1 ];
    int p3 = idmap[ 3 ];

    nii.resize( 3 );

    nii[ 0 ] = ilist[ p2 ] - ilist[ p1 ];
    nii[ 1 ] = jlist[ p2 ] - jlist[ p1 ];
    nii[ 2 ] = klist[ p2 ] - klist[ p1 ];

    njj.resize( 3 );

    njj[ 0 ] = ilist[ p3 ] - ilist[ p1 ];
    njj[ 1 ] = jlist[ p3 ] - jlist[ p1 ];
    njj[ 2 ] = klist[ p3 ] - klist[ p1 ];

    for ( int i = 0; i < 3; ++ i )
    {
        if ( nii[ i ] != 0 )
        {
            a = i;
            va = nii[ i ];
        }

        if ( njj[ i ] != 0 )
        {
            b = i;
            vb = njj[ i ];
        }
    }

    smin[ 0 ] = imin;
    smax[ 0 ] = imax;
    smin[ 1 ] = jmin;
    smax[ 1 ] = jmax;
    smin[ 2 ] = kmin;
    smax[ 2 ] = kmax;

    if ( va < 0 )
    {
        SWAP( smin[ a ], smax[ a ] );
    }

    if ( vb < 0 )
    {
        SWAP( smin[ b ], smax[ b ] );
    }

    smin[ b ] = - smin[ b ];
    smax[ b ] = - smax[ b ];

    cout << smin[ 0 ] << " " << smax[ 0 ] << " " ;
    cout << smin[ 1 ] << " " << smax[ 1 ] << " " ;
    cout << smin[ 2 ] << " " << smax[ 2 ] << endl << endl;

    int kkk = 1;
}

void PatchBox::CalcIdMap( PatchBox * box )
{
    for ( int i = 0; i < 4; ++ i )
    {
        int ip = box->idlist[ i ];
        int pos = -1;
        for ( int j = 0; j < 4; ++ j )
        {
            int jp = this->idlist[ j ];
            if ( ip == jp )
            {
                pos = j;
                break;
            }
        }
        idmap[ i ] = pos;
    }
    int kkk = 1;
}


MultiDomain::MultiDomain()
{
    ;
}

MultiDomain::~MultiDomain()
{
    int nSize = boxlist1.size();
    for ( int i = 0; i < nSize; ++ i )
    {
        delete this->boxlist1[ i ];
        delete this->boxlist2[ i ];
    }
}

void MultiDomain::Add( int zid1, int fid1, int zid2, int fid2, PatchBox * box1, PatchBox * box2 )
{
    this->zoneid1.push_back( zid1 );
    this->zoneid2.push_back( zid2 );
    this->fid1.push_back( fid1 );
    this->fid2.push_back( fid2 );

    PatchBox * b1 = new PatchBox();
    PatchBox * b2 = new PatchBox();

    this->boxlist1.push_back( b1 );
    this->boxlist2.push_back( b2 );

    * b1 = * box1;
    * b2 = * box2;
}

PBlkSet::PBlkSet()
{
    ;
}

PBlkSet::~PBlkSet()
{
    ;
}

void PBlkSet::ReSize( int nSize )
{
    id.resize( nSize );
    pinfo.resize( nSize );
}

void PBlkSet::Add( int idx, PBlk * pblk )
{
    int nSize = this->id.size();
    if ( idx < nSize )
    {
        set< PBlk *, ComparePBlk > * pset = pinfo[ idx ];
        set< PBlk *, ComparePBlk >::iterator iter;
        iter = pset->find( pblk );
        if ( iter == pset->end() )
        {
            pset->insert( pblk );
        }
        else
        {
            delete pblk;
        }
    }
    else
    {
        this->ReSize( nSize + 1 );
        set< PBlk *, ComparePBlk > * pset = new set< PBlk *, ComparePBlk >;
        pinfo[ idx ] = pset;
        pset->insert( pblk );
        this->id[ idx ] = idx;
    }
}

void PBlkSet::Analysys()
{
    int nSize = this->id.size();
    int maxpt = 0;
    int ip = -1;
    for ( int i = 0; i < nSize; ++ i )
    {
        int nn = pinfo[ i ]->size();
        if ( maxpt < nn )
        {
        maxpt = nn;
        ip = i;
        }
    }
    cout << " maxpt = " << maxpt << " ip = " << ip << endl;

    IntField multi_point;
    for ( int i = 0; i < nSize; ++ i )
    {
        int nn = pinfo[ i ]->size();
        if ( nn > 1 )
        {
            multi_point.push_back( i );
        }
    }
    int kkk = 1;
}

void PBlkSet::CalcDomainPatch( int nZone, MultiDomain & md )
{
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        for ( int jZone = iZone + 1; jZone < nZone; ++ jZone )
        {
            CalcDomainPatch( iZone, jZone, md );
        }
    }
}

void PBlkSet::CalcDomainPatch( int iZone, int jZone, MultiDomain & md )
{
    int n = 6;
    for ( int i = 0; i < n; ++ i )
    {
        for ( int j = 0; j < n; ++ j )
        {
            PatchBox box_i;
            PatchBox box_j;

            if ( iZone == 0 && jZone == 2 && i == 1 && j == 1 )
            {
                int kkk = 1;
            }
            bool flag = CrossDomain( iZone, i , jZone, j, box_i, box_j );
            if ( flag )
            {
                md.Add( iZone, i, jZone, j, & box_i, & box_j );
                md.Add( jZone, j, iZone, i, & box_j, & box_i );
            }
        }
    }
}

bool PBlkSet::BlkDomainInSet( int zoneid, int idomain, set< PBlk *, ComparePBlk > * pset, PBlk *& pblk )
{
    set< PBlk *, ComparePBlk >::iterator iter;

    for ( iter = pset->begin(); iter != pset->end(); ++ iter )
    {
        int blk = ( * iter )->blk;
        int fid = ( * iter )->fid;

        if ( zoneid == blk && idomain == fid )
        {
            pblk = ( * iter );
            return true;
        }
    }
    return false;
}

bool PBlkSet::CheckPBlkList( HXVector< PBlk * > & pblk_list, PatchBox & box )
{
    int nSize = pblk_list.size();
    if ( nSize <= 0 ) return false;

    int imin = pblk_list[ 0 ]->i;
    int jmin = pblk_list[ 0 ]->j;
    int kmin = pblk_list[ 0 ]->k;
    int imax = pblk_list[ 0 ]->i;
    int jmax = pblk_list[ 0 ]->j;
    int kmax = pblk_list[ 0 ]->k;
    for ( int i = 0; i < pblk_list.size(); ++ i )
    {
        int ii = pblk_list[ i ]->i;
        int jj = pblk_list[ i ]->j;
        int kk = pblk_list[ i ]->k;
        imin = MIN( ii, imin );
        jmin = MIN( jj, jmin );
        kmin = MIN( kk, kmin );

        imax = MAX( ii, imax );
        jmax = MAX( jj, jmax );
        kmax = MAX( kk, kmax );
    }


    bool flag_i = ( imin == imax );
    bool flag_j = ( jmin == jmax );
    bool flag_k = ( kmin == kmax );

    int iv = 1;
    int jv = 1;
    int kv = 1;

    if ( flag_i ) iv = 0;
    if ( flag_j ) jv = 0;
    if ( flag_k ) kv = 0;

    int ijkv = iv + jv + kv;
    if ( ijkv < 2 )
    {
        cout << " not a domain " << endl;
        return false;
    }

    bool flag = false;
    for ( int k = kmin; k <= kmax; ++ k )
    {
        for ( int j = jmin; j <= jmax; ++ j )
        {
            for ( int i = imin; i <= imax; ++ i )
            {
                flag = false;
                for ( int m = 0; m < pblk_list.size(); ++ m )
                {
                    int ii = pblk_list[ m ]->i;
                    int jj = pblk_list[ m ]->j;
                    int kk = pblk_list[ m ]->k;
                    if ( i == ii && j == jj && k == kk )
                    {
                        flag = true;
                        break;
                    }
                }
                if ( ! flag ) break;
            }
            if ( ! flag ) break;
        }
        if ( ! flag ) break;
    }

    if ( ! flag )
    {
        cout << " not a domain " << endl;
        return false;
    }

    box.Set( imin, imax, jmin, jmax, kmin, kmax );
    cout << imin << " " << imax << endl;
    cout << jmin << " " << jmax << endl;
    cout << kmin << " " << kmax << endl;
    cout << " is a domain " << endl;
    return true;
}

bool PBlkSet::CrossDomain( int iZone, int idomain, int jZone, int jdomain, PatchBox & box_i, PatchBox & box_j )
{
    IntField cross;
    int nSize = this->id.size();
    HXVector< PBlk * > pblk1_list;
    HXVector< PBlk * > pblk2_list;
    for ( int i = 0; i < nSize; ++ i )
    {
        int nn = pinfo[ i ]->size();
        set< PBlk *, ComparePBlk > * pset = pinfo[ i ];
        set< PBlk *, ComparePBlk >::iterator iter;

        PBlk * pblk1 = 0;
        PBlk * pblk2 = 0;

        bool flag1 = BlkDomainInSet( iZone, idomain, pset, pblk1 );
        bool flag2 = BlkDomainInSet( jZone, jdomain, pset, pblk2 );

        if ( flag1 && flag2 )
        {
            cross.push_back( i );
            pblk1_list.push_back( pblk1 );
            pblk2_list.push_back( pblk2 );
        }
    }
    bool valid_i = CheckPBlkList( pblk1_list, box_i );
    bool valid_j = CheckPBlkList( pblk2_list, box_j );
    cout << "zone[" << iZone << "][" << jZone << "] = " << idomain << ":" << jdomain << " num = " << cross.size() << endl;
    int kkk = 1;
    return valid_i && valid_j;
}

IjkSlice::IjkSlice()
{
    ;
}

IjkSlice::~IjkSlice()
{
    ;
}

void IjkSlice::CalcIdir()
{
    if ( imin == imax )
    {
        idir = 0;
    }
    else if ( jmin == jmax )
    {
        idir = 1;
    }
    else
    {
        idir = 2;
    }
}

void IjkSlice::RunI()
{
    int n = 4;
    IntField iv( n );
    IntField jv( n );
    IntField kv( n );
    int i = imin;
    for ( int k = kmin; k <= kmax - 1; ++ k )
    {
        for ( int j = jmin; j <= jmax - 1; ++ j )
        {
            iv[ 0 ] = i;
            jv[ 0 ] = j;
            kv[ 0 ] = k;

            iv[ 1 ] = i;
            jv[ 1 ] = j + 1;
            kv[ 1 ] = k;

            iv[ 2 ] = i;
            jv[ 2 ] = j + 1;
            kv[ 2 ] = k + 1;

            iv[ 3 ] = i;
            jv[ 3 ] = j;
            kv[ 3 ] = k + 1;
        }
    }
}

void IjkSlice::Run()
{
    if ( idir == 0 )
    {
        RunI();
    }
}

void IjkBox::CreateBox( StrGrid * grid )
{
    int ni = grid->ni;
    int nj = grid->nj;
    int nk = grid->nk;

    this->CreateBox( ni, nj, nk );
}

void IjkBox::CreateBox( int ni, int nj, int nk )
{
    this->Add( 1 , 1 , 1 , nj, 1 , nk, 0 );
    this->Add( ni, ni, 1 , nj, 1 , nk, 1 );
    this->Add( 1 , ni, 1 , 1 , 1 , nk, 2 );
    this->Add( 1 , ni, nj, nj, 1 , nk, 3 );
    this->Add( 1 , ni, 1 , nj, 1 , 1 , 4 );
    this->Add( 1 , ni, 1 , nj, nk, nk, 5 );
}

void IjkBox::Add( int iminIn, int imaxIn, int jminIn, int jmaxIn, int kminIn, int kmaxIn, int idirIn )
{
    this->imin.push_back( iminIn );
    this->imax.push_back( imaxIn );

    this->jmin.push_back( jminIn );
    this->jmax.push_back( jmaxIn );

    this->kmin.push_back( kminIn );
    this->kmax.push_back( kmaxIn );

    this->idir.push_back( idirIn );
}

DomainInp::DomainInp()
{
}

DomainInp::~DomainInp()
{
}

void DomainInp::Run()
{
    this->GeneInp();
}

void DomainInp::GeneInp()
{
    CgnsFactory * cgnsFactory = new CgnsFactory();

    GridMediator * gridMediator = new GridMediator();
    gridMediator->gridFile = ONEFLOW::GetDataValue< string >( "sourceGridFileName" );
    gridMediator->bcFile   = ONEFLOW::GetDataValue< string >( "sourceGridBcName" );
    gridMediator->gridType = "plot3d";

    gridMediator->ReadPlot3DCoor();
    this->OutputInp( gridMediator );

    delete cgnsFactory;
}

void DomainInp::GetId( int zid, int i, int j, int k, int & id, GridMediator * gridMediator, PointSearch * pointSearch )
{
    StrGrid * grid = ONEFLOW::StrGridCast( gridMediator->gridVector[ zid ] );
    Field3D & xs = * grid->strx;
    Field3D & ys = * grid->stry;
    Field3D & zs = * grid->strz;

    Real xm = xs( i, j, k );
    Real ym = ys( i, j, k );
    Real zm = zs( i, j, k );

    id = pointSearch->AddPoint( xm ,ym, zm );

    int kkk = 1;
}

void DomainInp::DumpCoor( int zid, int i, int j, int k, GridMediator * gridMediator, fstream & file )
{
    StrGrid * grid = ONEFLOW::StrGridCast( gridMediator->gridVector[ zid ] );
    Field3D & xs = * grid->strx;
    Field3D & ys = * grid->stry;
    Field3D & zs = * grid->strz;

    RealField coor( 3 );

    coor[ 0 ] = xs( i, j, k );
    coor[ 1 ] = ys( i, j, k );
    coor[ 2 ] = zs( i, j, k );

    file << coor[ 0 ] << " " << coor[ 1 ] << " " << coor[ 2 ] << endl;

}

void DomainInp::FindPhysicalPatch( StrGrid * grid, MultiDomain * md, int zid, IjkBox * ijkBox )
{
    int ni = grid->ni;
    int nj = grid->nj;
    int nk = grid->nk;

    int nSize = md->boxlist1.size();
    IntSet fidlist;
    for ( int i = 0; i < nSize; ++ i )
    {
        int zzid = md->zoneid1[ i ];
        int ffid = md->fid1[ i ];
        if ( zid != zzid ) continue;
        fidlist.insert( ffid );
    }

    ijkBox->CreateBox( grid );
    ijkBox->valid.resize( 6 );

    for ( int i = 0; i < 6; ++ i )
    {
        ijkBox->valid[ i ] = 0;
        if ( fidlist.find( i ) == fidlist.end() )
        {
            ijkBox->valid[ i ] = 1;
        }
    }

    int kkk = 1;
}

void DomainInp::Dump( MultiDomain * md, GridMediator * gridMediator, PointSearch * pointSearch )
{
    fstream file;
    string fileName = "test.inp";
    ONEFLOW::OpenPrjFile( file, fileName, ios_base::out );

    Grids grids = gridMediator->gridVector;
    int nZone = grids.size();

    int width = 5;
    int solver_id = 1;
    file << setw( width ) << solver_id << endl;
    file << setw( width ) << nZone << endl;
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        StrGrid * grid = ONEFLOW::StrGridCast( grids[ iZone ] );
        int ni = grid->ni;
        int nj = grid->nj;
        int nk = grid->nk;

        file << setw( width ) << ni;
        file << setw( width ) << nj;
        file << setw( width ) << nk << endl;

        file << "zone" << iZone + 1 << endl;

        IjkBox ijkBox;
        FindPhysicalPatch( grid, md, iZone, & ijkBox );
        int nn = ijkBox.valid.size();
        int count = 0;
        for ( int i = 0; i < nn; ++ i )
        {
            if ( ijkBox.valid[ i ] == 1 )
            {
                count ++;
            }
        }

        int nSize = md->boxlist1.size();
        for ( int i = 0; i < nSize; ++ i )
        {
            int zid1 = md->zoneid1[ i ];
            if ( zid1 == iZone )
            {
                count ++;
            }
        }

        file << setw( width ) << count << endl;

        for ( int i = 0; i < nn; ++ i )
        {
            if ( ijkBox.valid[ i ] == 1 )
            {
                int bctype = 4;
                file << setw( width ) << ijkBox.imin[ i ] << setw( width ) << ijkBox.imax[ i ];
                file << setw( width ) << ijkBox.jmin[ i ] << setw( width ) << ijkBox.jmax[ i ];
                file << setw( width ) << ijkBox.kmin[ i ] << setw( width ) << ijkBox.kmax[ i ];
                file << setw( width ) << bctype << endl;
            }
        }

        //int nSize = md->boxlist1.size();
        for ( int i = 0; i < nSize; ++ i )
        {
            PatchBox * box1 = md->boxlist1[ i ];
            PatchBox * box2 = md->boxlist2[ i ];

            int zid1 = md->zoneid1[ i ];

            if ( zid1 != iZone ) continue;

            box1->GenList();
            box1->idlist.resize( 0 );
            for ( int p = 0; p < 4; ++ p )
            {
                int ip = box1->ilist[ p ];
                int jp = box1->jlist[ p ];
                int kp = box1->klist[ p ];
                int id = -1;
                GetId( zid1, ip, jp, kp, id, gridMediator, pointSearch );
                box1->idlist.push_back( id );
                //DumpCoor( zid1, ip, jp, kp, gridMediator, file );
            }

            box1->CalcNormal();

            int zid2 = md->zoneid2[ i ];

            box2->GenList();
            box2->idlist.resize( 0 );

            for ( int p = 0; p < 4; ++ p )
            {
                int ip = box2->ilist[ p ];
                int jp = box2->jlist[ p ];
                int kp = box2->klist[ p ];
                int id = -1;
                GetId( zid2, ip, jp, kp, id, gridMediator, pointSearch );
                box2->idlist.push_back( id );
                //DumpCoor( zid2, ip, jp, kp, gridMediator, file );
            }
            box2->CalcIdMap( box1 );
            box2->CalcNormal();

            file << setw( width ) << box1->smin[ 0 ] << setw( width ) << box1->smax[ 0 ];
            file << setw( width ) << box1->smin[ 1 ] << setw( width ) << box1->smax[ 1 ];
            file << setw( width ) << box1->smin[ 2 ] << setw( width ) << box1->smax[ 2 ];
            file << setw( width ) << -1 << endl;
            file << setw( width ) << box2->smin[ 0 ] << setw( width ) << box2->smax[ 0 ];
            file << setw( width ) << box2->smin[ 1 ] << setw( width ) << box2->smax[ 1 ];
            file << setw( width ) << box2->smin[ 2 ] << setw( width ) << box2->smax[ 2 ];
            file << setw( width ) << zid2 + 1 << endl;

            int kkk = 1;
        }
    }
    ONEFLOW::CloseFile( file );
    int kkk = 1;
}

void DomainInp::OutputInp( GridMediator * gridMediator )
{
    Grids grids = gridMediator->gridVector;

    PointSearch pointSearch;
    pointSearch.Initialize( grids );
    PBlkSet pblkSet;

    int nZone = grids.size();

    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        StrGrid * grid = ONEFLOW::StrGridCast( grids[ iZone ] );
        int ni = grid->ni;
        int nj = grid->nj;
        int nk = grid->nk;

        IjkBox ijkBox;
        ijkBox.CreateBox( ni, nj, nk );

        cout << " block = " << iZone + 1 << "\n";
        this->CalcFacePoint( grid, & pointSearch, & ijkBox, iZone, & pblkSet );
    }
    pblkSet.Analysys();
    MultiDomain md;
    pblkSet.CalcDomainPatch( nZone, md );
    Dump( & md, gridMediator, & pointSearch );
    int kkk = 1;
}

void DomainInp::CalcDomainPatch( int nZone, GridMediator * gridMediator )
{
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        for ( int jZone = iZone + 1; jZone < nZone; ++ jZone )
        {
            CalcDomainPatch( iZone, jZone, gridMediator );
        }
    }
}

void DomainInp::CalcDomainPatch( int iZone, int jZone, GridMediator * gridMediator )
{
    Grids grids = gridMediator->gridVector;

    StrGrid * grid_i = ONEFLOW::StrGridCast( grids[ iZone ] );
    StrGrid * grid_j = ONEFLOW::StrGridCast( grids[ jZone ] );

    IjkBox ijkBox_i;
    ijkBox_i.CreateBox( grid_i );

    IjkBox ijkBox_j;
    ijkBox_j.CreateBox( grid_j );

    for ( int i = 0; i < ijkBox_i.idir.size(); ++ i )
    {
        for ( int j = 0; j < ijkBox_j.idir.size(); ++ j )
        {
        }
    }
}

void DomainInp::CalcFacePoint( StrGrid * grid, PointSearch * pointSearch, IjkBox * ijkBox, int zId, PBlkSet * pblkSet )
{
    Field3D & xs = * grid->strx;
    Field3D & ys = * grid->stry;
    Field3D & zs = * grid->strz;

    RealField coor( 3 );

    int nDomain = ijkBox->imin.size();

    for ( int n = 0; n < nDomain; ++ n )
    {
        int imin = ijkBox->imin[ n ];
        int imax = ijkBox->imax[ n ];
        int jmin = ijkBox->jmin[ n ];
        int jmax = ijkBox->jmax[ n ];
        int kmin = ijkBox->kmin[ n ];
        int kmax = ijkBox->kmax[ n ];

        for ( int k = kmin; k <= kmax; ++ k )
        {
            for ( int j = jmin; j <= jmax; ++ j )
            {
                for ( int i = imin; i <= imax; ++ i )
                {
                    PBlk * pblk = new PBlk();
                    pblk->blk = zId;
                    pblk->fid = n;
                    pblk->i = i;
                    pblk->j = j;
                    pblk->k = k;

                    Real xm = xs( i, j, k );
                    Real ym = ys( i, j, k );
                    Real zm = zs( i, j, k );

                    int pid = pointSearch->AddPoint( xm, ym, zm );
                    pblkSet->Add( pid, pblk );
                    int kkk = 1;
                }
            }
        }
    }
}

EndNameSpace