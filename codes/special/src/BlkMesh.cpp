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

#include "BlkMesh.h"
#include "MLine.h"
#include "MDomain.h"
#include "FileUtil.h"
#include "Prj.h"
#include "Dimension.h"
#include "BlockFaceSolver.h"
#include "Transfinite.h"
#include "StrGrid.h"
#include "NodeMesh.h"
#include <fstream>
#include <iomanip>
using namespace std;

BeginNameSpace( ONEFLOW )

BlkBasic::BlkBasic()
{
    ;
}

BlkBasic::~BlkBasic()
{
    ;
}

void BlkBasic::AddLocalPt( int p1, int p2 )
{
    IntField face;
    face.push_back( p1 );
    face.push_back( p2 );
    this->localpt.push_back( face );
}

void BlkBasic::AddLocalPt( int p1, int p2, int p3, int p4 )
{
    IntField face;
    face.push_back( p1 );
    face.push_back( p2 );
    face.push_back( p3 );
    face.push_back( p4 );
    this->localpt.push_back( face );
}

void BlkBasic::Add( IntField &iList, IntField &jList, IntField &kList, int i, int j, int k )
{
    iList.push_back( i );
    jList.push_back( j );
    kList.push_back( k );
}

void BlkBasic::DumpInp( fstream & file )
{
    int width = 5;

    file << setw( width ) << ni;
    file << setw( width ) << nj;
    if ( Dim::dimension == ONEFLOW::THREE_D )
    {
        file << setw( width ) << nk;
    }
    file << endl;

    file << "zone" << this->blk_id + 1 << "\n";

    int nFace = this->GetNSubDomain();

    file << setw( width )  << nFace << "\n";

    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        Face2D * face2d = this->facelist[ iFace ];
        face2d->CalcRegion();

        int imin = face2d->st.i;
        int jmin = face2d->st.j;
        int kmin = face2d->st.k;

        int imax = face2d->ed.i;
        int jmax = face2d->ed.j;
        int kmax = face2d->ed.k;

        int bcType = face2d->bcType;

        file << setw( width ) << imin;
        file << setw( width ) << imax;
        file << setw( width ) << jmin;
        file << setw( width ) << jmax;
        if ( Dim::dimension == ONEFLOW::THREE_D )
        {
            file << setw( width ) << kmin;
            file << setw( width ) << kmax;
        }
        file << setw( width ) << bcType;
        file << "\n";

        if ( bcType < 0 )
        {
            imin = face2d->t->st.i;
            jmin = face2d->t->st.j;
            kmin = face2d->t->st.k;

            imax = face2d->t->ed.i;
            jmax = face2d->t->ed.j;
            kmax = face2d->t->ed.k;

            int zid = face2d->t->bcType;

            file << setw( width ) << imin;
            file << setw( width ) << imax;
            file << setw( width ) << jmin;
            file << setw( width ) << jmax;
            if ( Dim::dimension == ONEFLOW::THREE_D )
            {
                file << setw( width ) << kmin;
                file << setw( width ) << kmax;
            }
            file << setw( width ) << zid;
            file << "\n";
        }
    }
}

Block3D::Block3D()
{
    int nMDomain = 6;
    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MDomain * mDomain = new MDomain();
        mDomain->pos = iMDomain;
        mDomain->coorMap = & this->coorMap;
        mDomainList.push_back( mDomain );
    }

    this->AddLocalPt( 1, 4, 8, 5 );
    this->AddLocalPt( 2, 3, 7, 6 );
    this->AddLocalPt( 1, 2, 6, 5 );
    this->AddLocalPt( 4, 3, 7, 8 );
    this->AddLocalPt( 1, 2, 3, 4 );
    this->AddLocalPt( 5, 6, 7, 8 );
}

Block3D::~Block3D()
{
    DeletePointer( mDomainList );
    DeletePointer( facelist );
}

void Block3D::Alloc()
{
    AllocateVector( x3d, ni, nj, nk );
    AllocateVector( y3d, ni, nj, nk );
    AllocateVector( z3d, ni, nj, nk );
}

int Block3D::GetNSubDomain()
{
    int nSubDomain = 0;
    int nMDomain = mDomainList.size();
    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MDomain * mDomain = mDomainList[ iMDomain ];
        nSubDomain += mDomain->GetNsubDomain();
    }
    return nSubDomain;
}

void Block3D::ConstructTopo()
{
    int nMDomain = mDomainList.size();
    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MDomain * mDomain = mDomainList[ iMDomain ];
        mDomain->ConstructMultiDomainTopo();
    }
    map< int, IntSet > p2dMap;
    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MDomain * mDomain = mDomainList[ iMDomain ];
        int nSize = mDomain->candidate_bcpoints.size();
        for ( int i = 0; i < nSize; ++ i )
        {
            int pt = mDomain->candidate_bcpoints[ i ];
            map< int, IntSet >::iterator iter;
            iter = p2dMap.find( pt );
            if ( iter == p2dMap.end() )
            {
                IntSet iset;
                iset.insert( iMDomain );
                p2dMap.insert( pair< int, IntSet >( pt, iset ) );
            }
            else
            {
                iter->second.insert( iMDomain );
            }
        }
        int kkk = 1;
    }

    IntField ctrl_points;
    map< int, IntSet >::iterator iter;
    for ( iter = p2dMap.begin(); iter != p2dMap.end(); ++ iter )
    {
        if ( iter->second.size() == 3 )
        {
            ctrl_points.push_back( iter->first );
        }
    }

    //for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    //{
    //    MDomain * mDomain = mDomainList[ iMDomain ];
    //    mDomain->CalcDomainCtrlPoints( ctrl_points );
    //}

    int p1, p2, p3, p4, p5, p6, p7, p8;

    GetCornerPoint( p1, 0, 2, 4 );
    GetCornerPoint( p2, 1, 2, 4 );
    GetCornerPoint( p3, 1, 3, 4 );
    GetCornerPoint( p4, 0, 3, 4 );
    GetCornerPoint( p5, 0, 2, 5 );
    GetCornerPoint( p6, 1, 2, 5 );
    GetCornerPoint( p7, 1, 3, 5 );
    GetCornerPoint( p8, 0, 3, 5 );

    this->controlpoints.push_back( p1 );
    this->controlpoints.push_back( p2 );
    this->controlpoints.push_back( p3 );
    this->controlpoints.push_back( p4 );
    this->controlpoints.push_back( p5 );
    this->controlpoints.push_back( p6 );
    this->controlpoints.push_back( p7 );
    this->controlpoints.push_back( p8 );

    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MDomain * mDomain = mDomainList[ iMDomain ];
        mDomain->CalcDomainCtrlPoints( this->controlpoints, this->localpt[iMDomain ] );
    }
    this->CalcBlkDim();
}

void Block3D::GetCornerPoint( int & pt, int id1, int id2, int id3 )
{
    MDomain * d1 = mDomainList[ id1 ];
    MDomain * d2 = mDomainList[ id2 ];
    MDomain * d3 = mDomainList[ id3 ];

    int nSize = d1->candidate_ctrlpoints.size();
    for ( int i = 0; i < nSize; ++ i )
    {
        int ip = d1->candidate_ctrlpoints[ i ];
        bool flag1 = InArray( ip, d2->candidate_ctrlpoints );
        if ( ! flag1 ) continue;
        bool flag2 = InArray( ip, d3->candidate_ctrlpoints );
        if ( ! flag2 ) continue;
        pt = ip;
        break;
    }
}

void Block3D::SetInterfaceBc()
{
    int nFace = this->facelist.size();
    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        Face2D * face2d = this->facelist[ iFace ];
        int domain_id = face2d->face_id;
        if ( face2d->bcType == -1 )
        {
            face2d->t = new Face2D();

            BlkF2C & face_struct = blkFaceSolver.myFaceSolver.face2Block[ domain_id ];
            int n_neibor = face_struct.cellList.size();
            int blk1 = face_struct.cellList[ 0 ] - 1;
            int blk2 = face_struct.cellList[ 1 ] - 1;

            int tblk = -1;
            if ( this->blk_id == blk1 )
            {
                tblk = blk2;
            }
            else
            {
                tblk = blk1;
            }
            Face2D * facet = blkFaceSolver.GetBlkFace2D( tblk, domain_id );
            face2d->t->bcType = tblk + 1;
            face2d->t->st = facet->st;
            face2d->t->ed = facet->ed;
        }
    }
}

void Block3D::CalcBlkDim()
{
    int nMDomain = mDomainList.size();
    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MDomain * mDomain = mDomainList[ iMDomain ];
        mDomain->CalcDim2D();
    }

    MDomain * d0 = mDomainList[ 0 ];
    MDomain * d2 = mDomainList[ 2 ];
    //d0: nj,nk
    //d2: ni,nk

    this->ni = d2->ni;
    this->nj = d0->ni;
    this->nk = d0->nj;

    IntField iList, jList, kList;
    Add( iList, jList, kList, 1, 1, 1 );
    Add( iList, jList, kList, ni, 1, 1 );
    Add( iList, jList, kList, ni, nj, 1 );
    Add( iList, jList, kList, 1, nj, 1 );
    Add( iList, jList, kList, 1, 1, nk );
    Add( iList, jList, kList, ni, 1, nk );
    Add( iList, jList, kList, ni, nj, nk );
    Add( iList, jList, kList, 1, nj, nk );

    for ( int iPoint = 0; iPoint < iList.size(); ++ iPoint )
    {
        int pt = this->controlpoints[ iPoint ];
        int i = iList[ iPoint ];
        int j = jList[ iPoint ];
        int k = kList[ iPoint ];
        CalcCoor c;
        c.SetCoor( i, j, k );
        this->coorMap.insert( pair<int, CalcCoor>( pt, c ) );
    }

    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MDomain * mDomain = mDomainList[ iMDomain ];
        mDomain->CalcCoor();
    }

    CreateFaceList();

    int kkk = 1;
}

void Block3D::CreateFaceList()
{
    int nMDomain = mDomainList.size();
    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MDomain * mDomain = mDomainList[ iMDomain ];
        mDomain->CreateInpFaceList( facelist );
    }

    int kkk = 1;

}


void Block3D::CreateBlockMesh()
{
    int nMDomain = mDomainList.size();
    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MDomain * mDomain = mDomainList[ iMDomain ];
        mDomain->SetBlkBcMesh( this );
    }

    this->GenerateBlockMesh();
}

void Block3D::GenerateBlockMesh()
{
    fstream file;
    OpenPrjFile( file, "grid/blkfaceplot.dat", ios_base::out );
    file << " VARIABLES = \"X\", \"Y\", \"Z\" \n";
    file << " ZONE I = " << ni << ", J = " << nj << " F = POINT \n";
    for ( int j = 0; j < nj; ++ j )
    {
        for ( int i = 0; i < ni; ++ i )
        {
            file << x3d[ i ][ j ][ 0 ] << " " << y3d[ i ][ j ][ 0 ] << " " << z3d[ i ][ j ][ 0 ] << "\n";
        }
    }
    file << " ZONE I = " << ni << ", J = " << nj << " F = POINT \n";
    for ( int j = 0; j < nj; ++ j )
    {
        for ( int i = 0; i < ni; ++ i )
        {
            file << x3d[ i ][ j ][ nk-1 ] << " " << y3d[ i ][ j ][ nk-1 ] << " " << z3d[ i ][ j ][ nk-1 ] << "\n";
        }
    }
    file << " ZONE I = " << ni << ", J = " << nk << " F = POINT \n";
    for ( int k = 0; k < nk; ++ k )
    {
        for ( int i = 0; i < ni; ++ i )
        {
            file << x3d[ i ][ 0 ][ k ] << " " << y3d[ i ][ 0 ][ k ] << " " << z3d[ i ][ 0 ][ k ] << "\n";
        }
    }
    file << " ZONE I = " << ni << ", J = " << nk << " F = POINT \n";
    for ( int k = 0; k < nk; ++ k )
    {
        for ( int i = 0; i < ni; ++ i )
        {
            file << x3d[ i ][ nj - 1 ][ k ] << " " << y3d[ i ][ nj - 1 ][ k ] << " " << z3d[ i ][ nj - 1 ][ k ] << "\n";
        }
    }
    file << " ZONE I = " << nj << ", J = " << nk << " F = POINT \n";
    for ( int k = 0; k < nk; ++ k )
    {
        for ( int j = 0; j < nj; ++ j )
        {
            file << x3d[ 0 ][ j ][ k ] << " " << y3d[ 0 ][ j ][ k ] << " " << z3d[ 0 ][ j ][ k ] << "\n";
        }
    }
    file << " ZONE I = " << nj << ", J = " << nk << " F = POINT \n";
    for ( int k = 0; k < nk; ++ k )
    {
        for ( int j = 0; j < nj; ++ j )
        {
            file << x3d[ ni - 1 ][ j ][ k ] << " " << y3d[ ni-1 ][ j ][ k ] << " " << z3d[ ni-1 ][ j ][ k ] << "\n";
        }
    }
    CloseFile( file );
    TransfiniteInterpolation( x3d, ni, nj, nk );
    TransfiniteInterpolation( y3d, ni, nj, nk );
    TransfiniteInterpolation( z3d, ni, nj, nk );
    //for ( int k = 1; k < nk - 1; ++ k )
    //{
    //    for ( int j = 1; j < nj - 1; ++ j )
    //    {
    //        for ( int i = 1; i < ni - 1; ++ i )
    //        {
    //            Real d = ( x3d[ i ][ j ][ nk - 1 ] - x3d[ i ][ j ][ 0 ] ) / ( nk - 1 );
    //            x3d[ i ][ j ][ k ] = x3d[ i ][ j ][ 0 ] + k * d;
    //            y3d[ i ][ j ][ k ] = y3d[ i ][ j ][ 0 ] + k * d;
    //            z3d[ i ][ j ][ k ] = z3d[ i ][ j ][ 0 ] + k * d;
    //        }
    //    }
    //}


    OpenPrjFile( file, "grid/blkplot.dat", ios_base::out );
    file << " VARIABLES = \"X\", \"Y\", \"Z\" \n";
    file << " ZONE I = " << ni << ", J = " << nj << ", K = " << nk << " F = POINT \n";
    for ( int k = 0; k < nk; ++ k )
    {
        for ( int j = 0; j < nj; ++ j )
        {
            for ( int i = 0; i < ni; ++ i )
            {
                file << x3d[ i ][ j ][ k ] << " " << y3d[ i ][ j ][ k ] << " " << z3d[ i ][ j ][ k ] << "\n";
            }
        }
    }
    CloseFile( file );

    OpenPrjFile( file, "grid/blkfaceplot111.dat", ios_base::out );
    file << " VARIABLES = \"X\", \"Y\", \"Z\" \n";
    file << " ZONE I = " << ni << ", J = " << nj << " F = POINT \n";
    for ( int j = 0; j < nj; ++ j )
    {
        for ( int i = 0; i < ni; ++ i )
        {
            file << x3d[ i ][ j ][ 0 ] << " " << y3d[ i ][ j ][ 0 ] << " " << z3d[ i ][ j ][ 0 ] << "\n";
        }
    }
    file << " ZONE I = " << ni << ", J = " << nj << " F = POINT \n";
    for ( int j = 0; j < nj; ++ j )
    {
        for ( int i = 0; i < ni; ++ i )
        {
            file << x3d[ i ][ j ][ nk-1 ] << " " << y3d[ i ][ j ][ nk-1 ] << " " << z3d[ i ][ j ][ nk-1 ] << "\n";
        }
    }
    file << " ZONE I = " << ni << ", J = " << nk << " F = POINT \n";
    for ( int k = 0; k < nk; ++ k )
    {
        for ( int i = 0; i < ni; ++ i )
        {
            file << x3d[ i ][ 0 ][ k ] << " " << y3d[ i ][ 0 ][ k ] << " " << z3d[ i ][ 0 ][ k ] << "\n";
        }
    }
    file << " ZONE I = " << ni << ", J = " << nk << " F = POINT \n";
    for ( int k = 0; k < nk; ++ k )
    {
        for ( int i = 0; i < ni; ++ i )
        {
            file << x3d[ i ][ nj - 1 ][ k ] << " " << y3d[ i ][ nj - 1 ][ k ] << " " << z3d[ i ][ nj - 1 ][ k ] << "\n";
        }
    }
    file << " ZONE I = " << nj << ", J = " << nk << " F = POINT \n";
    for ( int k = 0; k < nk; ++ k )
    {
        for ( int j = 0; j < nj; ++ j )
        {
            file << x3d[ 0 ][ j ][ k ] << " " << y3d[ 0 ][ j ][ k ] << " " << z3d[ 0 ][ j ][ k ] << "\n";
        }
    }
    file << " ZONE I = " << nj << ", J = " << nk << " F = POINT \n";
    for ( int k = 0; k < nk; ++ k )
    {
        for ( int j = 0; j < nj; ++ j )
        {
            file << x3d[ ni - 1 ][ j ][ k ] << " " << y3d[ ni-1 ][ j ][ k ] << " " << z3d[ ni-1 ][ j ][ k ] << "\n";
        }
    }
    CloseFile( file );
    int kkk = 1;
}

void Block3D::FillStrGrid( Grid * gridIn, int iZone )
{
    StrGrid * grid = StrGridCast( gridIn );
    int ni = this->ni;
    int nj = this->nj;
    int nk = this->nk;

    grid->id = iZone;
    grid->ni = ni;
    grid->nj = nj;
    grid->nk = nk;
    grid->SetBasicDimension();
    grid->nodeMesh->CreateNodes( grid->nNode );

    grid->SetLayout();

    SetGridXYZ3D( grid, x3d, y3d, z3d );
}

void SetGridXYZ3D( StrGrid * grid, RealField3D & x3d, RealField3D & y3d, RealField3D & z3d )
{
    int ni = grid->ni;
    int nj = grid->nj;
    int nk = grid->nk;
    int pos = 0;

    for ( int k = 0; k < nk; ++ k )
    {
        for ( int j = 0; j < nj; ++ j )
        {
            for ( int i = 0; i < ni; ++ i )
            {
                grid->nodeMesh->xN[ pos ] = x3d[ i ][ j ][ k ];
                grid->nodeMesh->yN[ pos ] = y3d[ i ][ j ][ k ];
                grid->nodeMesh->zN[ pos ] = z3d[ i ][ j ][ k ];
                ++ pos;
            }
        }
    }
}

Block2D::Block2D()
{
    int nMLine = 4;
    for ( int iMLine = 0; iMLine < nMLine; ++ iMLine )
    {
        MLine * mLine = new MLine();
        mLine->pos = iMLine;
        mLine->coorMap = & this->coorMap;
        mLineList.push_back( mLine );
    }

    this->AddLocalPt( 1, 2 );
    this->AddLocalPt( 2, 3 );
    this->AddLocalPt( 3, 4 );
    this->AddLocalPt( 4, 1 );
}

Block2D::~Block2D()
{
    DeletePointer( mLineList );
    DeletePointer( facelist );
}

int Block2D::GetNSubDomain()
{
    int nSubDomain = 0;
    int nMDomain = mLineList.size();
    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MLine * mLine = mLineList[ iMDomain ];
        nSubDomain += mLine->slineList.size();
    }
    return nSubDomain;
}

void Block2D::ConstructTopo()
{
    int nMDomain = mLineList.size();
    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MLine * mLine = mLineList[ iMDomain ];
        mLine->ConstructDomainTopo();
    }

    int p1, p2, p3, p4;

    GetCornerPoint( p1, 3, 0 );
    GetCornerPoint( p2, 0, 1 );
    GetCornerPoint( p3, 1, 2 );
    GetCornerPoint( p4, 2, 3 );

    this->controlpoints.push_back( p1 );
    this->controlpoints.push_back( p2 );
    this->controlpoints.push_back( p3 );
    this->controlpoints.push_back( p4 );

    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MLine * mLine = mLineList[ iMDomain ];
        mLine->CalcDomainCtrlPoints( this->controlpoints, this->localpt[iMDomain ] );
    }

    this->CalcBlkDim();
}

void Block2D::GetCornerPoint( int & pt, int id1, int id2 )
{
    MLine * d1 = mLineList[ id1 ];
    MLine * d2 = mLineList[ id2 ];

    int nSize = d1->candidate_ctrlpoints.size();
    for ( int i = 0; i < nSize; ++ i )
    {
        int ip = d1->candidate_ctrlpoints[ i ];
        bool flag1 = InArray( ip, d2->candidate_ctrlpoints );
        if ( ! flag1 ) continue;
        pt = ip;
        break;
    }
}

void Block2D::SetInterfaceBc()
{
    int nFace = this->facelist.size();
    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        Face2D * face2d = this->facelist[ iFace ];
        int domain_id = face2d->face_id;
        if ( face2d->bcType == -1 )
        {
            face2d->t = new Face2D();

            BlkF2C & face_struct = blkFaceSolver.myFaceSolver.face2Block[ domain_id ];
            int n_neibor = face_struct.cellList.size();
            int blk1 = face_struct.cellList[ 0 ] - 1;
            int blk2 = face_struct.cellList[ 1 ] - 1;

            int tblk = -1;
            if ( this->blk_id == blk1 )
            {
                tblk = blk2;
            }
            else
            {
                tblk = blk1;
            }
            Face2D * facet = blkFaceSolver.GetBlkFace2D( tblk, domain_id );
            face2d->t->bcType = tblk + 1;
            face2d->t->st = facet->st;
            face2d->t->ed = facet->ed;
        }
    }
}

void Block2D::CalcBlkDim()
{
    int nMDomain = mLineList.size();
    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MLine * mLine = mLineList[ iMDomain ];
        mLine->CalcDim1D();
    }

    MLine * d0 = mLineList[ 0 ];
    MLine * d1 = mLineList[ 1 ];
    this->ni = d0->ni;
    this->nj = d1->ni;
    this->nk = 1;

    IntField iList, jList, kList;
    Add( iList, jList, kList, 1, 1, 1 );
    Add( iList, jList, kList, ni, 1, 1 );
    Add( iList, jList, kList, ni, nj, 1 );
    Add( iList, jList, kList, 1, nj, 1 );

    for ( int iPoint = 0; iPoint < iList.size(); ++ iPoint )
    {
        int pt = this->controlpoints[ iPoint ];
        int i = iList[ iPoint ];
        int j = jList[ iPoint ];
        int k = kList[ iPoint ];
        CalcCoor c;
        c.SetCoor( i, j, k );
        this->coorMap.insert( pair<int, CalcCoor>( pt, c ) );
    }

    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MLine * mLine = mLineList[ iMDomain ];
        mLine->CalcCoor( & this->coorMap );
    }

    CreateFaceList();

    int kkk = 1;
}

void Block2D::CreateFaceList()
{
    int nMDomain = mLineList.size();
    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MLine * mLine = mLineList[ iMDomain ];
        mLine->CreateInpFaceList( facelist );
    }

    int kkk = 1;

}

EndNameSpace