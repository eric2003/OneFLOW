/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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
#include "Block2D.h"
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

Block2D::Block2D()
{
    int nMDomain = 1;
    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MDomain * mDomain = new MDomain();
        mDomain->pos = iMDomain;
        mDomain->coorMap = & this->coorMap;
        mDomainList.push_back( mDomain );
    }

    this->AddLocalPt( 1, 2, 3, 4 );
}

Block2D::~Block2D()
{
    DeletePointer( mLineList );
    DeletePointer( facelist );
}

void Block2D::Alloc()
{
    AllocateVector( x2d, ni, nj );
    AllocateVector( y2d, ni, nj );
    AllocateVector( z2d, ni, nj );
}

void Block2D::CreateBlockMesh2D()
{
    int nMDomain = mDomainList.size();
    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MDomain * mDomain = mDomainList[ iMDomain ];
        mDomain->SetBlkBcMesh( this );
    }

    this->GenerateBlockMesh2D();
}


void Block2D::GenerateBlockMesh2D()
{
}

void Block2D::DumpBlockMesh2D( std::fstream &file )
{
    file << " VARIABLES = \"X\", \"Y\", \"Z\" \n";
    file << " ZONE I = " << ni << ", J = " << nj << " F = POINT \n";
    for ( int j = 0; j < nj; ++ j )
    {
        for ( int i = 0; i < ni; ++ i )
        {
            file << x2d[ i ][ j ] << " " << y2d[ i ][ j ] << " " << z2d[ i ][ j ] << "\n";
        }
    }
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
    int nMDomain = mDomainList.size();
    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MDomain * mDomain = mDomainList[ iMDomain ];
        mDomain->ConstructMultiDomainTopo();
    }

    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MDomain * mDomain = mDomainList[ iMDomain ];
        mDomain->CalcDomainCtrlPoints();
    }

    MDomain * mDomain = mDomainList[ 0 ];

    this->controlpoints = mDomain->ctrlpoints;

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
    int nFaces = this->facelist.size();
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        Face2D * face2d = this->facelist[ iFace ];
        int domain_id = face2d->face_id;
        if ( face2d->bcType == -1 )
        {
            face2d->t = new Face2D();

            BlkF2C & face_struct = blkFaceSolver.line2Face[ domain_id - 1 ];
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
    int nMDomain = mDomainList.size();
    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MDomain * mDomain = mDomainList[ iMDomain ];
        mDomain->CalcDim2D();
    }

    MDomain * d = mDomainList[ 0 ];

    this->ni = d->ni;
    this->nj = d->nj;
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
        this->coorMap.insert( std::pair<int, CalcCoor>( pt, c ) );
    }

    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MDomain * mDomain = mDomainList[ iMDomain ];
        mDomain->CalcCoor();
    }

    CreateFaceList();

    int kkk = 1;
}

void Block2D::CreateFaceList()
{
    int nMDomain = mDomainList.size();
    for ( int iMDomain = 0; iMDomain < nMDomain; ++ iMDomain )
    {
        MDomain * mDomain = mDomainList[ iMDomain ];
        mDomain->CreateInpFaceList1D( facelist );
    }

    int kkk = 1;
}

void Block2D::FillStrGrid( Grid * gridIn, int iZone )
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
    grid->nodeMesh->CreateNodes( grid->nNodes );

    grid->SetLayout();

    SetGridXYZ2D( grid, x2d, y2d, z2d );
}


EndNameSpace
