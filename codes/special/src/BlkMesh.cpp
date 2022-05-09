/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2022 He Xin and the OneFLOW contributors.
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

#include "Prj.h"
#include "Dimension.h"
#include "BlockFaceSolver.h"
#include "Transfinite.h"
#include "StrGrid.h"
#include "NodeMesh.h"
#include <fstream>
#include <iomanip>


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

void BlkBasic::DumpInp( std::fstream & file )
{
    int width = 5;

    file << std::setw( width ) << ni;
    file << std::setw( width ) << nj;
    if ( Dim::dimension == ONEFLOW::THREE_D )
    {
        file << std::setw( width ) << nk;
    }
    file << std::endl;

    file << "zone" << this->blk_id + 1 << "\n";

    int nFaces = this->facelist.size();

    file << std::setw( width )  << nFaces << "\n";

    for ( int iFace = 0; iFace < nFaces; ++ iFace )
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

        file << std::setw( width ) << imin;
        file << std::setw( width ) << imax;
        file << std::setw( width ) << jmin;
        file << std::setw( width ) << jmax;
        if ( Dim::dimension == ONEFLOW::THREE_D )
        {
            file << std::setw( width ) << kmin;
            file << std::setw( width ) << kmax;
        }
        file << std::setw( width ) << bcType;
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

            file << std::setw( width ) << imin;
            file << std::setw( width ) << imax;
            file << std::setw( width ) << jmin;
            file << std::setw( width ) << jmax;
            if ( Dim::dimension == ONEFLOW::THREE_D )
            {
                file << std::setw( width ) << kmin;
                file << std::setw( width ) << kmax;
            }
            file << std::setw( width ) << zid;
            file << "\n";
        }
    }
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

void SetGridXYZ2D( StrGrid * grid, RealField2D & x2d, RealField2D & y2d, RealField2D & z2d )
{
    int ni = grid->ni;
    int nj = grid->nj;
    int nk = grid->nk;
    int pos = 0;

    for ( int j = 0; j < nj; ++ j )
    {
        for ( int i = 0; i < ni; ++ i )
        {
            grid->nodeMesh->xN[ pos ] = x2d[ i ][ j ];
            grid->nodeMesh->yN[ pos ] = y2d[ i ][ j ];
            grid->nodeMesh->zN[ pos ] = z2d[ i ][ j ];
            ++ pos;
        }
    }
}

EndNameSpace
