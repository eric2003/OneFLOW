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

#include "BlockElem.h"
#include "HXCgns.h"
#include "HXMath.h"

BeginNameSpace( ONEFLOW )

BlkElemHome bbElemHome;

BlkElem::BlkElem()
{
}

BlkElem::~BlkElem()
{
}

void BlkElem::Init( int eType )
{
    this->eType = eType;
    int nNode = 4;
    if ( eType == ONEFLOW::QUAD_4 )
    {
        nNode = 4;
    }
    else if ( eType == ONEFLOW::HEXA_8 )
    {
        nNode = 8;
    }

    nodeId.resize( nNode );

    for ( int iNode = 0; iNode < nNode; ++ iNode )
    {
        nodeId[ iNode ] = iNode;
    }

    if ( eType == ONEFLOW::QUAD_4 )
    {
        // QUAD_4
        // 1....2....3 ...4....1

        // Face 4
        // 1 2
        // 2 3
        // 3 4
        // 4 1
        this->PushElementFace( ONEFLOW::BAR_2, 1, 2 );
        this->PushElementFace( ONEFLOW::BAR_2, 2, 3 );
        this->PushElementFace( ONEFLOW::BAR_2, 3, 4 );
        this->PushElementFace( ONEFLOW::BAR_2, 4, 1 );
    }
    else if ( eType == ONEFLOW::HEXA_8 )
    {
        // HEXA_8

        // 1 2 3 4 5 6 7 8
        // Face 6
        // 1, 4, 8, 5 //0
        // 2, 3, 7, 6 //1
        // 1, 2, 6, 5 //2
        // 4, 3, 7, 8 //3
        // 1, 2, 3, 4 //4
        // 5, 6, 7, 8 //5
        this->PushElementFace( ONEFLOW::QUAD_4, 1, 4, 8, 5 );
        this->PushElementFace( ONEFLOW::QUAD_4, 2, 3, 7, 6 );
        this->PushElementFace( ONEFLOW::QUAD_4, 1, 2, 6, 5 );
        this->PushElementFace( ONEFLOW::QUAD_4, 4, 3, 7, 8 );
        this->PushElementFace( ONEFLOW::QUAD_4, 1, 2, 3, 4 );
        this->PushElementFace( ONEFLOW::QUAD_4, 5, 6, 7, 8 );
    }
}

void BlkElem::PushElementFace( int faceType, int p1, int p2 )
{
    faceTypeList.push_back( faceType );

    IntField face;
    face.push_back( p1 );
    face.push_back( p2 );
    for ( int i = 0; i < face.size(); ++ i )
    {
        face[ i ] -= 1;
    }
    faceList.push_back( face );
}

void BlkElem::PushElementFace( int faceType, int p1, int p2, int p3, int p4 )
{
    faceTypeList.push_back( faceType );

    IntField face;
    face.push_back( p1 );
    face.push_back( p2 );
    face.push_back( p3 );
    face.push_back( p4 );
    for ( int i = 0; i < face.size(); ++ i )
    {
        face[ i ] -= 1;
    }

    faceList.push_back( face );
}

BlkElemHome::BlkElemHome()
{
    initFlag = false;
}

BlkElemHome::~BlkElemHome()
{
    Free();
}

void BlkElemHome::Init()
{
    if ( initFlag ) return;
    initFlag = true;

    int n1 = ONEFLOW::QUAD_4;
    int n2 = ONEFLOW::HEXA_8;
    int nSize = MAX( n1, n2 ) + 1;

    elems.resize( nSize );
    IntField a;
    a.push_back( n1 );
    a.push_back( n2 );

    for ( int i = 0; i < a.size(); ++ i )
    {
        int eType = a[ i ];
        BlkElem * elem = new BlkElem();
        elem->Init( eType );
        elems[ eType ] = elem;
    }
}

void BlkElemHome::Free()
{
    int nSize = elems.size();
    for ( int i = 0; i < nSize; ++ i )
    {
        delete elems[ i ];
    }
}

BlkElem * BlkElemHome::GetBlkElem( int eType )
{
    return elems[ eType ];
}

BlkFace::BlkFace()
{
    ;
}

BlkFace::~BlkFace()
{
    ;
}

EndNameSpace