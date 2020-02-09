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

#include "Plate.h"
#include "NodeField.h"
#include "Zone.h"
#include "UnsGrid.h"
#include "NodeMesh.h"
#include "FaceTopo.h"
#include "BcRecord.h"
#include "PointSearch.h"
#include "HXMath.h"
#include "Parallel.h"
#include "DataBook.h"
#include "DataBase.h"
#include "NsCom.h"
#include "Boundary.h"
#include <algorithm>
#include <iomanip>
using namespace std;

BeginNameSpace( ONEFLOW )

SliceInfo::SliceInfo()
{
    ;
}

SliceInfo::~SliceInfo()
{
    ;
}

void SliceInfo::AddSlice( Real pos, int dir1, int dir2 )
{
    this->slicepos.push_back( pos );
    this->dir1.push_back( dir1 );
    this->dir2.push_back( dir2 );
}

PlaneData::PlaneData()
{
    ;
}

PlaneData::~PlaneData()
{
    ;
}

void PlaneData::Init()
{
    int nEqu = nodedata->GetNEqu();
    var1.resize( nEqu );
    var2.resize( nEqu );
    var .resize( nEqu );
    slicedata.resize( nEqu );
}

void PlaneData::AddPoint( Real xm, Real ym, Real zm )
{
    this->x.push_back( xm );
    this->y.push_back( ym );
    this->z.push_back( zm );
}

void PlaneData::AddVar( RealField & var )
{
    for ( int iEqu = 0; iEqu < var.size(); ++ iEqu )
    {
        this->slicedata[ iEqu ].push_back( var[ iEqu ] );
    }
}

void PlaneData::AddVar( int p1, int p2, Real c1, Real c2 )
{
    int nEqu = var.size();

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        var1[ iEqu ] = ( * nodedata )[ iEqu ][ p1 ];
        var2[ iEqu ] = ( * nodedata )[ iEqu ][ p2 ];
    }

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        var[ iEqu ] = c1 * var1[ iEqu ] + c2 * var2[ iEqu ];
    }

    this->AddVar( var );
}

void PlaneData::Write( DataBook * dataBook )
{
    int nNode = x.size();
    HXWrite( dataBook, nNode );
    HXWrite( dataBook, x );
    HXWrite( dataBook, y );
    HXWrite( dataBook, z );
    int nEqu = var.size();
    HXWrite( dataBook, nEqu );
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        int nElem = this->slicedata[ iEqu ].size();
        HXWrite( dataBook, nElem );
        HXWrite( dataBook, this->slicedata[ iEqu ] );
    }
}

void PlaneData::Read( DataBook * dataBook )
{
    int nNode = 0;
    RealField xx( nNode ), yy( nNode ), zz( nNode );
    HXRead( dataBook, nNode );
    HXRead( dataBook, xx );
    HXRead( dataBook, yy );
    HXRead( dataBook, zz );
    for ( int i = 0; i < nNode; ++ i )
    {
        x.push_back( xx[ i ] );
        y.push_back( yy[ i ] );
        z.push_back( zz[ i ] );
    }

    int nEqu = 0;
    HXRead( dataBook, nEqu );
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        int nElem = 0;
        HXRead( dataBook, nElem );
        RealField sdata( nElem );
        HXRead( dataBook, sdata );
        for ( int i = 0; i < nElem; ++ i )
        {
            this->slicedata[ iEqu ].push_back( sdata[ i ] );
        }
    }
}

void PlaneData::SortData( RealField & varList )
{
    HXVector< HXSort< Real > > sortList;
    HXSort< Real > svar;

    int nNode = this->x.size();
    if ( nNode <= 0 ) return;

    for ( int iNode = 0; iNode < nNode; ++ iNode )
    {
        svar.value = varList[ iNode ];
        svar.index = iNode;
        sortList.push_back( svar );
    }
    sort( sortList.begin(), sortList.end() );

    IntField indexList;
    for ( int iNode = 0; iNode < nNode; ++ iNode )
    {
        indexList.push_back( sortList[ iNode ].index );
    }

    ReorderList( x, indexList );
    ReorderList( y, indexList );
    ReorderList( z, indexList );

    int nEqu = slicedata.size();

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        ReorderList( slicedata[ iEqu ], indexList );
    }
}

void PlaneData::SortDataByAxis( int axis )
{
    if ( axis == X_DIR )
    {
        this->SortData( this->x );
    }
    else if ( axis == Y_DIR )
    {
        this->SortData( this->y );
    }
    else if ( axis == Z_DIR )
    {
        this->SortData( this->z );
    }
}

int PlaneData::FindYIndex()
{
    //找到壁面第一层的点
    RealField yList = this->y;
    sort( yList.begin(), yList.end() );

    int yIndex = 1;
    for ( int iNode = 0; iNode < yList.size(); ++ iNode )
    {
        if ( yList[ iNode ] )
        {
            yIndex = iNode;
            break;
        }
    }

    for ( int iNode = 0; iNode < yList.size(); ++ iNode )
    {
        if ( yList[ yIndex ] == this->y[ iNode ] )
        {
            yIndex = iNode;
            break;
        }
    }

    return yIndex;
}

LamData::LamData()
{
}

LamData::~LamData()
{
    int nField = data.size();
    for ( int i = 0; i < nField; ++ i )
    {
        delete data[ i ];
    }
}

void LamData::Init()
{
    int nField = data.size();
    for ( int i = 0; i < nField; ++ i )
    {
        data[ i ]->Init();
    }
}

void LamData::AddVar( RealField & point, int p1, int p2, Real c1, Real c2 )
{
    int nField = data.size();
    for ( int i = 0; i < nField; ++ i )
    {
        data[ i ]->AddPoint( point[0], point[1], point[2] );
        data[ i ]->AddVar( p1, p2, c1, c2 );
    }
}

int LamData::GetNNode()
{
    return data[ 0 ]->GetNNode();
}

void LamData::Write( DataBook * dataBook )
{
    int nField = data.size();
    HXWrite( dataBook, nField );
    for ( int i = 0; i < nField; ++ i )
    {
        data[ i ]->Write( dataBook );
    }
}

void LamData::Read( DataBook * dataBook )
{
    int nField = 0;
    HXRead( dataBook, nField );

    for ( int i = 0; i < nField; ++ i )
    {
        data[ i ]->Read( dataBook );
    }
}

void LamData::SortDataByAxis( int axis )
{
    int nField = data.size();
    for ( int i = 0; i < nField; ++ i )
    {
        data[ i ]->SortDataByAxis( axis );
    }
}

int LamData::FindYIndex()
{
    return data[ 0 ]->FindYIndex();
}

CuttingClass::CuttingClass()
{
    ;
}

CuttingClass::~CuttingClass()
{
    int nSlice = sliceData.size();
    for ( int i = 0; i < nSlice; ++ i )
    {
        delete sliceData[ i ];
    }
}

void CuttingClass::Init()
{
    int nSlice = sliceInfo.slicepos.size();
    sliceData.resize( nSlice );
    for ( int i = 0; i < nSlice; ++ i )
    {
        sliceData[ i ] = new LamData();
    }
}

void CuttingClass::Slice()
{
    int nField = nameList.size();
    HXVector< MRField * > fields( nField );

    for ( int i = 0; i < nField; ++ i )
    {
        fields[ i ] = CreateNodeVar( nameList[ i ] );
    }

    int nSlice = sliceData.size();
    for ( int iSlice = 0; iSlice < nSlice; ++ iSlice )
    {
        LamData * lam = sliceData[ iSlice ];

        for ( int j = 0; j < nField; ++ j )
        {
            PlaneData * pd = new PlaneData();
            pd->nodedata = fields[ j ];
            lam->data.push_back( pd );
        }
        lam->Init();
    }

    for ( int iSlice = 0; iSlice < nSlice; ++ iSlice )
    {
        this->CutPlane( sliceInfo.slicepos[ iSlice ], sliceInfo.dir1[ iSlice ], this->sliceData[ iSlice ] );
    }

    for ( int i = 0; i < nField; ++ i )
    {
        delete fields[ i ];
    }
}

void CuttingClass::Write( DataBook * dataBook )
{
    int nSlice = sliceData.size();
    HXWrite( dataBook, nSlice );

    for ( int i = 0; i < nSlice; ++ i )
    {
        HXWrite( dataBook, sliceInfo.slicepos[ i ] );
        HXWrite( dataBook, sliceInfo.dir1[ i ] );
        HXWrite( dataBook, sliceInfo.dir2[ i ] );

        LamData * lamData = sliceData[ i ];
        lamData->Write( dataBook );
    }
}

void CuttingClass::Read( DataBook * dataBook )
{
    int nSlice = 0;
    HXRead( dataBook, nSlice );

    for ( int i = 0; i < nSlice; ++ i )
    {
        Real pos;
        int dir1, dir2;
        HXRead( dataBook, pos  );
        HXRead( dataBook, dir1 );
        HXRead( dataBook, dir2 );

        LamData * lamData = sliceData[ i ];
        lamData->Read( dataBook );
    }
}

void CuttingClass::Swap()
{
    int tag = 0;
    if ( Parallel::pid != Parallel::serverid )
    {
        DataBook * dataBook = new DataBook();
        this->Write( dataBook );
        dataBook->Send( Parallel::serverid, tag );
        delete dataBook;
    }
    else
    {
        for ( int pid = 0; pid < Parallel::nProc; ++ pid )
        {
            if ( pid == Parallel::serverid ) continue;
            DataBook * dataBook = new DataBook();
            dataBook->Recv( pid, tag );
            this->Read( dataBook );
            delete dataBook;
        }
    }
}

RealField & CuttingClass::GetCoor( UnsGrid * grid, int cutAxis )
{
    Real * xyz = 0;

    if ( cutAxis == X_DIR )
    {
        return grid->nodeMesh->xN;
    }
    else if ( cutAxis == Y_DIR )
    {
        return grid->nodeMesh->yN;
    }
    else
    {
        return grid->nodeMesh->zN;
    }
}

void CuttingClass::CutPlane( Real cutPosition, int cutAxis, LamData * lamData )
{
    UnsGrid * grid = Zone::GetUnsGrid();
    RealField & x = grid->nodeMesh->xN;
    RealField & y = grid->nodeMesh->yN;
    RealField & z = grid->nodeMesh->zN;

    PointSearch * point_search = new PointSearch();
    point_search->InitializeSpecial( grid, 1.0e-8 );

    Real eps = half * point_search->GetTol();

    RealField & xyz = this->GetCoor( grid, cutAxis );

    int nFace = grid->nFace;
    LinkField & f2n = grid->faceTopo->f2n;

    RealField point( 3 );
    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        Real coorMin =   LARGE;
        Real coorMax = - LARGE;

        int nNode = f2n[ iFace ].size();
        for ( int iNode = 0; iNode < nNode; ++ iNode )
        {
            int nId = f2n[ iFace ][ iNode ];
            coorMin = MIN( xyz[ nId ], coorMin );
            coorMax = MAX( xyz[ nId ], coorMax );
        }

        if ( cutPosition < coorMin - eps || cutPosition > coorMax + eps )
        {
            continue;
        }

        for ( int iNode = 0; iNode < nNode; ++ iNode )
        {
            int p1 = f2n[ iFace ][ iNode     ];
            int p2 = f2n[ iFace ][ ( iNode + 1 ) % nNode ];

            coorMin = MIN( xyz[ p1 ], xyz[ p2 ] );
            coorMax = MAX( xyz[ p1 ], xyz[ p2 ] );

            if ( cutPosition < coorMin - eps ||
                 cutPosition > coorMax + eps )
            {
                continue;
            }

            Real cr = ( cutPosition - xyz[ p1 ] ) / ( xyz[ p2 ] - xyz[ p1 ] + SMALL );
            Real cl = 1 - cr;

            if ( ABS( xyz[ p2 ] - xyz[ p1 ] ) < 1.0e-10 )
            {
                cr = 0.5;
                cl = 0.5;
            }

            Real xl = x[ p1 ];
            Real yl = y[ p1 ];
            Real zl = z[ p1 ];

            Real xr = x[ p2 ];
            Real yr = y[ p2 ];
            Real zr = z[ p2 ];

            Real xm = cl * xl + cr * xr;
            Real ym = cl * yl + cr * yr;
            Real zm = cl * zl + cr * zr;

            point[ 0 ] = xm;
            point[ 1 ] = ym;
            point[ 2 ] = zm;

            if ( point_search->FindPoint( xm, ym, zm ) == INVALID_INDEX )
            {
                int ptId = point_search->AddPoint( xm, ym, zm );
                lamData->AddVar( point, p1, p2, cl, cr );
            }
        }
    }
    delete point_search;
}

EndNameSpace