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

#include "BlockFaceSolver.h"
#include "MLine.h"
#include "SDomain.h"
#include "DataBaseIO.h"
#include "StrGrid.h"
#include "BgGrid.h"
#include "GridState.h"
#include "NodeMesh.h"
#include "MDomain.h"
#include "SimpleDomain.h"
#include "LineMachine.h"
#include "DomainMachine.h"
#include "BlkMesh.h"
#include "HXPointer.h"
#include "CurveInfo.h"
#include "SegmentCtrl.h"
#include "BlockElem.h"
#include "FileUtil.h"
#include "Prj.h"
#include "HXCgns.h"
#include "Dimension.h"
#include <algorithm>
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

BlkFaceSolver blkFaceSolver;

MyFaceSolver::MyFaceSolver()
{
    init_flag = false;
}

MyFaceSolver::~MyFaceSolver()
{
    DeletePointer( sDomainList );
    DeletePointer( slineList );
}

IntField & MyFaceSolver::GetLine( int line_id )
{
    int id = line_id - 1;
    return lineList[ id ];
}

void MyFaceSolver::Alloc()
{
    if ( init_flag ) return;
    init_flag = true;
    int nLine = line_Machine.curveInfoList.size();
    for ( int i = 0; i < nLine; ++ i )
    {
        CurveInfo * curveInfo = line_Machine.GetCurveInfo( i + 1 );
        IntField face;
        face.push_back( curveInfo->p1 );
        face.push_back( curveInfo->p2 );
        this->lineList.push_back( face );
    }
    this->line2Face.resize( nLine );

    for ( int i = 0; i < nLine; ++ i )
    {
        IntField & face = this->lineList[ i ];
        Mid<int> fMid( face.size(), i );
        fMid.id = i;
        fMid.data = face;
        std::sort( fMid.data.begin(), fMid.data.end() );
        this->refLines.insert( fMid );
    }
}

void MyFaceSolver::AddLineToFace( int faceid, int pos, int lineid )
{
    this->Alloc();

    int id = lineid - 1;
    BlkF2C & line_struct = this->line2Face[ id ];

    line_struct.id  = lineid;
    line_struct.type = -1;
    line_struct.bctype = -1;
    line_struct.cellList.push_back( faceid );
    line_struct.posList.push_back( pos );

    faceset.insert( faceid );
}

void MyFaceSolver::CreateFaceList()
{
    int nFace = faceset.size();
    this->faceList.resize( nFace );
    this->facePosList.resize( nFace );

    int nLine = this->line2Face.size();
    for ( int iLine = 0; iLine < nLine; ++ iLine )
    {
        BlkF2C & line_struct = this->line2Face[ iLine ];
        int n = line_struct.cellList.size();
        for ( int i = 0; i < n; ++ i )
        {
            int face_id = line_struct.cellList[ i ] - 1;
            int pos = line_struct.posList[ i ];
            faceList[ face_id ].push_back( line_struct.id );
            this->facePosList[ face_id ].push_back( pos );
        }
    }

    this->face2Block.resize( nFace );
    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        IntField & face = this->faceList[ iFace ];
        Mid<int> fMid( face.size(), iFace );
        fMid.id = iFace;
        fMid.data = face;
        std::sort( fMid.data.begin(), fMid.data.end() );
        this->refFaces.insert( fMid );
    }
}

void MyFaceSolver::SetBoundary()
{
    int nFace = this->face2Block.size();
    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        int bcType = domain_Machine.bctypeList[ iFace ];
        BlkF2C & face_struct = this->face2Block[ iFace ];
        face_struct.bctype = bcType;
    }
}

int MyFaceSolver::FindLineId( IntField & line )
{
    return this->FindId( line, this->refLines );
}

int MyFaceSolver::FindId( IntField & varlist, set< Mid<int> > &refSets )
{
    Mid<int> fMid( varlist.size(), 0 );
    fMid.data = varlist;
    std::sort( fMid.data.begin(), fMid.data.end() );

    set< Mid<int> >::iterator iter = refSets.find( fMid );
    if ( iter == refSets.end() )
    {
        return ONEFLOW::INVALID_INDEX;
    }
    return iter->id;
}

void MyFaceSolver::BuildSDomainList()
{
    int nFace = this->face2Block.size();
    this->sDomainList.resize( nFace );
    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        SDomain * sDomain = new SDomain();
        sDomain->domain_id = iFace;
        this->sDomainList[ iFace ] = sDomain;
    }

    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        IntField & lineList = this->faceList[ iFace ];
        IntField & posList = this->facePosList[ iFace ];

        SDomain * sDomain = this->sDomainList[ iFace ];
        sDomain->SetDomain( iFace, lineList, posList );
        sDomain->ConstructSDomainCtrlPoint();
        sDomain->ConstructDomainTopo();
        sDomain->CalcDim2D();
        sDomain->ConstructLocalTopoAsBlk2D();
    }
    int kkk = 1;
}

void MyFaceSolver::GenerateFaceMesh()
{
    int nFace = this->faceList.size();
    fstream file;
    OpenPrjFile( file, "grid/strtecplot.dat", ios_base::out );
    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        SDomain * sDomain = this->sDomainList[ iFace ];
        sDomain->Alloc();
        sDomain->SetDomainBcMesh();
        sDomain->GenerateSDomainMesh( file );
    }
    CloseFile( file );
}

void MyFaceSolver::GenerateLineMesh()
{
    int nLine = line_Machine.curveInfoList.size();
    slineList.resize( nLine );
    for ( int iSLine = 0; iSLine < nLine; ++ iSLine )
    {
        SLine * sLine = new SLine();
        slineList[ iSLine ] = sLine;
        sLine->line_id = iSLine + 1;
        sLine->ni = line_Machine.dimList[ iSLine ];
        sLine->Alloc();
        sLine->CopyMesh();
    }
}

BlkFaceSolver::BlkFaceSolver()
{
    this->flag = false;
}

BlkFaceSolver::~BlkFaceSolver()
{
    DeletePointer( blkList );
}

Face2D * BlkFaceSolver::GetBlkFace2D( int blk, int face_id )
{
    Block3D * blk3d = this->blkList[ blk ];
    int nFace = blk3d->facelist.size();
    for ( int i = 0; i < nFace; ++ i )
    {
        Face2D * face2d = blk3d->facelist[ i ];
        int fid = face2d->face_id;
        if ( fid == face_id )
        {
            return face2d;
        }
    }
    return 0;
}

int BlkFaceSolver::FindFace( Mid<int> & face )
{
    set< Mid<int> >::iterator iter = this->myFaceSolver.refFaces.find( face );
    if ( iter == this->myFaceSolver.refFaces.end() )
    {
        return ONEFLOW::INVALID_INDEX;
    }
    return iter->id;
}

int BlkFaceSolver::FindFaceId( IntField & face )
{
    Mid<int> fMid( face.size(), 0 );
    fMid.data = face;
    std::sort( fMid.data.begin(), fMid.data.end() );

    set< Mid<int> >::iterator iter = this->myFaceSolver.refFaces.find( fMid );
    if ( iter == this->myFaceSolver.refFaces.end() )
    {
        return ONEFLOW::INVALID_INDEX;
    }
    return iter->id;
}

void BlkFaceSolver::Alloc()
{
    if ( flag ) return;
    flag = true;
    myFaceSolver.CreateFaceList();
}

void BlkFaceSolver::AddLineToFace( int faceid, int pos, int lineid )
{
    this->myFaceSolver.AddLineToFace( faceid, pos, lineid );
}

void BlkFaceSolver::AddFace2Block( int blockid, int pos, int faceid )
{
    this->Alloc();
    this->blkset.insert( blockid );

    int fid = faceid - 1;
    IntField & face = this->myFaceSolver.faceList[ fid ];

    BlkF2C & face_struct = this->myFaceSolver.face2Block[ fid ];

    face_struct.id  = fid;
    face_struct.type = -1;
    face_struct.bctype = -1;
    face_struct.cellList.push_back( blockid );
    face_struct.posList.push_back( pos );
}

void BlkFaceSolver::SetBoundary()
{
    this->myFaceSolver.SetBoundary();
}

void BlkFaceSolver::BuildBlkFace()
{
    int nBlock = this->blkset.size();
    this->blkList.resize( nBlock );
    for ( int iBlk = 0; iBlk < nBlock; ++ iBlk )
    {
        Block3D * blk3d = new Block3D();
        this->blkList[ iBlk ] = blk3d;
    }

    int nFace = this->myFaceSolver.face2Block.size();
    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        BlkF2C & face_struct = this->myFaceSolver.face2Block[ iFace ];

        IntField & lineList = this->myFaceSolver.faceList[ iFace ];
        IntField & posList = this->myFaceSolver.facePosList[ iFace ];

        int n_neibor = face_struct.cellList.size();
        for ( int i = 0; i < n_neibor; ++ i )
        {
            int blk_id = face_struct.cellList[ i ] - 1;
            int pos = face_struct.posList[ i ];

            Block3D * blk3d = this->blkList[ blk_id ];
            blk3d->blk_id = blk_id;
            MDomain * mDomain = blk3d->mDomainList[ pos ];
            mDomain->AddSubDomain( iFace, lineList, posList );
        }
    }
}

void BlkFaceSolver::BuildSDomainList()
{
    this->myFaceSolver.BuildSDomainList();
}

void BlkFaceSolver::ConstructBlockInfo()
{
    int nBlock = this->blkList.size();

    for ( int iBlk = 0; iBlk < nBlock; ++ iBlk )
    {
        Block3D * blk3d = this->blkList[ iBlk ];
        blk3d->ConstructTopo();
    }

    for ( int iBlk = 0; iBlk < nBlock; ++ iBlk )
    {
        Block3D * blk3d = this->blkList[ iBlk ];
        blk3d->SetInterfaceBc();
    }
}

void BlkFaceSolver::DumpBcInp()
{
    int nBlock = this->blkList.size();
    int flowSolverIndex = 1;
    int width = 5;

    fstream file;
    OpenPrjFile( file, "grid/strplot3d.inp", ios_base::out );

    file << setw( width ) << flowSolverIndex << endl;
    file << setw( width ) << nBlock << endl;
    for ( int iBlk = 0; iBlk < nBlock; ++ iBlk )
    {
        Block3D * blk3d = this->blkList[ iBlk ];
        blk3d->DumpInp( file );
    }

    CloseFile( file );
    int kkk = 1;
}

void BlkFaceSolver::DumpBlkScript()
{
    fstream file;
    OpenPrjFile( file, "grid/blkscript.txt", ios_base::out );
    IntField ctrlpoints;
    ctrlpoints.push_back( 1 );
    ctrlpoints.push_back( 2 );
    ctrlpoints.push_back( 3 );
    ctrlpoints.push_back( 4 );
    ctrlpoints.push_back( 5 );
    ctrlpoints.push_back( 6 );
    ctrlpoints.push_back( 7 );
    ctrlpoints.push_back( 8 );
    bbElemHome.Init();
    BlkElem * blkHexa = bbElemHome.GetBlkElem( ONEFLOW::HEXA_8 );
    DumpBlkScript( file, blkHexa, ctrlpoints );
    CloseFile( file );
    int kkk = 1;
}

void BlkFaceSolver::DumpBlkScript( fstream & file, BlkElem * blkHexa, IntField & ctrlpoints )
{
    //int p1 = ctrlpoints[ 0 ];
    //int p2 = ctrlpoints[ 1 ];
    //file << "Line(1)=" << "{" << p1 << "," << p2 << "}\n";

    int nFace = blkHexa->faceList.size();
    for ( int i = 0; i < nFace; ++ i )
    {
        IntField localid = blkHexa->faceList[ i ];
        this->DumpBlkScript( file, localid, ctrlpoints );
    }
}

void BlkFaceSolver::DumpBlkScript( fstream & file, IntField & localid, IntField & ctrlpoints )
{
    int npoint = localid.size();
    for ( int i = 0; i < npoint; ++ i )
    {
        int i1 = localid[ i ];
        int i2 = localid[ ( i + 1 ) % npoint ];
        int p1 = ctrlpoints[ i1 ];
        int p2 = ctrlpoints[ i2 ];
        int id = i + 1;
        file << "Line(" << id << ")=" << "{" << p1 << "," << p2 << "}\n";
    }
}

void BlkFaceSolver::GenerateBlkMesh()
{
    int nBlock = this->blkList.size();
    for ( int iBlk = 0; iBlk < nBlock; ++ iBlk )
    {
        Block3D * blk3d = this->blkList[ iBlk ];
        blk3d->Alloc();
        blk3d->CreateBlockMesh();
    }
}

void BlkFaceSolver::GenerateFaceMesh()
{
    BuildSDomainList();
    myFaceSolver.GenerateFaceMesh();
}

void BlkFaceSolver::GenerateLineMesh()
{
    line_Machine.GenerateAllLineMesh();
    myFaceSolver.GenerateLineMesh();
}

void BlkFaceSolver::DumpStandardGrid()
{
    int nBlock = this->blkList.size();

    Grids strGridList( nBlock );
    strGridList.SetDeleteFlag( true );

    for ( int iBlk = 0; iBlk < nBlock; ++ iBlk )
    {
        Grid * gridstr = ONEFLOW::CreateGrid( ONEFLOW::SMESH );
        StrGrid * grid = ONEFLOW::StrGridCast( gridstr );

        strGridList[ iBlk ] = grid;

        Block3D * blk3d = this->blkList[ iBlk ];
        blk3d->FillStrGrid( grid, iBlk );

    }

    DumpStandardGrid( strGridList );
}

void BlkFaceSolver::DumpStandardGrid( Grids & strGridList )
{
    fstream file;
    OpenPrjFile( file, "grid/strplot3d.grd", ios_base::out | ios_base::binary );

    int nZone = strGridList.size();
    HXWrite( & file, nZone );
    for ( int iBlock = 0; iBlock < nZone; ++ iBlock )
    {
        Grid * gridstr = strGridList[ iBlock ];
        StrGrid * grid = ONEFLOW::StrGridCast( gridstr );
        int ni = grid->ni;
        int nj = grid->nj;
        int nk = grid->nk;
        HXWrite( & file, ni );
        HXWrite( & file, nj );
        HXWrite( & file, nk );
    }

    for ( int iBlock = 0; iBlock < nZone; ++ iBlock )
    {
        Grid * gridstr = strGridList[ iBlock ];
        StrGrid * grid = ONEFLOW::StrGridCast( gridstr );
        HXWrite( & file, grid->nodeMesh->xN );
        HXWrite( & file, grid->nodeMesh->yN );
        HXWrite( & file, grid->nodeMesh->zN );
    }

    CloseFile( file );

}

EndNameSpace