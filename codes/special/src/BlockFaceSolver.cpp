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

#include "BlockFaceSolver.h"
#include "DomainMachine.h"
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
#include "Block3D.h"
#include "Block2D.h"
#include "HXPointer.h"
#include "CurveInfo.h"
#include "SegmentCtrl.h"
#include "BlockElem.h"
#include "FileUtil.h"
#include "Prj.h"
#include "HXCgns.h"
#include "Dimension.h"
#include "GridPara.h"
#include <algorithm>
#include <iostream>


BeginNameSpace( ONEFLOW )

BlkFaceSolver blkFaceSolver;


BlkFaceSolver::BlkFaceSolver()
{
    this->flag = false;
    this->init_flag = false;
}

BlkFaceSolver::~BlkFaceSolver()
{
    DeletePointer( blkList );
    DeletePointer( sDomainList );
    DeletePointer( slineList );
}

Face2D * BlkFaceSolver::GetBlkFace( int blk, int face_id )
{
    Block3D * blk3d = this->blkList[ blk ];
    int nFaces = blk3d->facelist.size();
    for ( int i = 0; i < nFaces; ++ i )
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

Face2D * BlkFaceSolver::GetBlkFace2D( int blk, int face_id )
{
    Block2D * blk2d = this->blkList2d[ blk ];
    int nFaces = blk2d->facelist.size();
    for ( int i = 0; i < nFaces; ++ i )
    {
        Face2D * face2d = blk2d->facelist[ i ];
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
    std::set< Mid<int> >::iterator iter = this->refFaces.find( face );
    if ( iter == this->refFaces.end() )
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

    std::set< Mid<int> >::iterator iter = this->refFaces.find( fMid );
    if ( iter == this->refFaces.end() )
    {
        return ONEFLOW::INVALID_INDEX;
    }
    return iter->id;
}

void BlkFaceSolver::MyFaceBuildSDomainList()
{
    int nFaces = this->face2Block.size();
    this->sDomainList.resize( nFaces );
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        SDomain * sDomain = new SDomain();
        sDomain->domain_id = iFace;
        this->sDomainList[ iFace ] = sDomain;
    }

    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        IntField & lineList = this->faceList[ iFace ];
        IntField & posList = this->faceLinePosList[ iFace ];

        SDomain * sDomain = this->sDomainList[ iFace ];
        sDomain->SetDomain( iFace, lineList, posList );
        sDomain->ConstructSDomainCtrlPoint();
        sDomain->ConstructDomainTopo();
        sDomain->CalcDim2D();
        sDomain->ConstructLocalTopoAsBlk2D();
    }
    int kkk = 1;
}

void BlkFaceSolver::MyFaceGenerateFaceMesh()
{
    int nFaces = this->faceList.size();
    std::fstream file;
    Prj::OpenPrjFile( file, "grid/facemesh_tecplot.dat", std::ios_base::out );
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        SDomain * sDomain = this->sDomainList[ iFace ];
        sDomain->Alloc();
        sDomain->SetDomainBcMesh();
        sDomain->GenerateSDomainMesh( file );
    }
    CloseFile( file );
}

void BlkFaceSolver::MyFaceGenerateLineMesh()
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

void BlkFaceSolver::Alloc()
{
    if ( flag ) return;
    flag = true;
    this->CreateFaceList();
}

void BlkFaceSolver::CreateFaceList()
{
    //All cell faces have been built in faceset
    int nFaces = faceset.size();
    this->faceList.resize( nFaces );
    this->faceLinePosList.resize( nFaces );

    int nLine = this->line2Face.size();
    for ( int iLine = 0; iLine < nLine; ++ iLine )
    {
        BlkF2C & line_struct = this->line2Face[ iLine ];
        int n = line_struct.cellList.size();
        for ( int i = 0; i < n; ++ i )
        {
            int face_id = line_struct.cellList[ i ] - 1;
            //line position in face
            int face_line_pos = line_struct.posList[ i ];
            faceList[ face_id ].push_back( line_struct.id );
            this->faceLinePosList[ face_id ].push_back( face_line_pos );
        }
    }

    this->face2Block.resize( nFaces );
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        IntField & face = this->faceList[ iFace ];
        Mid<int> fMid( face.size(), iFace );
        fMid.id = iFace;
        fMid.data = face;
        std::sort( fMid.data.begin(), fMid.data.end() );
        this->refFaces.insert( fMid );
    }
}

IntField & BlkFaceSolver::GetLine( int line_id )
{
    int id = line_id - 1;
    return lineList[ id ];
}

int BlkFaceSolver::FindLineId( IntField & line )
{
    return this->FindId( line, this->refLines );
}

int BlkFaceSolver::FindId( IntField & varlist, std::set< Mid<int> > &refSets )
{
    Mid<int> fMid( varlist.size(), 0 );
    fMid.data = varlist;
    std::sort( fMid.data.begin(), fMid.data.end() );

    std::set< Mid<int> >::iterator iter = refSets.find( fMid );
    if ( iter == refSets.end() )
    {
        return ONEFLOW::INVALID_INDEX;
    }
    return iter->id;
}

void BlkFaceSolver::MyFaceAlloc()
{
    if ( init_flag ) return;
    init_flag = true;
    int nLine = line_Machine.curveInfoList.size();
    for ( int i = 0; i < nLine; ++ i )
    {
        CurveInfo * curveInfo = line_Machine.GetCurveInfo( i + 1 );
        IntField line;
        line.push_back( curveInfo->p1 );
        line.push_back( curveInfo->p2 );
        this->lineList.push_back( line );
    }
    this->line2Face.resize( nLine );

    //The reference line segment std::set is std::set up to facilitate the search
    for ( int i = 0; i < nLine; ++ i )
    {
        IntField & line = this->lineList[ i ];
        Mid<int> lMid( line.size(), i );
        lMid.id = i;
        lMid.data = line;
        std::sort( lMid.data.begin(), lMid.data.end() );
        this->refLines.insert( lMid );
    }
}


void BlkFaceSolver::AddLineToFace( int faceid, int pos, int lineid )
{
    this->MyFaceAlloc();

    int id = lineid - 1;
    BlkF2C & line_struct = this->line2Face[ id ];

    int bctype = domain_Machine.GetBcType( lineid );

    line_struct.id  = lineid;
    line_struct.type = -1;
    line_struct.bctype = bctype;
    line_struct.cellList.push_back( faceid );
    line_struct.posList.push_back( pos );

    faceset.insert( faceid );
}

void BlkFaceSolver::AddFace2Block( int blockid, int pos, int faceid )
{
    //Set all the cell faces
    this->Alloc();
    this->blkset.insert( blockid );

    int fid = faceid - 1;
    IntField & face = this->faceList[ fid ];

    BlkF2C & face_struct = this->face2Block[ fid ];

    face_struct.id  = fid;
    face_struct.type = -1;
    face_struct.bctype = -1;
    face_struct.cellList.push_back( blockid );
    face_struct.posList.push_back( pos );
}

void BlkFaceSolver::SetBoundary()
{
    int nFaces = this->face2Block.size();
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        int bcType = domain_Machine.bctypeList[ iFace ];
        BlkF2C & face_struct = this->face2Block[ iFace ];
        face_struct.bctype = bcType;
    }
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

    int nFaces = this->face2Block.size();
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        BlkF2C & face_struct = this->face2Block[ iFace ];

        IntField & lineList = this->faceList[ iFace ];
        IntField & lineposList = this->faceLinePosList[ iFace ];

        int n_neibor = face_struct.cellList.size();
        for ( int i = 0; i < n_neibor; ++ i )
        {
            int blk_id = face_struct.cellList[ i ] - 1;
            int face_pos_in_blk = face_struct.posList[ i ];

            Block3D * blk3d = this->blkList[ blk_id ];
            blk3d->blk_id = blk_id;
            MDomain * mDomain = blk3d->mDomainList[ face_pos_in_blk ];
            mDomain->AddSubDomain( iFace, lineList, lineposList );
        }
    }
}

void BlkFaceSolver::BuildBlkFace2D()
{
    int nBlock = this->blkset.size();
    this->blkList2d.resize( nBlock );
    for ( int iBlk = 0; iBlk < nBlock; ++ iBlk )
    {
        Block2D * blk2d = new Block2D();
        this->blkList2d[ iBlk ] = blk2d;
    }

    int nFaces = this->face2Block.size();
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        BlkF2C & face_struct = this->face2Block[ iFace ];

        IntField & lineList = this->faceList[ iFace ];
        IntField & lineposList = this->faceLinePosList[ iFace ];

        int n_neibor = face_struct.cellList.size();
        for ( int i = 0; i < n_neibor; ++ i )
        {
            int blk_id = face_struct.cellList[ i ] - 1;
            int face_pos_in_blk = face_struct.posList[ i ] - 1;

            Block2D * blk2d = this->blkList2d[ blk_id ];
            blk2d->blk_id = blk_id;

            MDomain * mDomain = blk2d->mDomainList[ face_pos_in_blk ];
            mDomain->AddSubDomain( iFace, lineList, lineposList );
        }
    }
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

void BlkFaceSolver::ConstructBlockInfo2D()
{
    int nBlock = this->blkList2d.size();

    for ( int iBlk = 0; iBlk < nBlock; ++ iBlk )
    {
        Block2D * blk2d = this->blkList2d[ iBlk ];
        blk2d->ConstructTopo();
    }

    for ( int iBlk = 0; iBlk < nBlock; ++ iBlk )
    {
        Block2D * blk2d = this->blkList2d[ iBlk ];
        blk2d->SetInterfaceBc();
    }
}

void BlkFaceSolver::DumpBcInp()
{
    int nBlock = this->blkList.size();
    int flowSolverIndex = 1;
    int width = 5;

    std::fstream file;
    Prj::OpenPrjFile( file, grid_para.bcFile, std::ios_base::out );

    file << std::setw( width ) << flowSolverIndex << std::endl;
    file << std::setw( width ) << nBlock << std::endl;
    for ( int iBlk = 0; iBlk < nBlock; ++ iBlk )
    {
        Block3D * blk3d = this->blkList[ iBlk ];
        blk3d->DumpInp( file );
    }

    CloseFile( file );
    int kkk = 1;
}

void BlkFaceSolver::DumpBcInp2D()
{
    int nBlock = this->blkList2d.size();
    int flowSolverIndex = 1;
    int width = 5;

    std::fstream file;
    Prj::OpenPrjFile( file, grid_para.bcFile, std::ios_base::out );

    file << std::setw( width ) << flowSolverIndex << std::endl;
    file << std::setw( width ) << nBlock << std::endl;
    for ( int iBlk = 0; iBlk < nBlock; ++ iBlk )
    {
        Block2D * blk2d = this->blkList2d[ iBlk ];
        blk2d->DumpInp( file );
    }

    CloseFile( file );
    int kkk = 1;
}

void BlkFaceSolver::DumpBlkScript()
{
    std::fstream file;
    Prj::OpenPrjFile( file, "grid/blkscript.txt", std::ios_base::out );
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

void BlkFaceSolver::DumpBlkScript( std::fstream & file, BlkElem * blkHexa, IntField & ctrlpoints )
{
    //int p1 = ctrlpoints[ 0 ];
    //int p2 = ctrlpoints[ 1 ];
    //file << "Line(1)=" << "{" << p1 << "," << p2 << "}\n";

    int nFaces = blkHexa->faceList.size();
    for ( int i = 0; i < nFaces; ++ i )
    {
        IntField localid = blkHexa->faceList[ i ];
        this->DumpBlkScript( file, localid, ctrlpoints );
    }
}

void BlkFaceSolver::DumpBlkScript( std::fstream & file, IntField & localid, IntField & ctrlpoints )
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

void BlkFaceSolver::GenerateBlkMesh2D()
{
    int nBlock = this->blkList2d.size();
    for ( int iBlk = 0; iBlk < nBlock; ++ iBlk )
    {
        Block2D * blk2d = this->blkList2d[ iBlk ];
        blk2d->Alloc();
        blk2d->CreateBlockMesh2D();
    }

    std::fstream file;
    Prj::OpenPrjFile( file, "grid/blkplot2d.dat", std::ios_base::out );

    for ( int iBlk = 0; iBlk < nBlock; ++ iBlk )
    {
        Block2D * blk2d = this->blkList2d[ iBlk ];
        blk2d->DumpBlockMesh2D( file );
    }
    CloseFile( file );
}

void BlkFaceSolver::GenerateFaceMesh()
{
    this->MyFaceBuildSDomainList();
    this->MyFaceGenerateFaceMesh();
}

void BlkFaceSolver::GenerateLineMesh()
{
    line_Machine.GenerateAllLineMesh();
    this->MyFaceGenerateLineMesh();
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


void BlkFaceSolver::DumpStandardGrid2D()
{
    int nBlock = this->blkList2d.size();

    Grids strGridList( nBlock );
    strGridList.SetDeleteFlag( true );

    for ( int iBlk = 0; iBlk < nBlock; ++ iBlk )
    {
        Grid * gridstr = ONEFLOW::CreateGrid( ONEFLOW::SMESH );
        StrGrid * grid = ONEFLOW::StrGridCast( gridstr );

        strGridList[ iBlk ] = grid;

        Block2D * blk2d = this->blkList2d[ iBlk ];
        blk2d->FillStrGrid( grid, iBlk );

    }

    DumpStandardGrid( strGridList );
}

void BlkFaceSolver::DumpStandardGrid( Grids & strGridList )
{
    std::fstream file;
    Prj::OpenPrjFile( file, grid_para.gridFile, std::ios_base::out | std::ios_base::binary );

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

void BlkFaceSolver::GenerateFaceBlockLink()
{
    if ( Dim::dimension == ONEFLOW::THREE_D )
    {
        this->DumpBlkScript();
        this->SetBoundary();
        this->BuildBlkFace();
        this->ConstructBlockInfo();
        this->DumpBcInp();
        this->GenerateLineMesh();
        this->GenerateFaceMesh();
        this->GenerateBlkMesh();
        this->DumpStandardGrid();
    }
    else
    {
        this->SetBoundary();
        this->BuildBlkFace2D();
        this->ConstructBlockInfo2D();
        this->DumpBcInp2D();
        this->GenerateLineMesh();
        this->GenerateFaceMesh();
        this->GenerateBlkMesh2D();
        this->DumpStandardGrid2D();
    }
}

IntField GlobalGetLine( int line_id )
{
    return blkFaceSolver.GetLine( line_id );
}

EndNameSpace
