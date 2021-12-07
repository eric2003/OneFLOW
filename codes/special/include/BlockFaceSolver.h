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


#pragma once
#include "HXDefine.h"
#include "Mid.h"
#include "CalcCoor.h"
#include "SimpleDomain.h"
#include "GridDef.h"
#include <set>
#include <map>
#include <fstream>
using namespace std;
BeginNameSpace( ONEFLOW )

class Block3D;
class BlkMesh;
class SDomain;
class MDomain;
class CalcCoor;
class Face2D;
class MLine;
class SLine;
class Block2D;
class BlkElem;

class BlkFaceSolver
{
public:
    BlkFaceSolver();
    ~BlkFaceSolver();
public:
    IntSet blkset;
    HXVector< Block3D * > blkList;
    HXVector< Block2D * > blkList2d;
    bool flag;
public:
    bool init_flag;
    LinkField lineList; 
    LinkField faceList;
    LinkField faceLinePosList;
    std::set< Mid< int > > refLines;
    std::set< Mid< int > > refFaces;
    IntSet faceset;
    HXVector< BlkF2C > line2Face;
    HXVector< BlkF2C > face2Block;
    HXVector< SDomain * > sDomainList;
    HXVector< SLine * > slineList;
public:
    Face2D * GetBlkFace( int blk, int face_id );
    Face2D * GetBlkFace2D( int blk, int face_id );
public:
    void Alloc();
    void MyFaceAlloc();
    void CreateFaceList();
    int  FindLineId( IntField & line );
    int  FindId( IntField & varlist, std::set< Mid<int> > &refSets );
    int  FindFace( Mid<int> & face );
    int  FindFaceId( IntField & face );
    IntField & GetLine( int line_id );
    void MyFaceBuildSDomainList();
    void MyFaceGenerateFaceMesh();
    void MyFaceGenerateLineMesh();
public:
    void AddLineToFace( int faceid, int pos, int lineid );
    void AddFace2Block( int blockid, int pos, int faceid );
    void BuildBlkFace();
    void BuildBlkFace2D();
    void SetBoundary();
    void DumpBcInp();
    void DumpBcInp2D();
    void ConstructBlockInfo();
    void ConstructBlockInfo2D();
    void GenerateBlkMesh();
    void GenerateBlkMesh2D();
    void GenerateFaceMesh();
    void GenerateLineMesh();
    void DumpStandardGrid();
    void DumpStandardGrid2D();
    void DumpStandardGrid( Grids & strGridList );
    void GenerateFaceBlockLink();
public:
    void DumpBlkScript();
    void DumpBlkScript( std::fstream & file, BlkElem * blkHexa, IntField & ctrlpoints );
    void DumpBlkScript( std::fstream & file, IntField & localid, IntField & ctrlpoints );
};

extern BlkFaceSolver blkFaceSolver;

IntField GlobalGetLine( int line_id );

EndNameSpace
