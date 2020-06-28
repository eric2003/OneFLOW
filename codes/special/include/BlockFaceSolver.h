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

class MyFaceSolver
{
public:
    MyFaceSolver();
    ~MyFaceSolver();
public:
    bool init_flag;
    LinkField lineList; 
    LinkField faceList;
    LinkField facePosList;
    set< Mid< int > > refLines;
    set< Mid< int > > refFaces;
    IntSet faceset;
    HXVector< BlkF2C > line2Face;
    HXVector< BlkF2C > face2Block;
    HXVector< SDomain * > sDomainList;
    HXVector< SLine * > slineList;
public:
    void Alloc();
    void AddLineToFace( int faceid, int pos, int lineid );
    void CreateFaceList();
    void SetBoundary();
    int FindLineId( IntField & line );
    int FindId( IntField & varlist, set< Mid<int> > &refSets );
    IntField & GetLine( int line_id );
public:
    void BuildSDomainList();
    void GenerateFaceMesh();
    void GenerateLineMesh();
};

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
    HXVector< Block2D * > blk2dList;
    bool flag;
    MyFaceSolver myFaceSolver;
public:
    Face2D * GetBlkFace2D( int blk, int face_id );
    int  FindFace( Mid<int> & face );
    int  FindFaceId( IntField & face );
public:
    void Alloc();
    void AddLineToFace( int faceid, int pos, int lineid );
    void AddFace2Block( int blockid, int pos, int faceid );
    void BuildBlkFace();
    void SetBoundary();
    void DumpBcInp();
    void ConstructBlockInfo();
    void GenerateBlkMesh();
    void GenerateFaceMesh();
    void GenerateLineMesh();
    void BuildSDomainList();
    void DumpStandardGrid();
    void DumpStandardGrid( Grids & strGridList );
public:
    void DumpBlkScript();
    void DumpBlkScript( fstream & file, BlkElem * blkHexa, IntField & ctrlpoints );
    void DumpBlkScript( fstream & file, IntField & localid, IntField & ctrlpoints );
};

extern BlkFaceSolver blkFaceSolver;

EndNameSpace