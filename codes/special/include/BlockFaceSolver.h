/*---------------------------------------------------------------------------*\\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\\*---------------------------------------------------------------------------*/


#pragma once
#include "HXDefine.h"
#include "HXLookup.h"
#include "CalcCoor.h"
#include "SimpleDomain.h"
#include "GridHandles.h"
#include <set>
#include <map>
#include <fstream>
#include <memory>

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
    HXVector< std::unique_ptr< Block3D > > blkList;
    HXVector< std::unique_ptr< Block2D > > blkList2d;
    bool flag;
public:
    bool init_flag;
    LinkField faceList;
    LinkField faceLinePosList;
    HXLookup<int> lineLookup;
    HXLookup<int> faceLookup;
    IntSet faceset;
public:
    void Reset();
    const Face2D * GetBlkFace( int blk, int face_id ) const;
    const Face2D * GetBlkFace2D( int blk, int face_id ) const;
    IntField & GetLine( int line_id );
    const IntField & GetLine( int line_id ) const;
    BlkF2C & GetLineToFace( int line_id );
    const BlkF2C & GetLineToFace( int line_id ) const;
    BlkF2C & GetFaceToBlock( int faceIndex );
    const BlkF2C & GetFaceToBlock( int faceIndex ) const;
    SDomain * GetSDomain( int domainIndex );
    SLine * GetSLine( int lineIndex );
    int FindLineId( const IntField & line ) const;
public:
    void AddLineToFace( int faceid, int pos, int lineid );
    void AddFace2Block( int blockid, int pos, int faceid );
    void GenerateGrid();

private:
    HXVector< std::unique_ptr< SDomain > > sDomainList;
    HXVector< std::unique_ptr< SLine > > slineList;
    LinkField lineList;
    HXVector< BlkF2C > line2Face;
    HXVector< BlkF2C > face2Block;

    void Alloc();
    void InitializeLineTopology();
    void CreateFaceList();
    void BuildSurfaceDomainList();
    void GenerateSurfaceFaceMesh();
    void GenerateSurfaceLineMesh();
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
    void DumpBlkScript();
    void DumpBlkScript( std::fstream & file, BlkElem * blkHexa, IntField & ctrlpoints );
    void DumpBlkScript( std::fstream & file, IntField & localid, IntField & ctrlpoints );
};

extern BlkFaceSolver blkFaceSolver;


EndNameSpace
