/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2023 He Xin and the OneFLOW contributors.
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
#include "HXType.h"
#include "HXDefine.h"
#include "HXCgns.h"
#include "metis.h"
#include <vector>
#include <fstream>


BeginNameSpace( ONEFLOW )

class RealList
{
public:
    RealList() ;
    ~RealList();
public:
    std::vector< Real > data;
public:
    size_t GetNElements();
    void AddData( Real value );
    Real operator [] ( int i ) const
    {
        return data[ i ];
    }

    Real & operator [] ( int i )
    {
        return data[ i ];
    }

    RealList & operator = ( Real const & value )
    {
        for ( int i = 0; i < data.size(); ++ i )
        {
            data[ i ] = value;
        }
        return *this;
    }

    void Resize( int new_size );
};

class IntList
{
public:
    IntList() ;
    ~IntList();
    IntList( const IntList & rhs );
public:
    std::vector< int > data;
public:
    size_t GetNElements();
    void AddData( int value );

    int operator [] ( int i ) const
    {
        return data[ i ];
    }

    int & operator [] ( int i )
    {
        return data[ i ];
    }

    void Resize( int new_size );

    void ReOrder( IntList & orderMap );
};

class EList
{
public:
    EList() ;
    ~EList();
public:
    std::vector< std::vector< int > > data;
public:
    size_t GetNElements();
    void AddElem( IntList &elem );
    void AddElem( std::vector< int > &elem );

    std::vector< int > & operator [] ( int i )
    {
        return data[ i ];
    }

    void ReOrder( IntList & orderMap );

    void Resize( int new_size );
};

class ScalarGrid;

class ScalarBcco
{
public:
    ScalarBcco();
    ~ScalarBcco();
public:
    std::string bcName;
    int bcType;

    EList elements;
    IntList eTypes;

    std::vector< int > vertexList;
    IntList local_globalIds;
public:
    void PushBoundaryFace( int pt, int eType );
    void AddBcPoint( int bcVertex );
    void ScanBcFace( ScalarGrid * grid );
    void ProcessVertexBc( IntSet & bcVertex );
};

class ScalarBccos
{
public:
    ScalarBccos();
    ~ScalarBccos();
public:
    std::vector< ScalarBcco * > bccos;
public:
    void AddBcco( ScalarBcco * scalarBcco );
    void ScanBcFace( ScalarGrid * grid );
};

class DataBase;
class DataBook;
class ScalarIFace;
class CgnsZbase;
class CgnsZone;
class SectionManager;

class ScalarGrid
{
public:
    ScalarGrid();
    ~ScalarGrid();
public:
    int nNodes, nCells, nBFaces, nFaces;
    int nTCells;
    RealList xn, yn, zn;
    RealList xfc, yfc, zfc;
    RealList xfn, yfn, zfn;
    RealList area;
    RealList xcc, ycc, zcc;
    RealList vol;

    IntList lc;
    IntList rc;
    IntList lpos;
    IntList rpos;

    std::vector<int> cell2faces;
    std::vector<int> c2fpos;

    //global faceid
    std::vector<int> global_faceid;

    EList faces;
    EList elements;
    //
    EList boundaryElements;
    IntList bcETypes;
    //
    IntList fTypes;
    IntList eTypes;
    IntList fBcTypes;
    IntList bcTypes;
    IntList bcNameIds;
    ScalarBccos * scalarBccos;
    DataBase * dataBase;
    ScalarIFace * scalarIFace;
    int type, level;
    int id;
    int localId;
    int volBcType;
public:
    DataBase * GetDataBase() { return dataBase; };
public:
    int GetNNodes();
    int GetNCells();
    int GetNFaces();
    int GetNBFaces();
    int GetNTCells();
    void GenerateGrid( int ni, Real xmin, Real xmax );
    void CalcTopology();
    void PushElement( int p1, int p2, int eType );
    void PushBoundaryFace( int pt, int eType );
    void ReorderFaces();
    void CalcOrderMap( IntList &orderMap );
    void SetBcGhostCell();
    void AllocGeom();
    void ScanBcFace();
    void ScanBcFace( IntSet& bcVertex, int bcType );
    bool CheckBcFace( IntSet & bcVertex, std::vector< int > & nodeId );
    void AllocateBc();
    void SetBcTypes();
public:
    void ReadFromCgnsZbase( CgnsZbase * cgnsZbase );
    void ReadFromCgnsZone( CgnsZone * cgnsZone );
    void PushElement( CgIntField & eNodeId, int eType );
    void GenerateGridFromCgns( const std::string & prjFileName );
    void DumpCgnsGrid();
    void SetCgnsZone( CgnsZone * cgnsZone );
public:
    void CalcVolumeSection( SectionManager * volumeSectionManager );
    void CalcBoundarySection( SectionManager * bcSectionManager );
public:
    void CalcMetrics1D();
    void CalcFaceCenter1D();
    void CalcCellCenter1D();
    void CalcCellCenterVol1D();
    void CalcCellVolume1D();
    void CalcFaceNormal1D();
    void CalcGhostCellCenterVol1D();
public:
    void CalcC2C( EList & c2c );
    void CalcInterfaceToBcFace();
    void Normalize();
public:
    void GetSId( int i_interface, int & sId );
    void GetTId( int i_interface, int & tId );
public:
    void DumpCalcGrid();
    void ReadCalcGrid();
    void WriteGrid( DataBook * databook );
    void ReadGrid( DataBook * databook );
    void WriteGrid( std::fstream & file );
    void ReadGrid( std::fstream & file );
    void CreateNodes( int numberOfNodes );
    void WriteGridFaceTopology( DataBook * databook );
    void WriteBoundaryTopology( DataBook * databook );
    void ReadGridFaceTopology( DataBook * databook );
    void ReadBoundaryTopology( DataBook * databook );
    void NormalizeBc();
public:
    //partition
    void AddFaceType( int fType );
    void AddInterface( int global_interface_id, int neighbor_zoneid, int neighbor_cellid );
    void AddPhysicalBcFace( int global_face_id, int bctype, int lcell, int rcell );
    void AddInnerFace( int global_face_id, int bctype, int lcell, int rcell );
    void AddInterfaceBcFace( int global_face_id, int bctype, int lcell, int rcell, int nei_zoneid, int nei_cellid );
    void ReconstructNode( ScalarGrid * ggrid );

};


EndNameSpace
