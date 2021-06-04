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
#include "Configure.h"
#include "HXType.h"
#include "HXDefine.h"
#include "metis.h"
#include <vector>
using namespace std;

BeginNameSpace( ONEFLOW )

class RealList
{
public:
    RealList() ;
    ~RealList();
public:
    vector< Real > data;
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
    vector< int > data;
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
    vector< vector< int > > data;
public:
    size_t GetNElements();
    void AddElem( IntList &elem );
    void AddElem( vector< int > &elem );

    vector< int > & operator [] ( int i )
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
    int bcType;
    vector< int > vertexList;
public:
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
    vector< ScalarBcco * > bccos;
public:
    void AddBcco( ScalarBcco * scalarBcco );
    void ScanBcFace( ScalarGrid * grid );
};

class DataBase;
class GridTopo;
class DataBook;
class ScalarIFace;

class ScalarGrid
{
public:
    ScalarGrid();
    ScalarGrid( int grid_id );
    ~ScalarGrid();
public:
    size_t nNodes, nCells, nBFaces, nFaces;
    size_t nTCells;
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

    EList faces;
    EList elements;
    IntList fTypes;
    IntList eTypes;
    IntList fBcTypes;
    IntList bcTypes;
    IntList bcNameIds;
    ScalarBccos * scalarBccos;
    DataBase * dataBase;
    GridTopo * gridTopo;
    ScalarIFace * scalarIFace;
    int grid_id;
    int volBcType;
public:
    DataBase * GetDataBase() { return dataBase; };
public:
    size_t GetNNodes();
    size_t GetNCells();
    size_t GetNFaces();
    size_t GetNBFaces();
    size_t GetNTCells();
    void GenerateGrid( int ni, Real xmin, Real xmax );
    void CalcTopology();
    void PushElement( int p1, int p2, int eType );
    void ReorderFaces();
    void CalcOrderMap( IntList &orderMap );
    void SetBcGhostCell();
    void AllocGeom();
    void ScanBcFace();
    void ScanBcFace( IntSet& bcVertex, int bcType );
    bool CheckBcFace( IntSet & bcVertex, vector< int > & nodeId );
    void AllocateBc();
    void SetBcTypes();
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
public:
    void GetSId( int i_interface, int & sId );
    void GetTId( int i_interface, int & tId );
public:
    void DumpCalcGrid();
    void ReadCalcGrid();
    void WriteGrid( DataBook * databook );
    void ReadGrid( DataBook * databook );
    void CreateNodes( int numberOfNodes );
    void WriteGridFaceTopology( DataBook * databook );
    void WriteBoundaryTopology( DataBook * databook );
};


EndNameSpace