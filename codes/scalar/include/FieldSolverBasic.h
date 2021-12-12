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
#include "ScalarGrid.h"
#include "HXArray.h"
#include "Task.h"
#include <vector>
#include <string>


BeginNameSpace( ONEFLOW )

class ScalarField;
class ScalarGrid;
class FieldPara;
class ScalarFieldManager;
class ScalarFieldRecord;

template < typename T >
class SortArray
{
public:
    std::vector< T > data;
    bool operator < ( const SortArray< T > & rhs ) const
    {
        return this->data[ 0 ] < rhs.data[ 0 ];
    }
};

class DataBook;

class FieldSolverBasic
{
public:
    FieldSolverBasic();
    ~FieldSolverBasic();
public:
    ScalarField * field;
    ScalarGrid * grid;
    FieldPara * para;
    ScalarFieldManager * scalarFieldManager;
    std::vector< ScalarField * > fields;
    std::vector< ScalarGrid * > grids;
    bool tmpflag_delete_grids;
public:
    //tmp
    void FillTmpGridVector();
public:
    virtual void Run();
    void Init();
    void LoadGrid();
    void InitCtrlParameter();
    void AddZoneGrid();
    void CalcGridMetrics();
    void InitFlowField();
    void InitFlowField_Basic();
    void CommParallelInfo();
    void UploadInterface();
    void DownloadInterface();
    void UpdateInterface( TaskFunction sendAction, TaskFunction recvAction );
    void SwapInterfaceData( int iZone, int jZone, TaskFunction sendAction, TaskFunction recvAction );
public:
    void Visualize();
    void ToTecplot( RealField & xList, RealField & varlist, std::string const & fileName );
    void Theory( ScalarGrid * grid, Real time, RealField & theory );
    void GetVisualData( DataBook * dataBook );
    void AddVisualData( RealField & qList, RealField & theoryList, RealField & xcoorList );
    void AddVisualData( DataBook * dataBook, RealField & qList, RealField & theoryList, RealField & xcoorList );
    void Reorder( RealField & a, RealField & b, RealField & c );
public:
    Real ScalarFun( Real xm );
    Real SquareFun( Real xm );
};

void PrepareFieldSendData();
void PrepareFieldRecvData();
ScalarFieldRecord * PrepareSendScalarFieldRecord();
ScalarFieldRecord * PrepareRecvScalarFieldRecord();

void PrepareGeomSendData();
void PrepareGeomRecvData();

void TestMPI();

EndNameSpace
