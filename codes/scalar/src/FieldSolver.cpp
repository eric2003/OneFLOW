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

#include "FieldSolver.h"
#include "FieldPara.h"
#include "DataBase.h"
#include "ScalarDataIO.h"
#include "ScalarGrid.h"
#include "MetisGrid.h"
#include "ScalarField.h"
#include "ScalarIFace.h"
#include "ZoneState.h"
#include "ActionState.h"
#include "GridState.h"
#include "ScalarFieldRecord.h"
#include "ScalarAlloc.h"
#include "SolverDef.h"
#include "Prj.h"
#include "FileUtil.h"
#include "Parallel.h"
#include "ScalarZone.h"
#include "HXCgns.h"
#include "HXMath.h"
#include <iostream>
#include <vector>
using namespace std;

BeginNameSpace( ONEFLOW )

FieldSolver::FieldSolver()
{
    this->grid = new ScalarGrid();
    this->field = new ScalarField();
    this->para = new FieldPara();
    this->tmpflag_delete_grids = true;
    this->scalarFieldManager = new ScalarFieldManager();
    ScalarZone::Allocate();
}

FieldSolver::~FieldSolver()
{
    delete this->grid;
    delete this->field;
    delete this->para;
    ScalarZone::DeAllocate();
    delete this->scalarFieldManager;
    if ( tmpflag_delete_grids )
    {
        for ( int i = 0; i < grids.size(); ++ i )
        {
            delete grids[ i ];
        }
    }

    for ( int i = 0; i < fields.size(); ++ i )
    {
        delete fields[ i ];
    }
}

void FieldSolver::Run()
{
    this->Init();

    this->SolveFlowField();
}

void FieldSolver::Init()
{
    this->InitCtrlParameter();

    this->InitGrid();

    this->InitFlowField();

    this->InitParallelInfo();

    this->CommParallelInfo();
}

void FieldSolver::InitCtrlParameter()
{
    this->para->Init();
}

void FieldSolver::InitGrid()
{
    this->grid->GenerateGrid( this->para->nx, 0, this->para->len );
    this->grid->CalcTopology();
    this->grid->CalcMetrics1D();

    Part part;
    //part.PartitionGrid( this->grid, 4, & this->grids );
    int npart = this->para->nx - 1;
    //int npart = 32;
    cout << " npart = " << npart << "\n";
    //part.PartitionGrid( this->grid, 4, & this->grids );
    part.PartitionGrid( this->grid, npart, & this->grids );

    int nZones = this->grids.size();
    ZoneState::nZones = nZones;
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        ScalarZone::AddGrid( iZone, this->grids[ iZone ] );
    }
    int kkk = 1;
}

//void FieldSolver::InitGrid()
//{
//    this->grid->GenerateGrid( this->para->nx, 0, this->para->len );
//    this->grid->CalcTopology();
//    this->grid->CalcMetrics1D();
//
//    this->tmpflag_delete_grids = false;
//
//    this->grids.push_back( this->grid );
//}

void FieldSolver::InitFlowField()
{
    this->scalarFieldManager->Init();
    int nZones = this->grids.size();
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        ZoneState::zid = iZone;
        cout << " iZone = " << iZone << "\n";
        this->scalarFieldManager->AllocateAllFields();
    }

    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        ScalarField * field = new ScalarField();
        this->fields.push_back( field );
    }

    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        ZoneState::zid = iZone;
        this->InitFlowField_Basic();
    }
}

void FieldSolver::InitFlowField_Basic()
{
    ScalarGrid * grid = ScalarZone::GetGrid();

    RealField & q   = GetFieldReference< MRField > ( grid, "q" ).AsOneD();

    int nTCells = grid->GetNTCells();

    for ( int iCell = 0; iCell < nTCells; ++ iCell )
    {
        Real xm = grid->xcc[ iCell ];
        q[ iCell ] = this->ScalarFun( xm );
    }
}

Real FieldSolver::ScalarFun( Real xm )
{
    return SquareFun( xm );
}

Real FieldSolver::SquareFun( Real xm )
{
    if ( xm >= 0.5 && xm <= 1.0 )
    {
        return 2.0;
    }
    return 1.0;
}

void FieldSolver::InitParallelInfo()
{
    int nZones = this->grids.size();
    ZoneState::nZones = nZones;
    ZoneState::pid.resize( nZones );
    ZoneState::zoneType.resize( nZones );

    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        ZoneState::pid[ iZone ] = iZone % Parallel::nProc;
        ZoneState::zoneType[ iZone ] = ONEFLOW::UMESH;
    }
}

void FieldSolver::UploadInterface()
{
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        ZoneState::zid = iZone;

        this->scalarFieldManager->UploadInterfaceField();
    }
}

void FieldSolver::DownloadInterface()
{
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        ZoneState::zid = iZone;

        this->scalarFieldManager->DownloadInterfaceField();
    }
}

void FieldSolver::UpdateInterface( TaskFunction sendAction, TaskFunction recvAction )
{
    ActionState::dataBook = new DataBook();
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        //Loop through each zone
        ZoneState::zid = iZone;
        //Find out the neighbors of each zone 
        //(that is, the blocks that have docking relationship and need to transfer interface information)

        //int nNei = interFaceTopo.data[ iZone ].size();
        ScalarGrid * grid = ScalarZone::GetGrid();
        GridTopo * gridTopo = grid->gridTopo;
        ScalarIFace * scalarIFace = gridTopo->scalarIFace;
        int nNei = scalarIFace->data.size();

        //For all neighbors of this block (zone = iZone), exchange information
        ZoneState::zid  = iZone;
        for ( int iNei = 0; iNei < nNei; ++ iNei )
        {
            //jZone is the block number of the neighbor block
            //int jZone = interFaceTopo.data[ iZone ][ iNei ];
            int jZone = scalarIFace->data[ iNei ].zonej;
            ZoneState::inei = iNei;
            //Exchange data between this block (iZone) and its neighbor (jZone)
            //But it is not necessarily in this process, because the process of block Zid
            //and the process of block jZone may not be in this process
            this->SwapInterfaceData( iZone, jZone, sendAction, recvAction );
        }
    }
    delete ActionState::dataBook;
}

void FieldSolver::SwapInterfaceData( int iZone, int jZone, TaskFunction sendAction, TaskFunction recvAction )
{
    int sPid = ZoneState::pid[ iZone ];
    int rPid = ZoneState::pid[ jZone ];
    //The process ID of this block (block number iZone) is sPid,
    //and the process of the target block (jZone, the neighbor of iZone) is rpid
    //If the current zone is in this process, it means that this zone needs to send (send before receive)
    if ( Parallel::pid == sPid )
    {
        ZoneState::zid  = iZone;
        ZoneState::rzid = jZone;

        //This sendaction is not a real send, it just packages the data that needs to be sent into actionstate:: databook,
        //or it is similar to delivering express to SF
        sendAction();
    }
    //Hxswapdata is the place where these data are really delivered.
    //Here, after the package is delivered to SF express, SF express will send the package (data) to the specified address (process)

    HXSwapData( ActionState::dataBook, sPid, rPid );

    //If the current process is not on this process, 
    //then the only action that this process can take is to participate in receiving data as the receiver
    if ( Parallel::pid == rPid )
    {
        ZoneState::zid  = jZone;
        ZoneState::szid = iZone;
        //The zoneid of the receiving action block is actually jzone

        recvAction();
    }
}

//void FieldSolver::PrepareSendData()
//{
//    ScalarFieldRecord * fieldRecord = this->PrepareSendScalarFieldRecord();
//
//    //By design, the current zone is the jth neighbor of zone I.
//    //How many neighbors of the current zone do you need to find out? This value is neiid.
//
//    ScalarGrid * grid = ScalarZone::GetGrid();
//    GridTopo * gridTopo = grid->gridTopo;
//    ScalarIFace * scalarIFace = gridTopo->scalarIFace;
//
//    int nNei = scalarIFace->data.size();
//    int iNei = ZoneState::inei;
//
//    ScalarIFaceIJ & sij = scalarIFace->data[ iNei ];
//    vector< int > & interfaceId = sij.ifaces;
//
//    ActionState::dataBook->MoveToBegin();
//
//    int nRecords = fieldRecord->nEquList.size();
//
//    for ( int fieldId = 0; fieldId < nRecords; ++ fieldId )
//    {
//        MRField * field  = fieldRecord->GetField( fieldId );
//        HXWriteField( ActionState::dataBook, field, interfaceId );
//    }
//
//    delete fieldRecord;
//}
//
//void FieldSolver::ReceiveSendData()
//{
//    ScalarFieldRecord * fieldRecord = this->PrepareRecvScalarFieldRecord();
//
//    //By design, the current zone is the jth neighbor of zone I.
//    //How many neighbors of the current zone do you need to find out? This value is neiid.
//
//    ScalarGrid * grid = ScalarZone::GetGrid();
//    GridTopo * gridTopo = grid->gridTopo;
//    ScalarIFace * scalarIFace = gridTopo->scalarIFace;
//
//    int nNei = scalarIFace->data.size();
//    int jNei = scalarIFace->FindINeibor( ZoneState::szid );
//
//    ScalarIFaceIJ & sij = scalarIFace->data[ jNei ];
//    vector< int > & interfaceId = sij.recv_ifaces;
//
//    ActionState::dataBook->MoveToBegin();
//
//    int nRecords = fieldRecord->nEquList.size();
//
//    for ( int fieldId = 0; fieldId < nRecords; ++ fieldId )
//    {
//        MRField * field  = fieldRecord->GetField( fieldId );
//        HXReadField( ActionState::dataBook, field, interfaceId );
//    }
//
//    delete fieldRecord;
//}

void FieldSolver::CommParallelInfo()
{
    this->UploadInterface();
    this->UpdateInterface( PrepareFieldSendData, PrepareFieldRecvData );
    this->DownloadInterface();

    this->UpdateInterface( PrepareGeomSendData, PrepareGeomRecvData );
}

void FieldSolver::SolveFlowField()
{
    for ( int n = 0; n < para->nt; ++ n )
    {
        cout << " iStep = " << n << " nStep = " << para->nt << "\n";

        this->SolveOneStep();
    }
}

void FieldSolver::SolveOneStep()
{
    this->Boundary();
    this->GetQLQR();
    this->CalcInvFlux();
    this->UpdateResidual();
    this->TimeIntergral();
    this->Update();
    this->CommParallelInfo();
    this->Visualize();
}

void FieldSolver::Boundary()
{
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        ZoneState::zid = iZone;
        this->ZoneBoundary();
    }
}

void FieldSolver::ZoneBoundary()
{
    ScalarGrid * grid = ScalarZone::GetGrid();
    int nBFaces = grid->GetNBFaces();

    RealField  & q = GetFieldReference< MRField > ( grid, "q" ).AsOneD();

    int nTCells = grid->GetNTCells();

    for ( int iFace = 0; iFace < nBFaces; ++ iFace )
    {
        int bcType = grid->bcTypes[ iFace ];
        int lc = grid->lc[ iFace ];
        int rc = grid->rc[ iFace ];
        if ( bcType == ONEFLOW::BCInflow )
        {
            Real xm = grid->xcc[ rc ];
            q[ rc ] = this->ScalarFun( xm );
        }
        else if ( bcType == ONEFLOW::BCOutflow )
        {
            q[ rc ] = q[ lc ];
        }
    }
}

void FieldSolver::GetQLQR()
{
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        ZoneState::zid = iZone;
        this->ZoneGetQLQR();
    }
}

void FieldSolver::ZoneGetQLQR()
{
    ScalarGrid * grid = ScalarZone::GetGrid();
    int nFaces = grid->GetNFaces();

    RealField & q   = GetFieldReference< MRField > ( grid, "q" ).AsOneD();
    RealField & qf1 = GetFieldReference< MRField > ( grid, "qf1" ).AsOneD();
    RealField & qf2 = GetFieldReference< MRField > ( grid, "qf2" ).AsOneD();

    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        int lc = grid->lc[ iFace ];
        int rc = grid->rc[ iFace ];

        qf1[ iFace ] = q[ lc ];
        qf2[ iFace ] = q[ rc ];
    }
}

void FieldSolver::CalcInvFlux()
{
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        ZoneState::zid = iZone;
        this->ZoneCalcInvFlux();
    }
}

void FieldSolver::ZoneCalcInvFlux()
{
    ScalarGrid * grid = ScalarZone::GetGrid();

    RealField & invflux = GetFieldReference< MRField > ( grid, "invflux" ).AsOneD();
    RealField & qf1 = GetFieldReference< MRField > ( grid, "qf1" ).AsOneD();
    RealField & qf2 = GetFieldReference< MRField > ( grid, "qf2" ).AsOneD();

    int nFaces = grid->GetNFaces();
    Real vxl = 1.0;
    Real vyl = 0.0;
    Real vzl = 0.0;

    Real vxr = 1.0;
    Real vyr = 0.0;
    Real vzr = 0.0;
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        Real q_L = qf1[ iFace ];
        Real q_R = qf2[ iFace ];

        Real vnl  = grid->xfn[ iFace ] * vxl + grid->yfn[ iFace ] * vyl + grid->zfn[ iFace ] * vzl;
        Real vnr  = grid->xfn[ iFace ] * vxr + grid->yfn[ iFace ] * vyr + grid->zfn[ iFace ] * vzr;

        Real eigenL = vnl;
        Real eigenR = vnr;

        eigenL = half * ( eigenL + ABS( eigenL ) );
        eigenR = half * ( eigenR - ABS( eigenR ) );

        Real fL = q_L * eigenL;
        Real fR = q_R * eigenR;
        Real fM = fL + fR;

        Real area = grid->area[ iFace ];
        invflux[ iFace ] = fM * area;
    }
}

void FieldSolver::UpdateResidual()
{
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        ZoneState::zid = iZone;
        this->ZoneUpdateResidual();
    }
}

void FieldSolver::ZoneUpdateResidual()
{
    ScalarGrid * grid = ScalarZone::GetGrid();

    RealField & res = GetFieldReference< MRField > ( grid, "res" ).AsOneD();
    RealField & invflux = GetFieldReference< MRField > ( grid, "invflux" ).AsOneD();

    res = 0;
    this->AddF2CField( grid, res, invflux );
}

void FieldSolver::AddF2CField( ScalarGrid * grid, RealField & cField, RealField & fField )
{
    int nFaces = grid->GetNFaces();
    int nBFaces = grid->GetNBFaces();
    for ( int iFace = 0; iFace < nBFaces; ++ iFace )
    {
        int lc = grid->lc[ iFace ];
        cField[ lc ] -= fField[ iFace ];
    }

    for ( int iFace = nBFaces; iFace < nFaces; ++ iFace )
    {
        int lc = grid->lc[ iFace ];
        int rc = grid->rc[ iFace ];

        cField[ lc ] -= fField[ iFace ];
        cField[ rc ] += fField[ iFace ];
    }
}


void FieldSolver::TimeIntergral()
{
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        ZoneState::zid = iZone;
        this->ZoneTimeIntergral();
    }
}

void FieldSolver::ZoneTimeIntergral()
{
    ScalarGrid * grid = ScalarZone::GetGrid();
    RealField & res = GetFieldReference< MRField > ( grid, "res" ).AsOneD();

    int nCells = grid->GetNCells();
    for ( int iCell = 0; iCell < nCells; ++ iCell )
    {
        Real ovol = 1.0 / grid->vol[ iCell ];
        Real coef = para->dt * ovol;
        res[ iCell ] *= coef;
    }
}

void FieldSolver::Update()
{
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        ZoneState::zid = iZone;
        this->ZoneUpdate();
    }
}

void FieldSolver::ZoneUpdate()
{
    ScalarGrid * grid = ScalarZone::GetGrid();
    RealField & q = GetFieldReference< MRField > ( grid, "q" ).AsOneD();
    RealField & res = GetFieldReference< MRField > ( grid, "res" ).AsOneD();

    int nCells = grid->GetNCells();

    for ( int iCell = 0; iCell < nCells; ++ iCell )
    {
        q[ iCell ] += res[ iCell ];
    }
}

void FieldSolver::Visualize()
{
    RealField q;
    RealField theory;
    RealField xcoor;
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        ZoneState::zid = iZone;
        this->AddVisualData( q, theory, xcoor );
    }

    this->Reorder( xcoor, q, theory );

    this->ToTecplot( xcoor, q, "OneFLOW.plt" );
    this->ToTecplot( xcoor, theory, "theory.plt" );
}


void FieldSolver::Reorder( RealField & a, RealField & b, RealField & c )
{
    int nElements = a.size();

    vector< SortArray<Real> > xlist;
    for ( int iElem = 0; iElem < nElements; ++ iElem )
    {
        SortArray<Real> ab;
        ab.data.push_back( a[ iElem ] );
        ab.data.push_back( b[ iElem ] );
        ab.data.push_back( c[ iElem ] );
        xlist.push_back( ab );
    }

    std::sort( xlist.begin(), xlist.end() );

    for ( int iElem = 0; iElem < nElements; ++ iElem )
    {
        a[ iElem ] = xlist[ iElem ].data[ 0 ];
        b[ iElem ] = xlist[ iElem ].data[ 1 ];
        c[ iElem ] = xlist[ iElem ].data[ 2 ];
    }
}

void FieldSolver::AddVisualData( RealField & qList, RealField & theoryList, RealField & xcoorList )
{
    ScalarGrid * grid = ScalarZone::GetGrid();

    RealField & q = GetFieldReference< MRField > ( grid, "q" ).AsOneD();

    int nCells = grid->GetNCells();

    Real time = para->dt * para->nt;
    Real xs = para->c * time;

    RealField theory;
    Theory( grid, time, theory );

    for ( int iCell = 0; iCell < nCells; ++ iCell )
    {
        Real xm = grid->xcc[ iCell ];
        qList.push_back( q[ iCell ] );
        theoryList.push_back( theory[ iCell ] );
        xcoorList.push_back( xm );
    }
}

void FieldSolver::Theory( ScalarGrid * grid, Real time, RealField & theory )
{
    int nCells = grid->GetNCells();
    theory.resize( nCells );

    Real xs = para->c * time;

    for ( int iCell = 0; iCell < nCells; ++ iCell )
    {
        Real xm = grid->xcc[ iCell ];
        Real xm_new = xm - xs;
        theory[ iCell ] = this->ScalarFun( xm_new );
    }
}

void FieldSolver::ToTecplot( RealField & xList, RealField & varlist, string const & fileName )
{
    fstream file;
    OpenPrjFile( file, fileName, ios_base::out );

    int nSize = xList.size();
    file << "TITLE = " << "\"OneFLOW X-Y Plot\"" << "\n";
    file << "VARIABLES = " << "\"x\", " << "\"u\" " << "\n";
    file << "ZONE T = " << "\"Scalar Results\"," << " I = " << nSize << "\n";

    for ( int iCell = 0; iCell < nSize; ++ iCell )
    {
        Real xm = xList[ iCell ];
        Real fm = varlist[ iCell ];
        file << xm << " " << fm << "\n";
    }

    CloseFile( file );
}

void PrepareFieldSendData()
{
    ScalarFieldRecord * fieldRecord = PrepareSendScalarFieldRecord();

    //By design, the current zone is the jth neighbor of zone I.
    //How many neighbors of the current zone do you need to find out? This value is neiid.

    ScalarGrid * grid = ScalarZone::GetGrid();
    GridTopo * gridTopo = grid->gridTopo;
    ScalarIFace * scalarIFace = gridTopo->scalarIFace;

    int nNei = scalarIFace->data.size();
    int iNei = ZoneState::inei;

    ScalarIFaceIJ & sij = scalarIFace->data[ iNei ];
    vector< int > & interfaceId = sij.ifaces;

    ActionState::dataBook->MoveToBegin();

    int nRecords = fieldRecord->nEquList.size();

    for ( int fieldId = 0; fieldId < nRecords; ++ fieldId )
    {
        MRField * field  = fieldRecord->GetField( fieldId );
        HXWriteField( ActionState::dataBook, field, interfaceId );
    }

    delete fieldRecord;
}

void PrepareFieldRecvData()
{
    ScalarFieldRecord * fieldRecord = PrepareRecvScalarFieldRecord();

    //By design, the current zone is the jth neighbor of zone I.
    //How many neighbors of the current zone do you need to find out? This value is neiid.

    ScalarGrid * grid = ScalarZone::GetGrid();
    GridTopo * gridTopo = grid->gridTopo;
    ScalarIFace * scalarIFace = gridTopo->scalarIFace;

    int nNei = scalarIFace->data.size();
    int jNei = scalarIFace->FindINeibor( ZoneState::szid );

    ScalarIFaceIJ & sij = scalarIFace->data[ jNei ];
    vector< int > & interfaceId = sij.recv_ifaces;

    ActionState::dataBook->MoveToBegin();

    int nRecords = fieldRecord->nEquList.size();

    for ( int fieldId = 0; fieldId < nRecords; ++ fieldId )
    {
        MRField * field  = fieldRecord->GetField( fieldId );
        HXReadField( ActionState::dataBook, field, interfaceId );
    }

    delete fieldRecord;
}

ScalarFieldRecord * PrepareSendScalarFieldRecord()
{
    ScalarFieldRecord * fieldRecord = new ScalarFieldRecord();

    ScalarGrid * grid = ScalarZone::GetGrid();
    GridTopo * gridTopo = grid->gridTopo;
    ScalarIFace * scalarIFace = gridTopo->scalarIFace;

    StringField fieldNameList;
    fieldNameList.push_back( "q" );

    fieldRecord->AddFieldRecord( scalarIFace->dataSend, fieldNameList );

    return fieldRecord;
}

ScalarFieldRecord *  PrepareRecvScalarFieldRecord()
{
    ScalarFieldRecord * fieldRecord = new ScalarFieldRecord();

    ScalarGrid * grid = ScalarZone::GetGrid();
    GridTopo * gridTopo = grid->gridTopo;
    ScalarIFace * scalarIFace = gridTopo->scalarIFace;

    StringField fieldNameList;
    fieldNameList.push_back( "q" );

    fieldRecord->AddFieldRecord( scalarIFace->dataRecv, fieldNameList );

    return fieldRecord;
}

void PrepareGeomSendData()
{
    //By design, the current zone is the jth neighbor of zone I.
    //How many neighbors of the current zone do you need to find out? This value is neiid.

    ScalarGrid * grid = ScalarZone::GetGrid();
    GridTopo * gridTopo = grid->gridTopo;
    ScalarIFace * scalarIFace = gridTopo->scalarIFace;

    int nNei = scalarIFace->data.size();
    int iNei = ZoneState::inei;

    ScalarIFaceIJ & sij = scalarIFace->data[ iNei ];
    vector< int > & interfaceId = sij.ifaces;

    ActionState::dataBook->MoveToBegin();

    for ( int iLocalFace = 0; iLocalFace < interfaceId.size(); ++ iLocalFace )
    {
        int s1;
        int interface_id = interfaceId[ iLocalFace ];

        grid->GetSId( interface_id, s1 );

        HXWrite( ActionState::dataBook, grid->xcc[ s1 ] );
        HXWrite( ActionState::dataBook, grid->ycc[ s1 ] );
        HXWrite( ActionState::dataBook, grid->zcc[ s1 ] );
        HXWrite( ActionState::dataBook, grid->vol[ s1 ] );
    }
}

void PrepareGeomRecvData()
{
    //By design, the current zone is the jth neighbor of zone I.
    //How many neighbors of the current zone do you need to find out? This value is neiid.

    ScalarGrid * grid = ScalarZone::GetGrid();
    GridTopo * gridTopo = grid->gridTopo;
    ScalarIFace * scalarIFace = gridTopo->scalarIFace;

    int nNei = scalarIFace->data.size();
    int jNei = scalarIFace->FindINeibor( ZoneState::szid );

    ScalarIFaceIJ & sij = scalarIFace->data[ jNei ];
    vector< int > & interfaceId = sij.recv_ifaces;

    ActionState::dataBook->MoveToBegin();

    for ( int iLocalFace = 0; iLocalFace < interfaceId.size(); ++ iLocalFace )
    {
        int interface_id = interfaceId[ iLocalFace ];
        int t1;
        grid->GetTId( interface_id, t1 );

        HXRead( ActionState::dataBook, grid->xcc[ t1 ] );
        HXRead( ActionState::dataBook, grid->ycc[ t1 ] );
        HXRead( ActionState::dataBook, grid->zcc[ t1 ] );
        HXRead( ActionState::dataBook, grid->vol[ t1 ] );
    }
}


EndNameSpace