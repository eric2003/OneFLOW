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

#include "FieldSolverCuda.h"
#include "DataBase.h"
#include "Dimension.h"
#include "TimeTest.h"
#include "ScalarMetis.h"
#include "FieldPara.h"
#include "ZoneState.h"
#include "ScalarZone.h"
#include "HXMath.h"
#include "AAASolverCuda.h"
#ifdef ENABLE_CUDA
#include <cuda_runtime.h>
#endif

BeginNameSpace( ONEFLOW )

FieldSolverCuda::FieldSolverCuda()
{
}

FieldSolverCuda::~FieldSolverCuda()
{
}

void FieldSolverCuda::Run()
{
    int scalar_flag = ONEFLOW::GetDataValue< int >("scalar_flag");
    Dim::SetDimension( ONEFLOW::GetDataValue< int >( "dimension" ) );

    TimeTest ts;
    ts.RunTest();

    if ( scalar_flag == 0 )
    {
        ScalarMetis::Create1DMesh();
    }
    else if ( scalar_flag == 1 )
    {
        ScalarMetis::Create1DMeshFromCgns();
    }
    else if ( scalar_flag == 2 )
    {
        ScalarMetis::Run();
    }
    else 
    {
        this->Init();

        this->SolveFlowField();
    }
}

void FieldSolverCuda::SolveFlowField()
{
    TimeTest ts;
    for ( int n = 0; n < para->nt; ++ n )
    {
        if ( ( ( n + 1 ) % 200 ) == 0 )
        {
            std::cout << " iStep = " << n + 1 << " nStep = " << para->nt << "\n";
        }

        this->SolveOneStep();
        if ( ( ( n + 1 ) % 200 ) == 0 )
        {
            ts.ShowTimeSpan();
        }


    }
    this->Visualize();
}

void FieldSolverCuda::SolveOneStep()
{
    this->Boundary();
    this->GetQLQR();
    this->CalcInvFlux();
    this->UpdateResidual();
    this->TimeIntergral();
    this->Update();
    this->CommParallelInfo();
}

void FieldSolverCuda::Boundary()
{
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        if ( ! ZoneState::IsValidZone( iZone ) ) continue;
        ZoneState::zid = iZone;

        this->ZoneBoundary();
    }
}

void FieldSolverCuda::ZoneBoundary()
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

void FieldSolverCuda::GetQLQR()
{
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        if ( ! ZoneState::IsValidZone( iZone ) ) continue;
        ZoneState::zid = iZone;
        this->ZoneGetQLQR();
    }
}

void TestAdd()
{
    const int arraySize = 5;
    int a[arraySize] = { 1, 2, 3, 4, 5 };
    int b[arraySize] = { 10, 20, 30, 40, 50 };
    int c[arraySize] = { 0 };

    // Add vectors in parallel.
    addWithCuda( a, b, c, arraySize);

    printf("{1,2,3,4,5} + {10,20,30,40,50} = {%d,%d,%d,%d,%d}\n",
        c[0], c[1], c[2], c[3], c[4]);

    cudaDeviceReset();
}

void TestAddReal()
{
    const int arraySize = 5;
    Real a[arraySize] = { 1.1, 2.2, 3.3, 4.4, 5.5 };
    Real b[arraySize] = { 10, 20, 30, 40, 50 };
    Real c[arraySize] = { 0 };

    // Add vectors in parallel.
    addRealWithCuda( a, b, c, arraySize);

    printf("{1.1, 2.2, 3.3, 4.4, 5.5} + {10,20,30,40,50} = {%f,%f,%f,%f,%f}\n",
        c[0], c[1], c[2], c[3], c[4]);

    cudaDeviceReset();
}

void TestAddRealSwap()
{
    const int arraySize = 5;
    Real a[arraySize] = { 1.1, 2.2, 3.3, 4.4, 5.5 };
    Real b[arraySize] = { 10, 20, 30, 40, 50 };
    int  id[arraySize] = { 4, 3, 2, 1, 0 };
    Real c[arraySize] = { 0 };

    // Add vectors in parallel.
    addRealSwapWithCuda( a, b, id, c, arraySize);
    //51.1,42.2,33.3,24.4,15.5
    printf("{1.1, 2.2, 3.3, 4.4, 5.5} + {10,20,30,40,50} = {%f,%f,%f,%f,%f}\n",
        c[0], c[1], c[2], c[3], c[4]);

    cudaDeviceReset();
}

void TestVectorAddRealSwap()
{
    const int arraySize = 5;
    std::vector<Real> a = { 1.1, 2.2, 3.3, 4.4, 5.5 };
    std::vector<Real> b = { 10, 20, 30, 40, 50 };
    std::vector<int>  id = { 4, 3, 2, 1, 0 };
    std::vector<Real> c = { 0,0,0,0,0 };

    // Add vectors in parallel.
    addRealSwapWithCuda( &a[0], &b[0], &id[0], &c[0], arraySize);
    //51.1,42.2,33.3,24.4,15.5
    printf("TestVectorAddRealSwap {1.1, 2.2, 3.3, 4.4, 5.5} + {10,20,30,40,50} = {%f,%f,%f,%f,%f}\n",
        c[0], c[1], c[2], c[3], c[4]);

    cudaDeviceReset();
}

void TestVectorSetRealSwap()
{
    const int arraySize = 5;
    std::vector<Real> a = { 1.1, 2.2, 3.3, 4.4, 5.5 };
    std::vector<int>  id = { 4, 3, 2, 1, 0 };
    std::vector<Real> c = { 0,0,0,0,0 };

    // Add vectors in parallel.
    setRealSwapWithCuda( &a[0], &id[0], &c[0], arraySize);
    printf("TestVectorSetRealSwap  = {%f,%f,%f,%f,%f}\n",
        c[0], c[1], c[2], c[3], c[4]);

    cudaDeviceReset();
}

void TestVectorSetRealSwapNew()
{
    const int arraySize = 5;
    std::vector<int>  id = { 4, 3, 2, 1, 0 };
    std::vector<Real> a = { 0,0,0,0,0 };
    std::vector<Real> b = { 1.1, 2.2, 3.3, 4.4, 5.5 };

    // Add vectors in parallel.
    setRealSwapWithCudaNew( &a[0], &b[0], &id[0], arraySize);
    printf("TestVectorSetRealSwapNew a  = {%f,%f,%f,%f,%f}\n",
        a[0], a[1], a[2], a[3], a[4]);

    cudaDeviceReset();
}

void SetRealSwapWithCudaNewRealProblemCpu(Real *a, Real *b, int * id, unsigned int nFaces, unsigned int nCells)
{
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        int j = id[ iFace ];

        a[ iFace ] = b[ j ];
    }
}

void TestVectorSetRealSwapNewRealProblem()
{
    ScalarGrid * grid = ScalarZone::GetGrid();
    int nFaces = grid->GetNFaces();
    int nCells = grid->GetNCells();

    //RealField & q   = GetFieldReference< MRField > ( grid, "q" ).AsOneD();
    //RealField & qf1 = GetFieldReference< MRField > ( grid, "qf1" ).AsOneD();
    //RealField & qf2 = GetFieldReference< MRField > ( grid, "qf2" ).AsOneD();

    RealField q   = GetFieldReference< MRField > ( grid, "q" ).AsOneD();
    RealField qf1 = GetFieldReference< MRField > ( grid, "qf1" ).AsOneD();
    RealField qf2 = GetFieldReference< MRField > ( grid, "qf2" ).AsOneD();


    //std::vector<Real> qq =

    setRealSwapWithCudaNewRealProblem( &qf1[0], &q[0], &grid->lc.data[0], nFaces, nCells );
    printf("TestVectorSetRealSwapNewRealProblem a  = {%f,%f,%f,%f,%f}\n",
        qf1[0], qf1[1], qf1[2], qf1[3], qf1[4]);

    cudaDeviceReset();

    SetRealSwapWithCudaNewRealProblemCpu( &qf1[0], &q[0], &grid->lc.data[0], nFaces, nCells );
    printf("TestVectorSetRealSwapNewRealProblem a111  = {%f,%f,%f,%f,%f}\n",
        qf1[0], qf1[1], qf1[2], qf1[3], qf1[4]);
    cudaDeviceReset();

    std::vector<int> id(nFaces);
    for( int i = 0; i < nFaces; ++ i )
    {
        id[i] = i;
    }
    std::vector<Real> xx(nFaces, 0);
    std::vector<Real> yy(nFaces, 1);
    setRealSwapWithCudaNewRealProblem( &xx[0], &yy[0], &id[0], nFaces, nFaces );
    printf("TestVectorSetRealSwapNewRealProblem xx  = {%f,%f,%f,%f,%f}\n",
        xx[0], xx[1], xx[2], xx[3], xx[4]);
    cudaDeviceReset();

    setRealSwapWithCudaNew( &xx[0], &yy[0], &id[0], nFaces);
    printf("TestVectorSetRealSwapNew xx111  = {%f,%f,%f,%f,%f}\n",
        xx[0], xx[1], xx[2], xx[3], xx[4]);
    cudaDeviceReset();
    TestVectorSetRealSwapNew();
    cudaDeviceReset();
    //int nnn = 6;
    //int nnn = 10; //ok
    //int nnn = 2000; //error
    //int nnn = 1000; 
    int nnn = 5000; 
    std::vector<int>  iddd = { 4, 3, 2, 1, 0 };
    std::vector<Real> a(nnn, 0);
    std::vector<Real> b = { 1.1, 2.2, 3.3, 4.4, 5.5 };
    b.resize( nnn );
    iddd.resize( nnn );

    // Add vectors in parallel.
    setRealSwapWithCudaNew( &a[0], &b[0], &iddd[0], nnn);
    printf("TestVectorSetRealSwapNew a1  = {%f,%f,%f,%f,%f}\n",
        a[0], a[1], a[2], a[3], a[4]);
    int kkk = 1;
}


void FieldSolverCuda::ZoneGetQLQR()
{
    ScalarGrid * grid = ScalarZone::GetGrid();
    int nFaces = grid->GetNFaces();

    RealField & q   = GetFieldReference< MRField > ( grid, "q" ).AsOneD();
    RealField & qf1 = GetFieldReference< MRField > ( grid, "qf1" ).AsOneD();
    RealField & qf2 = GetFieldReference< MRField > ( grid, "qf2" ).AsOneD();

#ifdef ENABLE_CUDA
    //TestAdd();
    //TestAddReal();
    //TestAddRealSwap();
    //TestVectorAddRealSwap();
    //TestVectorSetRealSwap();
    //TestVectorSetRealSwapNew();
    //TestVectorSetRealSwapNewRealProblem();
    int nBFaces = grid->GetNBFaces();
    int nCells = grid->GetNCells();
    int nTCells = nCells + nBFaces;
    //int nTCells = nCells;
    SetValueWithCuda(&qf1[0], &q[0], &grid->lc.data[0], nFaces, nTCells);
    SetValueWithCuda(&qf2[0], &q[0], &grid->rc.data[0], nFaces, nTCells);
    int kkk = 1;
    //for ( int iFace = 0; iFace < nFaces; ++ iFace )
    //{
    //    int lc = grid->lc[ iFace ];
    //    int rc = grid->rc[ iFace ];

    //    qf1[ iFace ] = q[ lc ];
    //    qf2[ iFace ] = q[ rc ];
    //}
#else
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        int lc = grid->lc[ iFace ];
        int rc = grid->rc[ iFace ];

        qf1[ iFace ] = q[ lc ];
        qf2[ iFace ] = q[ rc ];
    }
#endif
}

//void FieldSolverCuda::ZoneGetQLQR()
//{
//    ScalarGrid * grid = ScalarZone::GetGrid();
//    int nFaces = grid->GetNFaces();
//
//    RealField & q   = GetFieldReference< MRField > ( grid, "q" ).AsOneD();
//    RealField & qf1 = GetFieldReference< MRField > ( grid, "qf1" ).AsOneD();
//    RealField & qf2 = GetFieldReference< MRField > ( grid, "qf2" ).AsOneD();
//
//    //int nCells = grid->GetNCells();
//    //SetValueWithCuda(&qf1[0], &q[0], &grid->lc[0], nFaces, nCells);
//    //SetValueWithCuda(&qf2[0], &q[0], &grid->rc[0], nFaces, nCells);
//    for ( int iFace = 0; iFace < nFaces; ++ iFace )
//    {
//        int lc = grid->lc[ iFace ];
//        int rc = grid->rc[ iFace ];
//
//        qf1[ iFace ] = q[ lc ];
//        qf2[ iFace ] = q[ rc ];
//    }
//}

void FieldSolverCuda::CalcInvFlux()
{
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        if ( ! ZoneState::IsValidZone( iZone ) ) continue;
        ZoneState::zid = iZone;
        this->ZoneCalcInvFlux();
    }
}

void FieldSolverCuda::ZoneCalcInvFlux()
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

void FieldSolverCuda::ZoneCalcInvFluxCuda()
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

    MyCalcInvFluxCuda(&qf1[0], &qf2[0], &invflux[0], &grid->xfn.data[0], &grid->yfn.data[0], &grid->zfn.data[0], &grid->area.data[0], nFaces);
}

void FieldSolverCuda::UpdateResidual()
{
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        if ( ! ZoneState::IsValidZone( iZone ) ) continue;
        ZoneState::zid = iZone;
        this->ZoneUpdateResidual();
    }
}

void FieldSolverCuda::ZoneUpdateResidual()
{
    ScalarGrid * grid = ScalarZone::GetGrid();

    RealField & res = GetFieldReference< MRField > ( grid, "res" ).AsOneD();
    RealField & invflux = GetFieldReference< MRField > ( grid, "invflux" ).AsOneD();

    res = 0;
    //this->AddF2CField( grid, res, invflux );
    //RealField tmp_res = res;
    //this->AddF2CField( grid, tmp_res, invflux );
    this->AddF2CFieldCuda( grid, res, invflux );
    int kkk = 1;
}

void FieldSolverCuda::AddF2CField( ScalarGrid * grid, RealField & cField, RealField & fField )
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

void FieldSolverCuda::AddF2CFieldCuda( ScalarGrid * grid, RealField & cField, RealField & fField )
{
    int nFaces = grid->GetNFaces();
    int nBFaces = grid->GetNBFaces();
    int nCells = grid->GetNCells();
    int nTCells = nCells + nBFaces;

    MyAddF2CFieldCuda(&fField[0], &cField[0], &grid->lc.data[0], &grid->rc.data[0], nBFaces, nFaces, nTCells);
}

void FieldSolverCuda::TimeIntergral()
{
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        if ( ! ZoneState::IsValidZone( iZone ) ) continue;
        ZoneState::zid = iZone;
        //this->ZoneTimeIntergral();
        this->ZoneTimeIntergralCuda();
    }
}

void FieldSolverCuda::ZoneTimeIntergral()
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

void FieldSolverCuda::ZoneTimeIntergralCuda()
{
    ScalarGrid * grid = ScalarZone::GetGrid();
    RealField & res = GetFieldReference< MRField > ( grid, "res" ).AsOneD();

    int nCells = grid->GetNCells();
    int nBFaces = grid->GetNBFaces();
    int nTCells = nCells + nBFaces;
    //MyZoneTimeIntergralCuda(&res[0], &grid->vol.data[0], para->dt, nCells, nTCells);
    MyZoneTimeIntergralCuda(&res[0], &grid->vol.data[0], para->dt, nCells);
    int kkk = 1;
}

void FieldSolverCuda::Update()
{
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        if ( ! ZoneState::IsValidZone( iZone ) ) continue;
        ZoneState::zid = iZone;
        //this->ZoneUpdate();
        this->ZoneUpdateCuda();
    }
}

void FieldSolverCuda::ZoneUpdate()
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

void FieldSolverCuda::ZoneUpdateCuda()
{
    ScalarGrid * grid = ScalarZone::GetGrid();
    RealField & q = GetFieldReference< MRField > ( grid, "q" ).AsOneD();
    RealField & res = GetFieldReference< MRField > ( grid, "res" ).AsOneD();

    int nCells = grid->GetNCells();
    MyZoneUpdateCuda(&q[0], &res[0], nCells);
}


EndNameSpace
