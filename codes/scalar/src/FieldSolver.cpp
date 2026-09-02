/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
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
#include "DataBase.h"
#include "Dimension.h"
#include "TimeTest.h"
#include "ScalarMetis.h"
#include "FieldPara.h"
#include "ZoneState.h"
#include "ScalarZone.h"
#include "HXMath.h"
#include "AccelRuntime.h"
#include "CpuFluxBackend.h"
#ifdef ONEFLOW_ENABLE_HIP
#include "HipFluxBackend.h"
#endif
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>


BeginNameSpace( ONEFLOW )

namespace
{

FluxBackend & GetScalarFluxBackend()
{
    static CpuFluxBackend cpuBackend;
#ifdef ONEFLOW_ENABLE_HIP
    static HipFluxBackend hipBackend;
    if ( AccelRuntime::Instance().IsAccelerator() ) return hipBackend;
#endif
    return cpuBackend;
}

bool IsScalarAccelValidationEnabled()
{
    const char * value = std::getenv( "ONEFLOW_ACCEL_VALIDATE" );
    return value != nullptr && std::string( value ) == "1";
}

void CheckScalarAccelError(
    const char * operation,
    const std::vector< Real > & reference,
    const Real * accelerated )
{
    Real maxError = 0.0;
    for ( std::size_t i = 0; i < reference.size(); ++ i )
    {
        maxError = std::max(
            maxError,
            std::abs( reference[ i ] - accelerated[ i ] ) );
    }
    if ( maxError > 1.0e-15 )
    {
        throw std::runtime_error(
            std::string( "OneFLOW scalar accelerator validation failed in " )
            + operation + ": max absolute error = "
            + std::to_string( maxError ) );
    }
}

}

FieldSolver::FieldSolver()
{
}

FieldSolver::~FieldSolver()
{
}

void FieldSolver::Run()
{
    TestMPI();

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

void FieldSolver::SolveFlowField()
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
            //this->Visualize();
        }


    }
    this->Visualize();
}

//void FieldSolver::SolveOneStep()
//{
//    TimeTest ts;
//    this->Boundary();
//    ts.ShowTimeSpan("Boundary");
//    this->GetQLQR();
//    ts.ShowTimeSpan("GetQLQR");
//    this->CalcInvFlux();
//    ts.ShowTimeSpan("CalcInvFlux");
//    this->UpdateResidual();
//    ts.ShowTimeSpan("UpdateResidual");
//    this->TimeIntergral();
//    ts.ShowTimeSpan("TimeIntergral");
//    this->Update();
//    ts.ShowTimeSpan("Update");
//    this->CommParallelInfo();
//    ts.ShowTimeSpan("CommParallelInfo");
//    //this->Visualize();
//    //ts.ShowTimeSpan("Visualize");
//}

void FieldSolver::SolveOneStep()
{
    this->Boundary();
    this->GetQLQR();
    this->CalcInvFlux();
    this->UpdateResidual();
    this->TimeIntergral();
    this->Update();
    this->CommParallelInfo();
}

void FieldSolver::Boundary()
{
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        if ( ! ZoneState::IsValidZone( iZone ) ) continue;
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
        if ( ! ZoneState::IsValidZone( iZone ) ) continue;
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
        if ( ! ZoneState::IsValidZone( iZone ) ) continue;
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

    const int nFaces = grid->GetNFaces();
    FaceStateView state;
    state.nFaces = nFaces;
    state.nEquations = 1;
    state.qLeft = &qf1[ 0 ];
    state.qRight = &qf2[ 0 ];
    state.xNormal = &grid->xfn[ 0 ];
    state.yNormal = &grid->yfn[ 0 ];
    state.zNormal = &grid->zfn[ 0 ];
    state.faceArea = &grid->area[ 0 ];

    FaceFluxView flux;
    flux.nFaces = nFaces;
    flux.nEquations = 1;
    flux.values = &invflux[ 0 ];
    GetScalarFluxBackend().CalcInvFlux( state, flux, 0 );

    if ( AccelRuntime::Instance().IsAccelerator()
         && IsScalarAccelValidationEnabled() )
    {
        std::vector< Real > referenceValues( nFaces, 0.0 );
        FaceFluxView referenceFlux;
        referenceFlux.nFaces = nFaces;
        referenceFlux.nEquations = 1;
        referenceFlux.values = referenceValues.data();
        CpuFluxBackend cpuBackend;
        cpuBackend.CalcInvFlux( state, referenceFlux, 0 );
        CheckScalarAccelError( "CalcInvFlux", referenceValues, &invflux[ 0 ] );
    }
}

void FieldSolver::UpdateResidual()
{
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        if ( ! ZoneState::IsValidZone( iZone ) ) continue;
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
    const int nFaces = grid->GetNFaces();
    FaceFluxView flux;
    flux.nFaces = nFaces;
    flux.nEquations = 1;
    flux.values = &invflux[ 0 ];

    FaceConnectivityView connectivity;
    connectivity.nFaces = nFaces;
    connectivity.nBoundaryFaces = grid->GetNBFaces();
    connectivity.leftCell = &grid->lc[ 0 ];
    connectivity.rightCell = &grid->rc[ 0 ];

    ResidualView residual;
    residual.nCells = grid->GetNCells();
    residual.nEquations = 1;
    residual.values = &res[ 0 ];
    GetScalarFluxBackend().AddFaceFlux( flux, connectivity, residual );

    if ( AccelRuntime::Instance().IsAccelerator()
         && IsScalarAccelValidationEnabled() )
    {
        std::vector< Real > referenceValues( residual.nCells, 0.0 );
        ResidualView referenceResidual;
        referenceResidual.nCells = residual.nCells;
        referenceResidual.nEquations = 1;
        referenceResidual.values = referenceValues.data();
        CpuFluxBackend cpuBackend;
        cpuBackend.AddFaceFlux( flux, connectivity, referenceResidual );
        CheckScalarAccelError(
            "AddFaceFlux", referenceValues, residual.values );
    }
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
        if ( ! ZoneState::IsValidZone( iZone ) ) continue;
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
        if ( ! ZoneState::IsValidZone( iZone ) ) continue;
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

EndNameSpace
