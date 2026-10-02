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

#include "UNsInvFlux.h"
#include "UNsGrad.h"
#include "Zone.h"
#include "Atmosphere.h"
#include "UnsGrid.h"
#include "DataBase.h"
#include "UCom.h"
#include "UNsCom.h"
#include "NsCom.h"
#include "NsIdx.h"
#include "HXMath.h"
#include "Fatal.h"
#include "Boundary.h"
#include "BcRecord.h"
#include "UNsLimiter.h"
#include "FieldManager.h"
#include "Iteration.h"
#include "TurbCom.h"
#include "UTurbCom.h"
#include "AccelRuntime.h"
#include "EulerCpuAdapter.h"
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>


BeginNameSpace( ONEFLOW )

UNsInvFlux::UNsInvFlux()
{
    limiter = std::make_unique<NsLimiter>();
}

UNsInvFlux::~UNsInvFlux()
{
}

void UNsInvFlux::CalcLimiter()
{
    limiter->CalcLimiter();
}

void UNsInvFlux::CalcInvFace()
{
    uns_grad.Init();
    uns_grad.CalcGrad();

    this->CalcLimiter();

    this->GetQlQrField();

    this->ReconstructFaceValueField();

    this->BoundaryQlQrFixField();
}

void UNsInvFlux::GetQlQrField()
{
    this->limiter->GetQlQr();
}

void UNsInvFlux::ReconstructFaceValueField()
{
    this->limiter->CalcFaceValue();
}

void UNsInvFlux::BoundaryQlQrFixField()
{
    this->limiter->BcQlQrFix();
}

void UNsInvFlux::CalcFlux()
{
    if ( nscom.icmpInv == 0 ) return;
    inv.Init();
    ug.Init();
    unsf.Init();

    invflux = std::make_unique<MRField>( nscom.nEqu, ug.nFaces );

    this->SetPointer( nscom.ischeme );

    this->CalcInvFace();
    this->CalcInvFlux();
    this->DumpInvFluxTrace();
    this->AddInvFlux();

    invflux.reset();
}

void UNsInvFlux::CalcInvFlux()
{
    if ( this->UseCpuBatchAdapter() )
    {
        this->CalcInvFluxCpuBatch();
        return;
    }

    for ( int fId = 0; fId < ug.nFaces; ++ fId )
    {
        ug.fId = fId;

        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        this->PrepareFaceValue();

        ( this->*invFluxPointer )();

        this->UpdateFaceInvFlux();
    }
}

bool UNsInvFlux::UseCpuBatchAdapter() const
{
    const char * enabled = std::getenv( "ONEFLOW_ENABLE_UNS_CPU_BATCH" );
    if ( enabled == nullptr || enabled[ 0 ] != '1' ) return false;
    if ( AccelRuntime::Instance().IsAccelerator() ) return false;
    return nscom.ischeme == ISCHEME_LAX_FRIEDRICHS
        && nscom.nEqu == 5 && !limiter->limfIsNullPtr() && limiter->GetNEquations() == 5;
}

void UNsInvFlux::CalcInvFluxCpuBatch()
{
    const int nFaces = ug.nFaces;
    const int nEquations = limiter->GetNEquations();
    std::vector< Real > primitiveLeft( nEquations * nFaces );
    std::vector< Real > primitiveRight( nEquations * nFaces );
    std::vector< Real > xNormal( nFaces );
    std::vector< Real > yNormal( nFaces );
    std::vector< Real > zNormal( nFaces );
    std::vector< Real > meshVelocityNormal( nFaces );
    std::vector< Real > faceArea( nFaces );
    std::vector< Real > faceFlux( nEquations * nFaces );

    MRField * qf1 = limiter->GetLeftField();
    MRField * qf2 = limiter->GetRightField();

    for ( int face = 0; face < nFaces; ++ face )
    {
        xNormal[ face ] = ( * ug.xfn )[ face ];
        yNormal[ face ] = ( * ug.yfn )[ face ];
        zNormal[ face ] = ( * ug.zfn )[ face ];
        meshVelocityNormal[ face ] = ( * ug.vfn )[ face ];
        faceArea[ face ] = ( * ug.farea )[ face ];
        for ( int equation = 0; equation < nEquations; ++ equation )
        {
            primitiveLeft[ equation * nFaces + face ] =
                ( * qf1 )[ equation ][ face ];
            primitiveRight[ equation * nFaces + face ] =
                ( * qf2 )[ equation ][ face ];
        }
    }

    PrimitiveFaceStateView primitiveState;
    primitiveState.nFaces = nFaces;
    primitiveState.nEquations = nEquations;
    primitiveState.primitiveLeft = primitiveLeft.data();
    primitiveState.primitiveRight = primitiveRight.data();
    primitiveState.xNormal = xNormal.data();
    primitiveState.yNormal = yNormal.data();
    primitiveState.zNormal = zNormal.data();
    primitiveState.meshVelocityNormal = meshVelocityNormal.data();
    primitiveState.faceArea = faceArea.data();
    primitiveState.gamma = nscom.gama_ref;

    FaceFluxView flux;
    flux.nFaces = nFaces;
    flux.nEquations = nEquations;
    flux.values = faceFlux.data();

    EulerCpuAdapter adapter;
    adapter.CalcInvFlux( primitiveState, flux, 1 );
    for ( int equation = 0; equation < nEquations; ++ equation )
    {
        for ( int face = 0; face < nFaces; ++ face )
        {
            ( * invflux )[ equation ][ face ] =
                faceFlux[ equation * nFaces + face ];
        }
    }
}

void UNsInvFlux::PrepareFaceValue()
{
    gcom.xfn   = ( * ug.xfn   )[ ug.fId ];
    gcom.yfn   = ( * ug.yfn   )[ ug.fId ];
    gcom.zfn   = ( * ug.zfn   )[ ug.fId ];
    gcom.vfn   = ( * ug.vfn   )[ ug.fId ];
    gcom.farea = ( * ug.farea )[ ug.fId ];

    nscom.gama1 = ( * unsf.gama )[ 0 ][ ug.lc ];
    nscom.gama2 = ( * unsf.gama )[ 0 ][ ug.rc ];
    nscom.gama  = half * ( nscom.gama1 + nscom.gama2 );

    inv.gama1 = nscom.gama1;
    inv.gama2 = nscom.gama2;
    inv.gama  = half * ( inv.gama1 + inv.gama2 );

    MRField * qf1 = limiter->GetLeftField();
    MRField * qf2 = limiter->GetRightField();

    int nEquations = limiter->GetNEquations();

    for ( int iEqu = 0; iEqu < nEquations; ++ iEqu )
    {
        inv.prim1[ iEqu ] = ( * qf1 )[ iEqu ][ ug.fId ];
        inv.prim2[ iEqu ] = ( * qf2 )[ iEqu ][ ug.fId ];
    }
}

void UNsInvFlux::UpdateFaceInvFlux()
{
    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        ( * invflux )[ iEqu ][ ug.fId ] = gcom.farea * inv.flux[ iEqu ];
    }
}

void UNsInvFlux::DumpInvFluxTrace()
{
    const char * traceFile = std::getenv( "ONEFLOW_UNS_TRACE_FILE" );
    if ( traceFile == nullptr || traceFile[ 0 ] == '\0' ) return;

    MRField * qf1 = limiter->GetLeftField();
    MRField * qf2 = limiter->GetRightField();

    if ( limiter->limfIsNullPtr() || qf1 == nullptr || qf2 == nullptr
         || invflux == nullptr )
    {
        throw std::runtime_error(
            "UNsInvFlux trace requested before face fields are available" );
    }

    std::ofstream output( traceFile, std::ios::binary | std::ios::trunc );
    if ( ! output )
    {
        throw std::runtime_error( "cannot open UNsInvFlux trace file" );
    }

    const char magic[ 8 ] = { 'O', 'F', 'T', 'R', 'C', '0', '1', '\0' };
    const std::uint64_t nFaces = static_cast< std::uint64_t >( ug.nFaces );
    const std::uint32_t nEquations =
        static_cast< std::uint32_t >( limiter->GetNEquations() );
    const std::uint32_t nArrays = 3;
    output.write( magic, sizeof( magic ) );
    output.write(
        reinterpret_cast< const char * >( & nFaces ), sizeof( nFaces ) );
    output.write(
        reinterpret_cast< const char * >( & nEquations ),
        sizeof( nEquations ) );
    output.write(
        reinterpret_cast< const char * >( & nArrays ), sizeof( nArrays ) );

    auto writeField = [&]( const MRField & field )
    {
        for ( std::uint32_t equation = 0; equation < nEquations; ++ equation )
        {
            const auto & values = field[ equation ];
            if ( values.size() < nFaces )
            {
                throw std::runtime_error(
                    "UNsInvFlux trace field has an invalid face extent" );
            }
            output.write(
                reinterpret_cast< const char * >( values.data() ),
                static_cast< std::streamsize >(
                    nFaces * sizeof( Real ) ) );
        }
    };

    writeField( *limiter->GetLeftField());
    writeField( *limiter->GetRightField());
    writeField( *invflux );
    if ( ! output )
    {
        throw std::runtime_error( "failed while writing UNsInvFlux trace" );
    }
}

void UNsInvFlux::AddInvFlux()
{
    UnsGrid * grid = Zone::GetUnsGrid();
    MRField * res = GetFieldPointer< MRField >( grid, "res" );

    ONEFLOW::AddF2CField( res, invflux.get() );
}

EndNameSpace
