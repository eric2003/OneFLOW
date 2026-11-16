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

#include "UNsVisFlux.h"
#include "NsVisFluxPolicy.h"
#include "Sutherland.h"
#include "HeatFlux.h"
#include "Zone.h"
#include "ZoneState.h"
#include "UnsGrid.h"
#include "DataBase.h"
#include "FieldWrap.h"
#include "SolverDef.h"
#include "NsCtrl.h"
#include "UCom.h"
#include "UNsCom.h"
#include "NsCom.h"
#include "VisGrad.h"
#include "UNsGrad.h"
#include "NsIdx.h"
#include "HXMath.h"
#include "Boundary.h"
#include "BcRecord.h"
#include "ULimiter.h"
#include "FieldManager.h"
#include "Iteration.h"
#include "StageProfiler.h"
#include "Ctrl.h"
#include "HipProductionPolicy.h"
#ifdef ONEFLOW_ENABLE_HIP
#include "HipFluxBackend.h"
#endif
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <fstream>



BeginNameSpace( ONEFLOW )


UNsVisFlux::UNsVisFlux()
{
    ;
}

UNsVisFlux::~UNsVisFlux()
{
    ;
}

void UNsVisFlux::SetVisPointer()
{
    //std::cout << " nscom.ivischeme = " << nscom.ivischeme << "\n";
    if ( nscom.ivischeme == VIS_STD )
    {
        this->visPointer = & UNsVisFlux::StdMethod;
    }
    else if ( nscom.ivischeme == VIS_TEST )
    {
        this->visPointer = & UNsVisFlux::TestMethod;
    }
    else if ( nscom.ivischeme == VIS_NEW1 )
    {
        this->visPointer = & UNsVisFlux::New1Method;
    }
    else if ( nscom.ivischeme == VIS_NEW2 )
    {
        this->visPointer = & UNsVisFlux::New2Method;
    }
    else if ( nscom.ivischeme == VIS_AVER )
    {
        this->visPointer = & UNsVisFlux::AverMethod;
    }
    else
    {
        this->visPointer = & UNsVisFlux::TestMethod;
    }
}

void UNsVisFlux::CalcFlux()
{
    ScopedStageTimer timer( "viscous_turbulence" );
    if ( vis_model.vismodel == 0 ) return;
    ug.Init();
    unsf.Init();
    visQ.Init( nscom.nEqu );
    visT.Init( nscom.nTModel );
    vis.Init();
    heat_flux.Init();

    visflux = std::make_unique<MRField>( nscom.nEqu, ug.nFaces );

    this->SetVisPointer();

    { ScopedStageTimer gradientTimer( "gradient" );
        this->PrepareField();
    }
    // One device face loop when the feature is enabled and the resident state is
    // bound; otherwise the unchanged host per-face loop below.  The choice is
    // made here, before either path produces output.
    if ( ! this->CalcVisFluxHip() )
    {
        this->CalcVisFlux();
    }
    this->AddVisFlux();

    visflux.reset();
}

void UNsVisFlux::PrepareField()
{
    ut_grad.Init();
    ut_grad.CalcGrad();
}

void UNsVisFlux::CalcVisFlux()
{
    // Optional whole-field dump for comparing against the device arm.  Same
    // path name as the backend's; whichever arm runs writes its own file.
    static const char * probePath = std::getenv( "ONEFLOW_VISCOUS_PROBE_FILE" );
    const bool probing = probePath != nullptr && probePath[ 0 ] != '\0';
    std::ofstream probeOut;
    if ( probing )
    {
        probeOut.open( probePath, std::ios::app );
    }

    for ( int fId = 0; fId < ug.nFaces; ++ fId )
    {
        ug.fId = fId;

        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        this->PrepareFaceValue();

        this->CalcFaceVisFlux();

        if ( probing && probeOut )
        {
            // vis.fvis is the unscaled flux; apply the same coefficient the
            // caller applies below so the two dumps are the same quantity.
            const Real coeff = ViscousFaceFluxCoefficient(
                nscom.oreynolds, gcom.farea );
            probeOut << fId;
            probeOut << ' ' << coeff * vis.fvis[ IDX::IR ];
            probeOut << ' ' << coeff * vis.fvis[ IDX::IRU ];
            probeOut << ' ' << coeff * vis.fvis[ IDX::IRV ];
            probeOut << ' ' << coeff * vis.fvis[ IDX::IRW ];
            probeOut << ' ' << coeff * vis.fvis[ IDX::IRE ];
            probeOut << '\n';
        }

        this->UpdateFaceVisFlux();
    }
}

bool UNsVisFlux::UseHipVisFlux() const
{
    return HipProductionFeatureEnabled( "ONEFLOW_ENABLE_UNS_HIP_VISFLUX" );
}

// Device face loop for the viscous flux, mirroring CalcVisFlux one-for-one.
//
// The host version above walks every face and, for each one, reads a dozen
// per-cell scalars out of the grid arrays.  On the device that becomes one
// thread per face reading the same arrays from resident storage, so the only
// per-call traffic is the two viscosity planes (Sutherland and the SA output
// are per-call inputs, not resident state) and the nFaces-wide result.
//
// Returns false when the device path is not selected, so the caller runs the
// host loop exactly as before.  It does not fall back after writing part of
// visflux: the choice is made before any output is produced.
bool UNsVisFlux::CalcVisFluxHip()
{
#ifdef ONEFLOW_ENABLE_HIP
    if ( ! this->UseHipVisFlux() ) return false;
    if ( visflux == nullptr || unsf.q == nullptr || unsf.visl == nullptr
         || unsf.vist == nullptr || unsf.tempr == nullptr
         || unsf.dtdx == nullptr || unsf.dtdy == nullptr || unsf.dtdz == nullptr )
    {
        return false;
    }
    const int nEqu = nscom.nEqu;
    if ( nEqu != 5 )
    {
        // The device kernel is the five-equation NS form.  Declining keeps a
        // three-equation case on the host path instead of failing it.
        return false;
    }
    UnsGrid * grid = Zone::GetUnsGrid();
    const int totalCells = ug.nTCell;

    // Contiguous equation-major staging for the device result.  visflux is an
    // MRField, i.e. one RealField per equation with no single data() spanning
    // them, so the device needs a flat buffer and the loop below scatters it
    // back.
    std::vector< Real > visfluxBuffer(
        static_cast< std::size_t >( nscom.nEqu ) * ug.nFaces, 0.0 );

    ViscousFluxView view;
    view.nCells = ug.nCells;
    view.nGhostCells = totalCells - ug.nCells;
    view.nFaces = ug.nFaces;
    view.nBoundaryFaces = ug.nBFaces;
    view.nEquations = nEqu;
    view.xFace = ug.xfc->data();
    view.yFace = ug.yfc->data();
    view.zFace = ug.zfc->data();
    view.xNormal = ug.xfn->data();
    view.yNormal = ug.yfn->data();
    view.zNormal = ug.zfn->data();
    view.faceArea = ug.farea->data();
    view.xCell = ug.xcc->data();
    view.yCell = ug.ycc->data();
    view.zCell = ug.zcc->data();
    view.cellVolume = ug.cvol->data();
    view.leftCell = ug.lcf->data();
    view.rightCell = ug.rcf->data();
    for ( int equation = 0; equation < nEqu; ++ equation )
    {
        view.q[ equation ] = ( * unsf.q )[ equation ].data();
        view.dqdx[ equation ] = ( * unsf.dqdx )[ equation ].data();
        view.dqdy[ equation ] = ( * unsf.dqdy )[ equation ].data();
        view.dqdz[ equation ] = ( * unsf.dqdz )[ equation ].data();
    }
    view.temperature = ( * unsf.tempr )[ 0 ].data();
    view.dTdx = ( * unsf.dtdx )[ 0 ].data();
    view.dTdy = ( * unsf.dtdy )[ 0 ].data();
    view.dTdz = ( * unsf.dtdz )[ 0 ].data();
    view.laminarViscosity = ( * unsf.visl )[ 0 ].data();
    view.turbulentViscosity = ( * unsf.vist )[ 0 ].data();
    view.gamma = nscom.gama;
    // oprl/oprt are the reciprocal Prandtl numbers (NsCom.cpp:106 sets
    // oprl = one/prl), which is what the conductivity helper multiplies by.
    // Filling these with prl/prt inverts the thermal term.
    view.inversePrandtlLaminar = nscom.oprl;
    view.inversePrandtlTurbulent = nscom.oprt;
    view.inverseReynolds = nscom.oreynolds;
    // cp comes from ConstCp(); the host reads nscom.const_cp through the
    // conductivity helper, so pass the same value rather than re-deriving it.
    view.specificHeatCp = nscom.const_cp;
    view.faceFlux = visfluxBuffer.data();
    view.faceFluxStride = nEqu;
    // grid, not nullptr: the resident gradient state was bound under this key by
    // UNsGrad::CalcGradHip, which is what lets the kernel reuse the device-side
    // primitive state and gradient instead of re-uploading them.
    view.cacheKey = grid;

    HipFluxBackend::Shared().CalcViscousFaceFlux( view );

    // MRField is Marray<Real>, a two-level vector with no contiguous data(),
    // so the device result has to be scattered back into the solver's layout.
    //
    // The kernel writes equation-major with the view's stride:
    //   faceFlux[ equation * faceFluxStride + face ]   ( faceFluxStride = nEqu )
    // This loop used to read visfluxBuffer[ equation * ug.nFaces + face ], a
    // different layout: with stride 5 the kernel only ever wrote the first
    // 4*5 + nFaces entries, so every plane read back here was almost entirely
    // zero and the viscous term silently vanished from the residual.  The
    // per-face flux itself was correct, which is why a whole-field dump showed
    // agreement while the run diverged to NaN in a hundred steps.
    for ( int equation = 0; equation < nEqu; ++ equation )
    {
        RealField & plane = ( * visflux )[ equation ];
        const Real * source =
            visfluxBuffer.data() + static_cast< std::size_t >( equation );
        for ( int face = 0; face < ug.nFaces; ++ face )
        {
            plane[ face ] = source[ face * nEqu ];
        }
    }
    return true;
#else
    return false;
#endif
}

void UNsVisFlux::CalcFaceVisFlux()
{
    this->CalcHeatFlux();

    this->CalcStress();

    this->CalcNsVisFlux();
}

void UNsVisFlux::CalcHeatFlux()
{
    this->ZeroHeatFlux();

    this->AddChemHeatFlux();

    this->AddHeatFlux();

    SaveHeatFlux();
}

void UNsVisFlux::SaveHeatFlux()
{
    if ( ug.fId >= ug.nBFaces ) return;
    if ( ug.bcRecord->bcType[ ug.fId ] != BC::SOLID_SURFACE ) return;

    // Non-owning view (.get()); unique_ptr keeps ownership.
    SurfaceValue * heat_sur = heat_flux.heatflux[ ZoneState::zid ].get();

    // Dereference unique_ptr for RealField reference (owner unchanged).
    RealField & hf = *(heat_sur->var);

    Real non_dim_heatflux = - nscom.oreynolds * vis.qNormal;
    heat_sur->var->push_back( non_dim_heatflux );
}

void UNsVisFlux::CalcStress()
{
    StokesStressTensor(
        nscom.vis,
        vis.dudx, vis.dudy, vis.dudz,
        vis.dvdx, vis.dvdy, vis.dvdz,
        vis.dwdx, vis.dwdy, vis.dwdz,
        vis.txx, vis.tyy, vis.tzz,
        vis.txy, vis.txz, vis.tyz );

    this->CalcAniStress();
}

void UNsVisFlux::CalcAniStress()
{
    if ( ctrl.nrokplus <= 0 ) return;
    Real two3rdRhok = two3rd * vis.rhok;
    vis.txx += vis.b11 - two3rdRhok;
    vis.tyy += vis.b22 - two3rdRhok;
    vis.tzz += vis.b33 - two3rdRhok;
    vis.txy += vis.b12;
    vis.txz += vis.b13;
    vis.tyz += vis.b23;
}

void UNsVisFlux::CalcNsVisFlux()
{
    vis.fvis[ IDX::IR  ] = 0.0;
    ProjectStressMomentumFlux(
        gcom.xfn, gcom.yfn, gcom.zfn,
        vis.txx, vis.tyy, vis.tzz,
        vis.txy, vis.txz, vis.tyz,
        vis.fvis[ IDX::IRU ], vis.fvis[ IDX::IRV ], vis.fvis[ IDX::IRW ] );
    vis.fvis[ IDX::IRE ] = ViscousEnergyFlux(
        vis.um, vis.vm, vis.wm,
        vis.fvis[ IDX::IRU ], vis.fvis[ IDX::IRV ], vis.fvis[ IDX::IRW ],
        vis.qNormal );
}

void UNsVisFlux::ZeroHeatFlux()
{
    vis.qNormal = 0.0;
    vis.qx      = 0.0;
    vis.qy      = 0.0;
    vis.qz      = 0.0;
}

void UNsVisFlux::AddChemHeatFlux()
{
    if ( nscom.chemModel == 1 )
    {
    }
}

void UNsVisFlux::AddHeatFlux()
{
    vis.qNormal = 0.0;
    vis.qx      = 0.0;
    vis.qy      = 0.0;
    vis.qz      = 0.0;

    nscom.kcp = EffectiveThermalConductivity(
        nscom.visl, nscom.vist, nscom.oprl, nscom.oprt, nscom.const_cp );
    vis.qNormal += FaceNormalDot( gcom.xfn, gcom.yfn, gcom.zfn, vis.qx, vis.qy, vis.qz );
    vis.qNormal += ConductiveHeatFluxNormal( nscom.kcp, vis.dtdn );
}

void UNsVisFlux::AddVisFlux()
{
    UnsGrid * grid = Zone::GetUnsGrid();
    MRField * res = FieldHome::GetResidualField( grid, NS_SOLVER );

    ONEFLOW::AddF2CField( res, visflux.get() );
    if ( ctrl.debugFieldDump != 0 )
    {
        Real mindiff = 1.0e-10;
        int idumpface = 1;
        int idumpcell = 0;
        HXDebug::DumpField( "VisFaceFlux.debug", visflux.get() );
        HXDebug::CompareFile( mindiff, idumpface );
        HXDebug::DumpResField( "VisResFlux.debug" );
        HXDebug::CompareFile( mindiff, idumpcell );
    }
}

void UNsVisFlux::PrepareFaceValue()
{
    gcom.xfn   = ( * ug.xfn   )[ ug.fId ];
    gcom.yfn   = ( * ug.yfn   )[ ug.fId ];
    gcom.zfn   = ( * ug.zfn   )[ ug.fId ];
    gcom.vfn   = ( * ug.vfn   )[ ug.fId ];
    gcom.farea = ( * ug.farea )[ ug.fId ];

    gcom.CalcTangent();

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        visQ.dqdx1[ iEqu ] = ( * unsf.dqdx )[ iEqu ][ ug.lc ];
        visQ.dqdy1[ iEqu ] = ( * unsf.dqdy )[ iEqu ][ ug.lc ];
        visQ.dqdz1[ iEqu ] = ( * unsf.dqdz )[ iEqu ][ ug.lc ];

        visQ.dqdx2[ iEqu ] = ( * unsf.dqdx )[ iEqu ][ ug.rc ];
        visQ.dqdy2[ iEqu ] = ( * unsf.dqdy )[ iEqu ][ ug.rc ];
        visQ.dqdz2[ iEqu ] = ( * unsf.dqdz )[ iEqu ][ ug.rc ];
    }

    for ( int iEqu = 0; iEqu < nscom.nTModel; ++ iEqu )
    {
        visT.dqdx1[ iEqu ] = ( * unsf.dtdx )[ iEqu ][ ug.lc ];
        visT.dqdy1[ iEqu ] = ( * unsf.dtdy )[ iEqu ][ ug.lc ];
        visT.dqdz1[ iEqu ] = ( * unsf.dtdz )[ iEqu ][ ug.lc ];

        visT.dqdx2[ iEqu ] = ( * unsf.dtdx )[ iEqu ][ ug.rc ];
        visT.dqdy2[ iEqu ] = ( * unsf.dtdy )[ iEqu ][ ug.rc ];
        visT.dqdz2[ iEqu ] = ( * unsf.dtdz )[ iEqu ][ ug.rc ];
    }

    nscom.visl1 = ( * unsf.visl )[ 0 ][ ug.lc ];
    nscom.visl2 = ( * unsf.visl )[ 0 ][ ug.rc ];

    nscom.vist1 = ( * unsf.vist )[ 0 ][ ug.lc ];
    nscom.vist2 = ( * unsf.vist )[ 0 ][ ug.rc ];

    nscom.visl = half * ( nscom.visl1 + nscom.visl2 );
    nscom.vist = half * ( nscom.vist1 + nscom.vist2 );
    nscom.vis  = nscom.visl + nscom.vist;

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        visQ.q1[ iEqu ] = ( * unsf.q )[ iEqu ][ ug.lc ];
        visQ.q2[ iEqu ] = ( * unsf.q )[ iEqu ][ ug.rc ];
    }

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        visQ.q11[ iEqu ] = visQ.q1[ iEqu ];
        visQ.q22[ iEqu ] = visQ.q2[ iEqu ];
    }

    for ( int iEqu = 0; iEqu < nscom.nTModel; ++ iEqu )
    {
        visT.q1[ iEqu ] = ( * unsf.tempr )[ iEqu ][ ug.lc ];
        visT.q2[ iEqu ] = ( * unsf.tempr )[ iEqu ][ ug.rc ];
    }

    for ( int iEqu = 0; iEqu < nscom.nTModel; ++ iEqu )
    {
        visT.q11[ iEqu ] = visT.q1[ iEqu ];
        visT.q22[ iEqu ] = visT.q2[ iEqu ];
    }

    this->AverGrad();
    this->CalcFaceWeight();

    ( this->* visPointer )();

    this->SaveFacePara();
}

void UNsVisFlux::SaveFacePara()
{
    vis.dudx  = visQ.dqdx[ IDX::IU ];
    vis.dudy  = visQ.dqdy[ IDX::IU ];
    vis.dudz  = visQ.dqdz[ IDX::IU ];

    vis.dvdx  = visQ.dqdx[ IDX::IV ];
    vis.dvdy  = visQ.dqdy[ IDX::IV ];
    vis.dvdz  = visQ.dqdz[ IDX::IV ];

    vis.dwdx  = visQ.dqdx[ IDX::IW ];
    vis.dwdy  = visQ.dqdy[ IDX::IW ];
    vis.dwdz  = visQ.dqdz[ IDX::IW ];

    vis.um  = visQ.q[ IDX::IU ];
    vis.vm  = visQ.q[ IDX::IV ];
    vis.wm  = visQ.q[ IDX::IW ];

    vis.dtdn = visT.dqdn[ IDX::ITT ];
    vis.tmid = visT.q[ IDX::ITT ];
}

void UNsVisFlux::CalcFaceWeight()
{
    vgg.CalcFaceWeight();
}

void UNsVisFlux::AverMethod()
{
    this->ZeroNormalGrad();

    this->AverFaceValue();

    this->AverGrad();

    this->CalcNormalGrad();
}

void UNsVisFlux::StdMethod()
{
    this->CalcGradCoef();

    this->ZeroNormalGrad();

    this->AverFaceValue();

    this->AverGrad();

    this->CorrectFaceGrad();

    this->CalcNormalGrad();
}

void UNsVisFlux::TestMethod()
{
    this->ZeroNormalGrad();

    this->AverFaceValue();

    this->PrepareCellGeom();

    this->CalcTestMethod();

    this->ModifyFaceGrad();
}

void UNsVisFlux::New1Method()
{
    this->ZeroNormalGrad();

    this->AccurateFaceValue();

    this->PrepareCellGeom();

    this->CalcNew1Method();

    this->ModifyFaceGrad();
}

void UNsVisFlux::New2Method()
{
    this->ZeroNormalGrad();

    this->AccurateFaceValue();

    this->PrepareCellGeom();

    this->CalcNew2Method();

    this->ModifyFaceGrad();
}

void UNsVisFlux::CalcGradCoef()
{
    vgg.CalcGradCoef();
}


void UNsVisFlux::PrepareCellGeom()
{
    vgg.PrepareCellGeom();
}

void UNsVisFlux::UpdateFaceVisFlux()
{
    const Real coeff = ViscousFaceFluxCoefficient( nscom.oreynolds, gcom.farea );
    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        ( * visflux )[ iEqu ][ ug.fId ] = coeff * vis.fvis[ iEqu ];
    }
}

void CalcLaminarViscosity( int flag )
{
    ug.Init();
    unsf.Init();
    ug.SetStEd( flag );

    Real minLimit = 0.0;

    Sutherland sutherland;
    sutherland.Init( nscom.tref_dim );

    for ( int cId = ug.ist; cId < ug.ied; ++ cId )
    {
        Real temperature = ( * unsf.tempr )[ IDX::ITT ][ cId ];
        Real visl = sutherland.CalcViscosity( temperature );
        ( * unsf.visl )[ 0 ][ cId ] = MAX( minLimit, visl );
    }
}

EndNameSpace
