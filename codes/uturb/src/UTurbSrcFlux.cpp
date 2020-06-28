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
#include "UTurbSrcFlux.h"
#include "UTurbGrad.h"
#include "UNsGrad.h"
#include "UTurbCom.h"
#include "NsCtrl.h"
#include "NsIdx.h"
#include "TurbCom.h"
#include "UTurbCom.h"
#include "Com.h"
#include "UCom.h"
#include "DataBase.h"
#include "UnsGrid.h"
#include "Zone.h"
#include "ZoneState.h"
#include "HXMath.h"
#include "CellMesh.h"
#include "BcRecord.h"
#include "Boundary.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

UTurbSrcFlux::UTurbSrcFlux()
{
    ;
}

UTurbSrcFlux::~UTurbSrcFlux()
{
    ;
}

void UTurbSrcFlux::Init()
{
    ug.Init();
    turbcom.Init();
    this->CalcGrad();
    if ( turbcom.nEqu == 1 )
    {
        this->CalcLengthScaleSa();
    }
    else if ( turbcom.nEqu >= 2 )
    {
        if ( vis_model.visname.substr( 0, 13 ) == "2eq-kw-menter"  )
        {
            this->CalcCrossingTerm();
            this->CalcBlendingTerm();
            this->CalcLengthScaleSst();
        }
        else if ( vis_model.visname.substr( 0, 12 ) == "easm-kw-2005" )
        {
            this->CalcCrossingTerm();
            this->CalcBlendingTerm();
        }
        else
        {
            this->CalcCrossingTerm();
        }
    }
    this->SetSrcFluxPointer();

}

void UTurbSrcFlux::CalcSrcFlux()
{
    //ReadTmp();
    Init();
    if ( turbcom.nEqu == 1 )
    {
        this->CalcSrcFlux1Equ();
    }
    else if ( turbcom.nEqu >= 2 )
    {
        this->CalcSrcFlux2Equ();
    }
}

void UTurbSrcFlux::CalcVist()
{
    InitVist();
    if ( turbcom.nEqu == 1 )
    {
        this->CalcVist1Equ();
    }
    else if ( turbcom.nEqu >= 2 )
    {
        this->CalcVist2Equ();
    }
    SetGhostCellVist();
    CalcVistMax();
}

void UTurbSrcFlux::SetGhostCellVist()
{
    for ( int fId = 0; fId < ug.nBFace; ++ fId )
    {
        int lc  = ( * ug.lcf ) [ fId ];
        int rc  = ( * ug.rcf ) [ fId ];

        ( * uturbf.vist )[ 0 ][ rc ] = ( * uturbf.vist )[ 0 ][ lc ];
    }
}

void UTurbSrcFlux::InitVist()
{
    ug.Init();
    turbcom.Init();
    uturbf.Init();
    this->CalcGrad();
}

void UTurbSrcFlux::ZeroSpectrum()
{
    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        ug.cId = cId;
        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            ( * uturbf.impsr )[ iEqu ][ ug.cId ] = 0.0;
        }
    }
}

void UTurbSrcFlux::ReadTmp()
{
    static int iii = 0;
    if ( iii ) return;
    iii = 1;
    fstream file;
    file.open( "turbflowsrc.dat", ios_base::in | ios_base::binary );
    if ( ! file ) exit( 0 );

       for ( int cId = 0; cId < ug.nTCell; ++ cId )
    {
        for ( int iEqu = 0; iEqu < 5; ++ iEqu )
        {
            file.read( reinterpret_cast< char * >( & ( * uturbf.q_ns )[ iEqu ][ cId ] ), sizeof( double ) );
        }
    }

       for ( int cId = 0; cId < ug.nTCell; ++ cId )
    {
        file.read( reinterpret_cast< char * >( & ( * uturbf.dqdx_ns )[ IDX::IU ][ cId ] ), sizeof( double ) );
        file.read( reinterpret_cast< char * >( & ( * uturbf.dqdy_ns )[ IDX::IU ][ cId ] ), sizeof( double ) );
        file.read( reinterpret_cast< char * >( & ( * uturbf.dqdz_ns )[ IDX::IU ][ cId ] ), sizeof( double ) );
        file.read( reinterpret_cast< char * >( & ( * uturbf.dqdx_ns )[ IDX::IV ][ cId ] ), sizeof( double ) );
        file.read( reinterpret_cast< char * >( & ( * uturbf.dqdy_ns )[ IDX::IV ][ cId ] ), sizeof( double ) );
        file.read( reinterpret_cast< char * >( & ( * uturbf.dqdz_ns )[ IDX::IV ][ cId ] ), sizeof( double ) );
        file.read( reinterpret_cast< char * >( & ( * uturbf.dqdx_ns )[ IDX::IW ][ cId ] ), sizeof( double ) );
        file.read( reinterpret_cast< char * >( & ( * uturbf.dqdy_ns )[ IDX::IW ][ cId ] ), sizeof( double ) );
        file.read( reinterpret_cast< char * >( & ( * uturbf.dqdz_ns )[ IDX::IW ][ cId ] ), sizeof( double ) );
    }

       for ( int cId = 0; cId < ug.nTCell; ++ cId )
    {
        file.read( reinterpret_cast< char * >( & ( * uturbf.visl )[ 0 ][ cId ] ), sizeof( double ) );
    }

       for ( int cId = 0; cId < ug.nTCell; ++ cId )
    {
        file.read( reinterpret_cast< char * >( & ( * uturbf.vist )[ 0 ][ cId ] ), sizeof( double ) );
    }

    for ( int iCell = 0; iCell < ug.nTCell; ++ iCell )
    {
        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            file.read( reinterpret_cast< char * >( & ( * uturbf.q )[ iEqu ][ iCell ] ), sizeof( double ) );
        }
    }
    file.close();
    file.clear();
}

void UTurbSrcFlux::CalcVistMax()
{
    turbcom.maxid = 0;
    turbcom.maxvist = -1.0;
    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        ug.cId = cId;

        turbcom.vist = ( * uturbf.vist )[ 0 ][ ug.cId ];

        if ( turbcom.maxvist < turbcom.vist )
        {
            turbcom.maxvist = turbcom.vist;
            turbcom.maxid = cId;
        }
    }
    //cout << "maxvist = " << turbcom.maxvist << " cell id = " << turbcom.maxid << endl;
}

void UTurbSrcFlux::CalcVist1Equ()
{
    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        ug.cId = cId;

        if ( cId == 6155 )
        {
            int kkk = 1;
        }
        if ( turbcom.rho < 0 || NotANumber( turbcom.rho ) )
        {
            cout << " zone = " << ZoneState::zid << " cId = " << cId << " rho = " << turbcom.rho << "\n";
            cin.get();
        }
        turbcom.rho  = ( * uturbf.q_ns  )[ IDX::IR  ][ ug.cId ];
        turbcom.nuet = ( * uturbf.q  )[ ISA ][ ug.cId ];
        turbcom.visl = ( * uturbf.visl )[ 0 ][ ug.cId ];

        this->CalcCellVist1Equ();

        ( * uturbf.vist )[ 0 ][ ug.cId ] = turbcom.vist;
    }
}

void UTurbSrcFlux::CalcVist2Equ()
{
    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        ug.cId = cId;

        turbcom.dudx = ( * uturbf.dqdx_ns )[ IDX::IU ][ ug.cId ];
        turbcom.dudy = ( * uturbf.dqdy_ns )[ IDX::IU ][ ug.cId ];
        turbcom.dudz = ( * uturbf.dqdz_ns )[ IDX::IU ][ ug.cId ];
        turbcom.dvdx = ( * uturbf.dqdx_ns )[ IDX::IV ][ ug.cId ];
        turbcom.dvdy = ( * uturbf.dqdy_ns )[ IDX::IV ][ ug.cId ];
        turbcom.dvdz = ( * uturbf.dqdz_ns )[ IDX::IV ][ ug.cId ];
        turbcom.dwdx = ( * uturbf.dqdx_ns )[ IDX::IW ][ ug.cId ];
        turbcom.dwdy = ( * uturbf.dqdy_ns )[ IDX::IW ][ ug.cId ];
        turbcom.dwdz = ( * uturbf.dqdz_ns )[ IDX::IW ][ ug.cId ];

        turbcom.rho  = ( * uturbf.q_ns  )[ IDX::IR  ][ ug.cId ];
        turbcom.ke   = ( * uturbf.q  )[ IKE ][ ug.cId ];
        turbcom.kw   = ( * uturbf.q  )[ IKW ][ ug.cId ];
        turbcom.visl = ( * uturbf.visl )[ 0 ][ ug.cId ];
        turbcom.vist = ( * uturbf.vist )[ 0 ][ ug.cId ];
        turbcom.dist  = ( * uturbf.dist )[ ug.cId ];

        this->CalcCellVist2Equ();

        ( * uturbf.vist )[ 0 ][ ug.cId ] = turbcom.vist;
    }
}

void UTurbSrcFlux::CalcSrcFlux1Equ()
{
    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        ug.cId = cId;

        if ( cId == 22 )
        {
            int kkk = 1;
        }

        this->PrepareCellValue1Equ();

        this->CalcSrcSa();

        this->Update1Equ();
    }
}

void UTurbSrcFlux::CalcSrcFlux2Equ()
{
    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        ug.cId = cId;

        if ( cId == 11 )
        {
            int kkk = 1;
        }

        this->PrepareCellValue();

        this->CalcSrc2Equ();

        this->Update2Equ();
    }
}

void UTurbSrcFlux::Update1Equ()
{
    ( * uturbf.impsr )[ ISA ][ ug.cId ] += turbcom.spec_sa * gcom.cvol;
    ( * uturbf.res   )[ ISA ][ ug.cId ] += turbcom.res_sa * gcom.cvol;
}

void UTurbSrcFlux::Update2Equ()
{
    ( * uturbf.res )[ IKE ][ ug.cId ] += ( turbcom.srck ) * gcom.cvol;
    ( * uturbf.res )[ IKW ][ ug.cId ] += ( turbcom.srcw ) * gcom.cvol;

    ( * uturbf.impsr )[ IKE ][ ug.cId ] += - ( turbcom.diak + turbcom.fskn * turbcom.oork ) * gcom.cvol;
    ( * uturbf.impsr )[ IKW ][ ug.cId ] += - ( turbcom.diaw + turbcom.fswn * turbcom.oorw ) * gcom.cvol;

    if ( turbcom.transition_model == ITReGama ) 
    {
        ( * uturbf.res )[ ITGama ][ ug.cId ] += turbcom.srcg * gcom.cvol;
        ( * uturbf.res )[ ITRect ][ ug.cId ] += turbcom.srcr * gcom.cvol;

        ( * uturbf.impsr )[ ITGama ][ ug.cId ] += turbcom.speg * gcom.cvol;
        ( * uturbf.impsr )[ ITRect ][ ug.cId ] += turbcom.sper * gcom.cvol;
    }
}

void UTurbSrcFlux::CalcGradDebug()
{
    uturb_grad.Init();
    uturb_grad.CalcGrad();
    uns_grad.Init();
    uns_grad.CalcGrad();
}

void UTurbSrcFlux::CalcGrad()
{
    uturb_grad.Init();
    uturb_grad.CalcGrad();
    //uturb_grad.CalcGradDebug();
}

void UTurbSrcFlux::PrepareCellValue()
{
    turbcom.dudx = ( * uturbf.dqdx_ns )[ IDX::IU ][ ug.cId ];
    turbcom.dudy = ( * uturbf.dqdy_ns )[ IDX::IU ][ ug.cId ];
    turbcom.dudz = ( * uturbf.dqdz_ns )[ IDX::IU ][ ug.cId ];
    turbcom.dvdx = ( * uturbf.dqdx_ns )[ IDX::IV ][ ug.cId ];
    turbcom.dvdy = ( * uturbf.dqdy_ns )[ IDX::IV ][ ug.cId ];
    turbcom.dvdz = ( * uturbf.dqdz_ns )[ IDX::IV ][ ug.cId ];
    turbcom.dwdx = ( * uturbf.dqdx_ns )[ IDX::IW ][ ug.cId ];
    turbcom.dwdy = ( * uturbf.dqdy_ns )[ IDX::IW ][ ug.cId ];
    turbcom.dwdz = ( * uturbf.dqdz_ns )[ IDX::IW ][ ug.cId ];

    turbcom.rho  = ( * uturbf.q_ns  )[ IDX::IR  ][ ug.cId ];
    turbcom.ke   = ( * uturbf.q  )[ IKE ][ ug.cId ];
    turbcom.kw   = ( * uturbf.q  )[ IKW ][ ug.cId ];
    turbcom.visl = ( * uturbf.visl )[ 0 ][ ug.cId ];
    turbcom.vist = ( * uturbf.vist )[ 0 ][ ug.cId ];
    gcom.cvol    = ( * ug.cvol )[ ug.cId ];
    turbcom.cross_term = ( * uturbf.cross )[ 0 ][ ug.cId ];

    if ( turbcom.ibld == 1 )
    {
        turbcom.bld = ( * uturbf.bld )[ 0 ][ ug.cId ];
    }

    if ( turbcom.des_model ) 
    {
        turbcom.len_scale = ( * uturbf.len_scale )[ 0 ][ ug.cId ];
    }
}

void UTurbSrcFlux::PrepareCellValue1Equ()
{
    turbcom.dudx = ( * uturbf.dqdx_ns )[ IDX::IU ][ ug.cId ];
    turbcom.dudy = ( * uturbf.dqdy_ns )[ IDX::IU ][ ug.cId ];
    turbcom.dudz = ( * uturbf.dqdz_ns )[ IDX::IU ][ ug.cId ];
    turbcom.dvdx = ( * uturbf.dqdx_ns )[ IDX::IV ][ ug.cId ];
    turbcom.dvdy = ( * uturbf.dqdy_ns )[ IDX::IV ][ ug.cId ];
    turbcom.dvdz = ( * uturbf.dqdz_ns )[ IDX::IV ][ ug.cId ];
    turbcom.dwdx = ( * uturbf.dqdx_ns )[ IDX::IW ][ ug.cId ];
    turbcom.dwdy = ( * uturbf.dqdy_ns )[ IDX::IW ][ ug.cId ];
    turbcom.dwdz = ( * uturbf.dqdz_ns )[ IDX::IW ][ ug.cId ];

    turbcom.rho  = ( * uturbf.q_ns  )[ IDX::IR  ][ ug.cId ];
    turbcom.nuet = ( * uturbf.q  )[ ISA ][ ug.cId ];
    turbcom.visl = ( * uturbf.visl )[ 0 ][ ug.cId ];
    turbcom.vist = ( * uturbf.vist )[ 0 ][ ug.cId ];
    gcom.cvol    = ( * ug.cvol )[ ug.cId ];
    turbcom.cross_term = ( * uturbf.cross )[ 0 ][ ug.cId ];

    turbcom.dqdxSa = ( * uturbf.dqdx )[ ISA ][ ug.cId ];
    turbcom.dqdySa = ( * uturbf.dqdy )[ ISA ][ ug.cId ];
    turbcom.dqdzSa = ( * uturbf.dqdz )[ ISA ][ ug.cId ];

    turbcom.dist  = ( * uturbf.dist )[ ug.cId ];

    if ( turbcom.des_model ) 
    {
        turbcom.len_scale = ( * uturbf.len_scale )[ 0 ][  ug.cId ];
    }
    else
    {
        turbcom.len_scale = ( * uturbf.dist )[ ug.cId ];
    }
}

void UTurbSrcFlux::CalcLengthScaleSa()
{
    if ( turbcom.des_model == DES ) 
    {
        this->CalcLengthScaleOfSaDes();
    }
    else if ( turbcom.des_model == DDES ) 
    {
        this->CalcLengthScaleOfSaDdes();
    }
    else if ( turbcom.des_model == IDDES ) 
    {
        this->CalcLengthScaleOfSaIddes();
    }
    else
    {
        this->CalcLengthScaleOfWallDist();
    }
}

void UTurbSrcFlux::CalcLengthScaleSst()
{
    if ( turbcom.des_model == DES ) 
    {
        this->CalcLengthScaleOfSstDes();
    }
    else if ( turbcom.des_model == DDES ) 
    {
        this->CalcLengthScaleOfSstDdes();
    }
    else if ( turbcom.des_model == IDDES ) 
    {
        this->CalcLengthScaleOfSstIddes();
    }
    else
    {
        this->CalcLengthScaleOfWallDist();
    }
}

void UTurbSrcFlux::CalcLengthScaleOfSaDes()
{
    UnsGrid * grid = Zone::GetUnsGrid();
    MRField * len_scale = GetFieldPointer< MRField > ( grid, "len_scale" );

    int numberOfCells = grid->nCell;

    RealField lesLength( numberOfCells );

    CalcLengthLesOfSa( lesLength );

    RealField & wall_dist = grid->cellMesh->dist;

    for ( int cId = 0; cId < numberOfCells; ++ cId ) 
    {
        Real len_Rans = wall_dist[ cId ];
        Real len_Les  = lesLength[ cId ];

        ( * len_scale )[ 0 ][ cId ] = MIN( len_Rans, len_Les );
    }
}

void UTurbSrcFlux::CalcLengthScaleOfSstDes()
{
    UnsGrid * grid = Zone::GetUnsGrid();
    MRField * len_scale = GetFieldPointer< MRField > ( grid, "len_scale" );

    int numberOfCells = grid->nCell;

    RealField lesLength( numberOfCells );

    CalcLengthLesOfSst( lesLength );

    RealField & wall_dist = grid->cellMesh->dist;

    for ( int cId = 0; cId < numberOfCells; ++ cId ) 
    {
        Real ke = ( * uturbf.q  )[ IKE ][ cId ];
        Real kw = ( * uturbf.q  )[ IKW ][ cId ];

        Real len_Rans = sqrt( ABS( ke ) ) / ( turbcom.betas * kw + SMALL );
        Real len_Les  = lesLength[ cId ];

        ( * len_scale )[ 0 ][ cId ] = MIN( len_Rans, len_Les );
    }
}

void UTurbSrcFlux::CalcLengthScaleOfSaDdes()
{
    UnsGrid * grid = Zone::GetUnsGrid();
    MRField * len_scale = GetFieldPointer< MRField > ( grid, "len_scale" );

    int nCell = grid->nCell;

    RealField & wall_dist = grid->cellMesh->dist;

    RealField lesLength( nCell );
    CalcLengthLesOfSa( lesLength );

    for ( int cId = 0; cId < nCell; ++ cId ) 
    {
        Real rho = ( * uturbf.q_ns  )[ IDX::IR ][ cId ];
        Real visl = ( * uturbf.visl  )[ 0 ][ cId ];
        Real vist = ( * uturbf.vist  )[ 0 ][ cId ];

        Real dist  = ( * uturbf.dist )[ cId ];
        Real dist2 = SQR( dist );
        Real len_Rans = dist;
        Real len_Les  = lesLength[ cId ];

        turbcom.dudx = ( * uturbf.dqdx_ns )[ IDX::IU ][ ug.cId ];
        turbcom.dudy = ( * uturbf.dqdy_ns )[ IDX::IU ][ ug.cId ];
        turbcom.dudz = ( * uturbf.dqdz_ns )[ IDX::IU ][ ug.cId ];
        turbcom.dvdx = ( * uturbf.dqdx_ns )[ IDX::IV ][ ug.cId ];
        turbcom.dvdy = ( * uturbf.dqdy_ns )[ IDX::IV ][ ug.cId ];
        turbcom.dvdz = ( * uturbf.dqdz_ns )[ IDX::IV ][ ug.cId ];
        turbcom.dwdx = ( * uturbf.dqdx_ns )[ IDX::IW ][ ug.cId ];
        turbcom.dwdy = ( * uturbf.dqdy_ns )[ IDX::IW ][ ug.cId ];
        turbcom.dwdz = ( * uturbf.dqdz_ns )[ IDX::IW ][ ug.cId ];

        Real gradu2 = SQR( turbcom.dudx, turbcom.dudy, turbcom.dudz );
        Real gradv2 = SQR( turbcom.dvdx, turbcom.dvdy, turbcom.dvdz );
        Real gradw2 = SQR( turbcom.dwdx, turbcom.dwdy, turbcom.dwdz );
        Real grad2  = gradu2 + gradv2 + gradw2;

        Real nueff = turbcom.oreynolds * ( visl + vist ) / rho;
        Real rd    = nueff / ( turbcom.karm2 * dist2 * MAX( sqrt( grad2 ), 1.0e-20 ) );
        Real fd    = 1.0 - tanh( POWER3( 8.0 * rd ) );

        ( * len_scale )[ 0 ][ cId ] = len_Rans - fd * MAX( 0.0, len_Rans - len_Les );
    }
}

void UTurbSrcFlux::CalcLengthScaleOfSstDdes()
{
    UnsGrid * grid = Zone::GetUnsGrid();
    MRField * len_scale = GetFieldPointer< MRField > ( grid, "len_scale" );

    int nCell = grid->nCell;

    RealField & wall_dist = grid->cellMesh->dist;

    RealField lesLength( nCell );
    CalcLengthLesOfSst( lesLength );

    for ( int cId = 0; cId < nCell; ++ cId ) 
    {
        Real rho = ( * uturbf.q_ns  )[ IDX::IR ][ cId ];
        Real visl = ( * uturbf.visl  )[ 0 ][ cId ];
        Real vist = ( * uturbf.vist  )[ 0 ][ cId ];

        Real ke = ( * uturbf.q  )[ IKE ][ cId ];
        Real kw = ( * uturbf.q  )[ IKW ][ cId ];

        Real dist  = wall_dist[ cId ];
        Real dist2 = SQR( dist );

        Real len_Rans = sqrt( ABS( ke ) ) / ( turbcom.betas * kw + SMALL );
        Real len_Les  = lesLength[ cId ];

        turbcom.dudx = ( * uturbf.dqdx_ns )[ IDX::IU ][ ug.cId ];
        turbcom.dudy = ( * uturbf.dqdy_ns )[ IDX::IU ][ ug.cId ];
        turbcom.dudz = ( * uturbf.dqdz_ns )[ IDX::IU ][ ug.cId ];
        turbcom.dvdx = ( * uturbf.dqdx_ns )[ IDX::IV ][ ug.cId ];
        turbcom.dvdy = ( * uturbf.dqdy_ns )[ IDX::IV ][ ug.cId ];
        turbcom.dvdz = ( * uturbf.dqdz_ns )[ IDX::IV ][ ug.cId ];
        turbcom.dwdx = ( * uturbf.dqdx_ns )[ IDX::IW ][ ug.cId ];
        turbcom.dwdy = ( * uturbf.dqdy_ns )[ IDX::IW ][ ug.cId ];
        turbcom.dwdz = ( * uturbf.dqdz_ns )[ IDX::IW ][ ug.cId ];

        Real s11 = turbcom.dudx;
        Real s22 = turbcom.dvdy;
        Real s33 = turbcom.dwdz;
        Real s12 = half * ( turbcom.dudy + turbcom.dvdx );
        Real s13 = half * ( turbcom.dudz + turbcom.dwdx );
        Real s23 = half * ( turbcom.dvdz + turbcom.dwdy );
        Real w12 = half * ( turbcom.dudy - turbcom.dvdx );
        Real w13 = half * ( turbcom.dudz - turbcom.dwdx );
        Real w23 = half * ( turbcom.dvdz - turbcom.dwdy );

        Real sij2  = two * ( SQR( s11, s22, s33 ) + two * SQR( s12, s13, s23 ) );
        Real vort2 = four * SQR( w12, w13, w23 );
        Real term  = MAX(  sqrt( half * ( sij2 + vort2 ) ),  1.0e-20  );
        Real nueff = turbcom.oreynolds * ( visl + vist ) / rho;

        Real rd = nueff / ( turbcom.karm2 * dist2 * term );
        Real fd = 1.0 - tanh( POWER3( 20.0 * rd ) );

        ( * len_scale )[ 0 ][ cId ] = len_Rans - fd * MAX( 0.0, len_Rans - len_Les );
    }
}

void UTurbSrcFlux::CalcLengthScaleOfSaIddes()
{
    UnsGrid * grid = Zone::GetUnsGrid();
    MRField * len_scale = GetFieldPointer< MRField > ( grid, "len_scale" );

    int nCell = grid->nCell;

    RealField & wall_dist = grid->cellMesh->dist;
    RealField & span = grid->cellMesh->span;

    RealField lesLength( nCell );

    CalcLengthLesOfSa( lesLength );

    RealField lowReynoldsCorr( nCell );
    CalcLowReynoldsCorrection( lowReynoldsCorr );

    Real ct   = 1.63;
    Real cl   = 3.55;

    for ( int cId = 0; cId < nCell; ++ cId ) 
    {
        Real rho = ( * uturbf.q_ns  )[ IDX::IR ][ cId ];
        Real visl = ( * uturbf.visl )[ 0 ][ cId ];
        Real vist = ( * uturbf.vist )[ 0 ][ cId ];

        Real dist  = wall_dist[ cId ];
        Real dist2 = SQR( dist );
        Real len_Rans = dist;
        Real len_Les  = lesLength[ cId ];

        Real cell_span = span[ cId ];

        Real alf  = 0.25 - dist / cell_span;
        Real fb   = MIN(  2.0 * exp( - 9.0 * alf * alf ),  1.0  );

        turbcom.dudx = ( * uturbf.dqdx_ns )[ IDX::IU ][ ug.cId ];
        turbcom.dudy = ( * uturbf.dqdy_ns )[ IDX::IU ][ ug.cId ];
        turbcom.dudz = ( * uturbf.dqdz_ns )[ IDX::IU ][ ug.cId ];
        turbcom.dvdx = ( * uturbf.dqdx_ns )[ IDX::IV ][ ug.cId ];
        turbcom.dvdy = ( * uturbf.dqdy_ns )[ IDX::IV ][ ug.cId ];
        turbcom.dvdz = ( * uturbf.dqdz_ns )[ IDX::IV ][ ug.cId ];
        turbcom.dwdx = ( * uturbf.dqdx_ns )[ IDX::IW ][ ug.cId ];
        turbcom.dwdy = ( * uturbf.dqdy_ns )[ IDX::IW ][ ug.cId ];
        turbcom.dwdz = ( * uturbf.dqdz_ns )[ IDX::IW ][ ug.cId ];

        Real gradu2 = SQR( turbcom.dudx, turbcom.dudy, turbcom.dudz );
        Real gradv2 = SQR( turbcom.dvdx, turbcom.dvdy, turbcom.dvdz );
        Real gradw2 = SQR( turbcom.dwdx, turbcom.dwdy, turbcom.dwdz );
        Real grad2  = gradu2 + gradv2 + gradw2;

        Real coef  = 1.0 / ( turbcom.karm2 * dist2 * MAX( sqrt( grad2 ), 1.0e-10 ) );

        Real nuet = turbcom.oreynolds * vist / rho;
        Real rdt  = nuet * coef;

        Real fdt  = 1.0 - tanh( POWER3( 8.0 * rdt ) );
        Real fdb  = MAX( ( 1.0 - fdt ), fb );

        Real nuel = turbcom.oreynolds * visl / rho;
        Real rdl  = nuel * coef;

        Real fl   = tanh(  pow( cl * cl * rdl, 10.0 )  );
        Real ft   = tanh(  POWER3( ct * ct * rdt )  );

        Real fe2  = 1.0 - MAX( ft, fl );

        Real fe1 = 0.0;
        if ( alf < 0.0 )
        {
            fe1 = 2.0 * exp( - 9.0 * alf * alf );
        }
        else
        {
            fe1 = 2.0 * exp( - 11.09 * alf * alf );
        }

        Real cfix = lowReynoldsCorr[ cId ];
        Real fe = fe2 * MAX( ( fe1 - 1.0 ), 0.0 ) * cfix;

        ( * len_scale )[ 0 ][ cId ]  =  fdb * ( 1.0 + fe ) * len_Rans  +  ( 1.0 - fdb ) * len_Les;
    }
}

void UTurbSrcFlux::CalcLengthScaleOfSstIddes()
{
    UnsGrid * grid = Zone::GetUnsGrid();
    MRField * len_scale = GetFieldPointer< MRField > ( grid, "len_scale" );

    int nCell = grid->nCell;

    RealField & wall_dist = grid->cellMesh->dist;
    RealField & span = grid->cellMesh->span;

    RealField lesLength( nCell );

    CalcLengthLesOfSst( lesLength );
    //sa model
    //Real ct   = 1.63;
    //Real cl   = 3.55;

    //sst model
    Real ct = 1.87;
    Real cl = 5.0;

    Real cdt1 = 20.0;
    Real cdt2 = 3.0;

    for ( int cId = 0; cId < nCell; ++ cId ) 
    {
        Real rho = ( * uturbf.q_ns  )[ IDX::IR ][ cId ];
        Real visl = ( * uturbf.visl )[ 0 ][ cId ];
        Real vist = ( * uturbf.vist )[ 0 ][ cId ];

        Real ke = ( * uturbf.q  )[ IKE ][ cId ];
        Real kw = ( * uturbf.q  )[ IKW ][ cId ];

        Real dist  = wall_dist[ cId ];
        Real dist2 = SQR( dist );

        Real cell_span = span[ cId ];

        Real len_Rans = sqrt( ABS( ke ) ) / ( turbcom.betas * kw + SMALL );
        Real len_Les  = lesLength[ cId ];

        turbcom.dudx = ( * uturbf.dqdx_ns )[ IDX::IU ][ ug.cId ];
        turbcom.dudy = ( * uturbf.dqdy_ns )[ IDX::IU ][ ug.cId ];
        turbcom.dudz = ( * uturbf.dqdz_ns )[ IDX::IU ][ ug.cId ];
        turbcom.dvdx = ( * uturbf.dqdx_ns )[ IDX::IV ][ ug.cId ];
        turbcom.dvdy = ( * uturbf.dqdy_ns )[ IDX::IV ][ ug.cId ];
        turbcom.dvdz = ( * uturbf.dqdz_ns )[ IDX::IV ][ ug.cId ];
        turbcom.dwdx = ( * uturbf.dqdx_ns )[ IDX::IW ][ ug.cId ];
        turbcom.dwdy = ( * uturbf.dqdy_ns )[ IDX::IW ][ ug.cId ];
        turbcom.dwdz = ( * uturbf.dqdz_ns )[ IDX::IW ][ ug.cId ];

        Real s11 = turbcom.dudx;
        Real s22 = turbcom.dvdy;
        Real s33 = turbcom.dwdz;
        Real s12 = half * ( turbcom.dudy + turbcom.dvdx );
        Real s13 = half * ( turbcom.dudz + turbcom.dwdx );
        Real s23 = half * ( turbcom.dvdz + turbcom.dwdy );
        Real w12 = half * ( turbcom.dudy - turbcom.dvdx );
        Real w13 = half * ( turbcom.dudz - turbcom.dwdx );
        Real w23 = half * ( turbcom.dvdz - turbcom.dwdy );

        Real sij2  = two * ( SQR( s11, s22, s33 ) + two * SQR( s12, s13, s23 ) );
        Real vort2 = four * SQR( w12, w13, w23 );
        Real term  = MAX( sqrt( half * ( sij2 + vort2 ) ),  1.0e-20  );

        Real coef = 1.0 / ( turbcom.karm2 * dist2 * term );

        Real alf  = 0.25 - dist / cell_span;
        Real fb   = MIN(  2.0 * exp( - 9.0 * alf * alf ),  1.0  );

        Real nuet = turbcom.oreynolds * vist / rho;
        Real rdt = nuet * coef;

        Real fdt = 1.0 - tanh(  pow( ( cdt1 * rdt ), cdt2 )  );
        //Real fdt  = 1.0 - tanh( POWER3( 8.0 * rdt ) );
        Real fdb  = MAX( ( 1.0 - fdt ), fb );

        Real nuel = turbcom.oreynolds * visl / rho;
        Real rdl  = nuel * coef;

        Real fl   = tanh(  pow( cl * cl * rdl, 10.0 )  );
        Real ft   = tanh(  POWER3( ct * ct * rdt )  );

        Real fe2  = 1.0 - MAX( ft, fl );

        Real fe1 = 0.0;
        if ( alf < 0.0 )
        {
            fe1 = 2.0 * exp( - 9.0 * alf * alf );
        }
        else
        {
            fe1 = 2.0 * exp( - 11.09 * alf * alf );
        }

        Real cfix = 1.0;
        Real fe = fe2 * MAX( ( fe1 - 1.0 ), 0.0 ) * cfix;

        ( * len_scale )[ 0 ][ cId ]  =  fdb * ( 1.0 + fe ) * len_Rans  +  ( 1.0 - fdb ) * len_Les;
    }
}


void UTurbSrcFlux::CalcLengthScaleOfWallDist()
{
    UnsGrid * grid = Zone::GetUnsGrid();
    MRField * len_scale = GetFieldPointer< MRField > ( grid, "len_scale" );

    RealField & wall_dist = grid->cellMesh->dist;
    int nCell = grid->nCell;

    for ( int cId = 0; cId < nCell; ++ cId ) 
    {
        ( * len_scale )[ 0 ][ cId ] = wall_dist[ cId ];
    }
}

void UTurbSrcFlux::CalcBlendingTerm()
{
    turbcom.ibld = 1;

    this->CalcCdkwmax();

    this->CalcBlendField();

    this->ModifyBlendingTerm();
}

void UTurbSrcFlux::CalcCdkwmax()
{
    //Calc maximum crossing diffusion term across flowfield
    turbcom.cdkwmax = 0.0;

    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        ug.cId = cId;
        turbcom.rho  = ( * uturbf.q_ns  )[ IDX::IR ][ ug.cId ];
        turbcom.ke   = ( * uturbf.q  )[ IKE ][ ug.cId ];
        turbcom.kw   = ( * uturbf.q  )[ IKW ][ ug.cId ];
        turbcom.cross_term = ( * uturbf.cross )[ 0 ][ ug.cId ];

        this->CalcCrossDiff();
    }

    turbcom.CalcCdkwmin();
}

void UTurbSrcFlux::CalcBlendField()
{
    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        ug.cId = cId;
        if ( cId == 11 )
        {
            int kkk = 1;
        }
        turbcom.rho  = ( * uturbf.q_ns  )[ IDX::IR ][ ug.cId ];
        turbcom.ke   = ( * uturbf.q  )[ IKE ][ ug.cId ];
        turbcom.kw   = ( * uturbf.q  )[ IKW ][ ug.cId ];

        turbcom.cross_term = ( * uturbf.cross )[ 0 ][ ug.cId ];

        turbcom.visl = ( * uturbf.visl )[ 0 ][ ug.cId ];
        turbcom.dist = ( * uturbf.dist )[ ug.cId ];

        turbcom.CalcCellBlendingTerm();

        ( * uturbf.bld )[ 0 ][ ug.cId ] = turbcom.bld;
    }
}

void UTurbSrcFlux::ModifyBlendingTerm()
{
    if ( turbcom.transition_model == ITReGama ) 
    {
    }

    IntField & bcType = ug.bcRecord->bcType;

    for ( int fId = 0; fId < ug.nBFace; ++ fId )
    {
        int bc_type = bcType[ fId ];

        if ( bc_type == BC::INTERFACE ) continue;

        int lc  = ( * ug.lcf ) [ fId ];
        int rc  = ( * ug.rcf ) [ fId ];

        ( * uturbf.bld )[ 0 ][ rc ] = ( * uturbf.bld )[ 0 ][ lc ];
    }
}

void UTurbSrcFlux::CalcCrossingTerm()
{
    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        ug.cId = cId;

        turbcom.dkedx = ( * uturbf.dqdx )[ IKE ][ ug.cId ];
        turbcom.dkedy = ( * uturbf.dqdy )[ IKE ][ ug.cId ];
        turbcom.dkedz = ( * uturbf.dqdz )[ IKE ][ ug.cId ];

        turbcom.dkwdx = ( * uturbf.dqdx )[ IKW ][ ug.cId ];
        turbcom.dkwdy = ( * uturbf.dqdy )[ IKW ][ ug.cId ];
        turbcom.dkwdz = ( * uturbf.dqdz )[ IKW ][ ug.cId ];

        turbcom.CalcCrossing();

        ( * uturbf.cross )[ 0 ][ cId ] = turbcom.cross_term;
    }
}

void UTurbSrcFlux::CalcCrossDiff()
{
    turbcom.CalcCrossDiff();
}

void CalcLengthLesOfSa( RealField & lesLengthField )
{
    UnsGrid * grid = Zone::GetUnsGrid();

    int nCell = grid->nCell;

    RealField lowReynoldsCorr( nCell );

    CalcLowReynoldsCorrection( lowReynoldsCorr );

    RealField lengthScale = GetLengthScale();

    for ( int cId = 0; cId < nCell; ++ cId ) 
    {
        Real term1 = lowReynoldsCorr[ cId ];
        Real term2 = lengthScale[ cId ];

        lesLengthField[ cId ] = turbcom.cdes * term1 * term2;
    }
}

void CalcLengthLesOfSst( RealField & lesLengthField )
{
    UnsGrid * grid = Zone::GetUnsGrid();
    int nCell = grid->nCell;

    RealField lengthScale = GetLengthScale();

    for ( int cId = 0; cId < nCell; ++ cId ) 
    {
        Real term1 = lengthScale[ cId ];
        turbcom.bld = ( * uturbf.bld )[ 0 ][ cId ];
        turbcom.cdes = ( 1.0 - turbcom.bld ) * turbcom.cdes_ke + turbcom.bld * turbcom.cdes_kw;

        lesLengthField[ cId ] = turbcom.cdes * term1;
    }
}

RealField GetLengthScale()
{
    if ( turbcom.des_model == IDDES )
    {
        RealField lenth_scale;
        CalcSubgridLengthScale( lenth_scale );
        return lenth_scale;
    }
    else
    {
        RealField & largestSpacing = GetLargestSpacing();
        return largestSpacing;
    }
}

void CalcLowReynoldsCorrection( RealField & lowReynoldsCorr )
{
    UnsGrid * grid = Zone::GetUnsGrid();

    for ( int cId = 0; cId < grid->nCell; ++ cId ) 
    {
        Real nuet = ( * uturbf.q )[ ISA ][ cId ];
        Real rho = ( * uturbf.q_ns  )[ IDX::IR ][ cId ];
        Real visl = ( * uturbf.visl  )[ 0 ][ cId ];
        Real olam = rho / ( visl + SMALL );
        Real xsi  = nuet * olam + SMALL;
        Real xsi3 = POWER3( xsi );
        Real xsi2 = SQR( xsi );

        turbcom.ft2 = 0.0; 
        if ( turbcom.ft2_flag )
        {
            turbcom.ft2 = turbcom.ct3 * exp( - turbcom.ct4 * xsi2 );
        }

        turbcom.fv1 = xsi3 / ( xsi3 + turbcom.cv13 );
        turbcom.fv2 = one - xsi / ( one + xsi * turbcom.fv1 );
        Real term1 = turbcom.cb1 / ( turbcom.cw1 * turbcom.karm2 * turbcom.fwStar );
        Real term2 = ( turbcom.ft2 + ( 1.0 - turbcom.ft2 ) * turbcom.fv2 );
        Real term3   =  1.0  -  term1 * term2;
        Real term4 =  turbcom.fv1 * MAX( 1.0e-10, ( 1.0 - turbcom.ft2 ) );

        lowReynoldsCorr[ cId ] = sqrt(  MIN( 1.0e2, term3 / term4 )  );
    }
}

void CalcSubgridLengthScale( RealField & lenth_scale )
{
    UnsGrid * grid = Zone::GetUnsGrid();

    CalcCellSpan( grid );

    RealField & wall_dist = grid->cellMesh->dist;
    RealField & largestSpacing = GetLargestSpacing();

    lenth_scale.resize( grid->nCell );
    
    for ( int cId = 0; cId < grid->nCell; ++ cId ) 
    {
        Real cw    = 0.15;
        Real dist = wall_dist[ cId ];
        Real span = largestSpacing[ cId ];

        Real subgrid_scale = MIN( cw * MAX( dist, span ),  span );

        lenth_scale[ cId ] = subgrid_scale;
    }
}

RealField & GetLargestSpacing()
{
    UnsGrid * grid = Zone::GetUnsGrid();
    CalcCellSpan( grid );
    return grid->cellMesh->span;
}

EndNameSpace