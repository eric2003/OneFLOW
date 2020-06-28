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

#include "UTurbSpectrum.h"
#include "NsCom.h"
#include "TurbCom.h"
#include "TurbLusgs.h"
#include "UTurbCom.h"
#include "UNsCom.h"
#include "UCom.h"
#include "Zone.h"
#include "Grid.h"
#include "UnsGrid.h"
#include "FaceTopo.h"
#include "FaceMesh.h"
#include "CellMesh.h"
#include "HXMath.h"
#include "DataBase.h"
#include "FieldBase.h"
#include "Unsteady.h"
#include "UsdData.h"
#include "Ctrl.h"
#include "NsIdx.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )


UTurbSpectrum::UTurbSpectrum()
{
}

UTurbSpectrum::~UTurbSpectrum()
{
}

void UTurbSpectrum::ReadTmp()
{
    static int iii = 0;
    if ( iii ) return;
    iii = 1;
    fstream file;
    file.open( "turbflowsrc.dat", ios_base::in | ios_base::binary );
    if ( ! file ) exit( 0 );

    uturbf.Init();

       for ( int cId = 0; cId < ug.nTCell; ++ cId )
    {
        for ( int iEqu = 0; iEqu < 5; ++ iEqu )
        {
            //Real tmp1;
            //file.read( reinterpret_cast< char * >( & tmp1 ), sizeof( double ) );
            //( * uturbf.q_ns )[ iEqu ][ cId ] = tmp1;
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

    vector< Real > tmp1( ug.nTCell ), tmp2( ug.nTCell );

       for ( int cId = 0; cId < ug.nTCell; ++ cId )
    {
        tmp1[ cId ] = ( * uturbf.timestep )[ 0 ][ cId ];
    }

       for ( int cId = 0; cId < ug.nTCell; ++ cId )
    {
        file.read( reinterpret_cast< char * >( & ( * uturbf.timestep )[ 0 ][ cId ] ), sizeof( double ) );
    }

       for ( int cId = 0; cId < ug.nTCell; ++ cId )
    {
        tmp2[ cId ] = ( * uturbf.timestep )[ 0 ][ cId ];
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


void UTurbSpectrum::CalcSpectrum()
{
    ug.Init();
    turblu.Init();
    turbsp.Init();
    uturbf.Init();
    //ReadTmp();
    this->ZeroSpectrum();
    if ( turbcom.nEqu == 1 )
    {
        this->CalcSpectrum1Equ();
    }
    else if ( turbcom.nEqu >= 2 )
    {
        this->CalcSpectrum2Equ();
    }
}

void UTurbSpectrum::ZeroSpectrum()
{
    ( * uturbf.impsr ) = 0;
    ( * uturbf.matrix_l ) = 0;
    ( * uturbf.matrix_r ) = 0;
}

void UTurbSpectrum::CalcSpectrum1Equ()
{
    this->CalcUnsteadySpectrum();

    for ( int fId = 0; fId < ug.nFace; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        this->PrepareFaceValue1Equ();

        this->CalcFaceSpectrum1Equ();

        this->UpdateSpectrumRadius();
    }
}

void UTurbSpectrum::CalcSpectrum2Equ()
{
    this->CalcUnsteadySpectrum();

    for ( int fId = 0; fId < ug.nFace; ++ fId )
    {
        ug.fId = fId;
        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        this->PrepareFaceValue2Equ();

        this->CalcFaceSpectrum2Equ();

        this->UpdateSpectrumRadius();
    }
}

void UTurbSpectrum::CalcUnsteadySpectrum()
{
    if ( ctrl.idualtime == 0 )//单时间步，注意:是usd.sp2!
    {
        for ( int cId = 0; cId < ug.nCell; ++ cId )
        {
            Real vol = ( * ug.cvol )[ cId ];
            Real ts  = ( * uturbf.timestep )[ 0 ][ cId ] * turbcom.turb_cfl_ratio;
            Real unsteadyTerm = (  usd.sp2 / ts ) * vol;

            for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
            {
                ( * uturbf.impsr )[ iEqu ][ cId ] += unsteadyTerm;
            }
        }
    }
    else
    {
        for ( int cId = 0; cId < ug.nCell; ++ cId )
        {
            Real vol = ( * ug.cvol )[ cId ];
            Real ts  = ( * unsf.timestep )[ 0 ][ cId ] * turbcom.turb_cfl_ratio;
            Real unsteadyTerm = (  usd.sp1 / ts + usd.sp2 / ctrl.pdt1 ) * vol;

            for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
            {
                ( * uturbf.impsr )[ iEqu ][ cId ] += unsteadyTerm;
            }

        }
    }
}


void UTurbSpectrum::UpdateSpectrumRadius()
{
    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        ( * uturbf.impsr )[ iEqu ][ ug.lc ] += turbsp.radius1[ iEqu ];
        ( * uturbf.impsr )[ iEqu ][ ug.rc ] += turbsp.radius2[ iEqu ];
    }

    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        ( * uturbf.matrix_l )[ iEqu ][ ug.fId ] += turbsp.matrix1[ iEqu ];
        ( * uturbf.matrix_r )[ iEqu ][ ug.fId ] += turbsp.matrix2[ iEqu ];
    }
}

void UTurbSpectrum::PrepareFaceValue1Equ()
{
    gcom.xfn   = ( * ug.xfn   )[ ug.fId ];
    gcom.yfn   = ( * ug.yfn   )[ ug.fId ];
    gcom.zfn   = ( * ug.zfn   )[ ug.fId ];
    gcom.vfn   = ( * ug.vfn   )[ ug.fId ];
    gcom.farea = ( * ug.farea )[ ug.fId ];

    gcom.cvol1 = ( * ug.cvol )[ ug.lc ];
    gcom.cvol2 = ( * ug.cvol )[ ug.rc ];

    for ( int iEqu = 0; iEqu < nscom.nEqu; ++ iEqu )
    {
        turbsp.q1_ns[ iEqu ] = ( * uturbf.q_ns )[ iEqu ][ ug.lc ];
        turbsp.q2_ns[ iEqu ] = ( * uturbf.q_ns )[ iEqu ][ ug.rc ];
    }

    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turbsp.q1[ iEqu ] = ( * uturbf.q )[ iEqu ][ ug.lc ];
        turbsp.q2[ iEqu ] = ( * uturbf.q )[ iEqu ][ ug.rc ];
    }

    turbcom.visl1 = ( * uturbf.visl )[ 0 ][ ug.lc ];
    turbcom.visl2 = ( * uturbf.visl )[ 0 ][ ug.rc ];
}

void UTurbSpectrum::PrepareFaceValue2Equ()
{
    gcom.xfn   = ( * ug.xfn   )[ ug.fId ];
    gcom.yfn   = ( * ug.yfn   )[ ug.fId ];
    gcom.zfn   = ( * ug.zfn   )[ ug.fId ];
    gcom.vfn   = ( * ug.vfn   )[ ug.fId ];
    gcom.farea = ( * ug.farea )[ ug.fId ];

    gcom.cvol1 = ( * ug.cvol )[ ug.lc ];
    gcom.cvol2 = ( * ug.cvol )[ ug.rc ];

    for ( int iEqu = 0; iEqu < nscom.nEqu; ++ iEqu )
    {
        turbsp.q1_ns[ iEqu ] = ( * uturbf.q_ns )[ iEqu ][ ug.lc ];
        turbsp.q2_ns[ iEqu ] = ( * uturbf.q_ns )[ iEqu ][ ug.rc ];
    }

    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turbsp.q1[ iEqu ] = ( * uturbf.q )[ iEqu ][ ug.lc ];
        turbsp.q2[ iEqu ] = ( * uturbf.q )[ iEqu ][ ug.rc ];
    }

    turbcom.visl1 = ( * uturbf.visl )[ 0 ][ ug.lc ];
    turbcom.visl2 = ( * uturbf.visl )[ 0 ][ ug.rc ];

    turbcom.vist1 = ( * uturbf.vist )[ 0 ][ ug.lc ];
    turbcom.vist2 = ( * uturbf.vist )[ 0 ][ ug.rc ];

    turbcom.bld1 = ( * uturbf.bld )[ 0 ][ ug.lc ];
    turbcom.bld2 = ( * uturbf.bld )[ 0 ][ ug.rc ];
}

EndNameSpace