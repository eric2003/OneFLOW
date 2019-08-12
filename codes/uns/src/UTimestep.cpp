/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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

#include "UTimestep.h"
#include "NsCom.h"
#include "UNsCom.h"
#include "UCom.h"
#include "TurbCom.h"
#include "UTurbCom.h"
#include "NsCtrl.h"
#include "Zone.h"
#include "Grid.h"
#include "UnsGrid.h"
#include "FaceTopo.h"
#include "FaceMesh.h"
#include "CellMesh.h"
#include "HXMath.h"
#include "DataBase.h"
#include "FieldBase.h"
#include "Iteration.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

void InitUnsField()
{
    UnsGrid * grid = Zone::GetUnsGrid();
    unsf.invsr = ONEFLOW::GetFieldPointer< MRField >( grid, "invsr" );
    unsf.vissr = ONEFLOW::GetFieldPointer< MRField >( grid, "vissr" );
    unsf.gama  = ONEFLOW::GetFieldPointer< MRField >( grid, "gama" );
    unsf.q  = ONEFLOW::GetFieldPointer< MRField >( grid, "q" );
    //unsf.q1 = ONEFLOW::GetFieldPointer< MRField >( grid, "q1" );
    //unsf.q2 = ONEFLOW::GetFieldPointer< MRField >( grid, "q2" );
    unsf.dq = ONEFLOW::GetFieldPointer< MRField >( grid, "dq" );
    unsf.dqdx = ONEFLOW::GetFieldPointer< MRField >( grid, "dqdx" );
    unsf.dqdy = ONEFLOW::GetFieldPointer< MRField >( grid, "dqdy" );
    unsf.dqdz = ONEFLOW::GetFieldPointer< MRField >( grid, "dqdz" );

    unsf.visl = ONEFLOW::GetFieldPointer< MRField >( grid, "visl" );
    unsf.vist = ONEFLOW::GetFieldPointer< MRField >( grid, "vist" );
    //unsf.timestep = ONEFLOW::GetFieldPointer< MRField >( grid, "timestep" );
    unsf.res = ONEFLOW::GetFieldPointer< MRField >( grid, "res" );
    unsf.res1 = ONEFLOW::GetFieldPointer< MRField >( grid, "res1" );
    unsf.res2 = ONEFLOW::GetFieldPointer< MRField >( grid, "res2" );

    //unsf.tempr = ONEFLOW::GetFieldPointer< MRField >( grid, "tempr" );
}

void InitTimeStepUns()
{
    UnsGrid * grid = Zone::GetUnsGrid();

    InitUnsField();

    FaceTopo * faceTopo = grid->faceTopo;
    ug.lcf = & faceTopo->lCell;
    ug.rcf = & faceTopo->rCell;

    FaceMesh * faceMesh = grid->faceMesh;
    CellMesh * cellMesh = grid->cellMesh;

    ug.fnx = & faceMesh->xfn;
    ug.fny = & faceMesh->yfn;
    ug.fnz = & faceMesh->zfn;
    ug.vfn = & faceMesh->vfn;
    ug.farea = & faceMesh->area;

    ug.vfx = & faceMesh->vfx;
    ug.vfy = & faceMesh->vfy;
    ug.vfz = & faceMesh->vfz;

    ug.xfc = & faceMesh->xfc;
    ug.yfc = & faceMesh->yfc;
    ug.zfc = & faceMesh->zfc;

    ug.ccx = & cellMesh->xcc;
    ug.ccy = & cellMesh->ycc;
    ug.ccz = & cellMesh->zcc;

    ug.cvol  = & cellMesh->vol;
    ug.cvol1 = & cellMesh->vol;
    ug.cvol2 = & cellMesh->vol;
    nscom.Init();
}

UTimestep::UTimestep()
{
    ;
}

UTimestep::~UTimestep()
{
    ;
}

void UTimestep::Init()
{
    ug.Init();
    unsf.Init();
}

void UTimestep::ReadTmp()
{
    static int iii = 0;
    if ( iii ) return;
    iii = 1;
    fstream file;
    file.open( "nsflow.dat", ios_base::in | ios_base::binary );
    if ( ! file ) exit( 0 );

    unsf.Init();

       for ( int cId = 0; cId < ug.nTCell; ++ cId )
    {
        for ( int iEqu = 0; iEqu < 5; ++ iEqu )
        {
            file.read( reinterpret_cast< char * >( & ( * unsf.q )[ iEqu ][ cId ] ), sizeof( double ) );
        }
    }

       for ( int cId = 0; cId < ug.nTCell; ++ cId )
    {
        file.read( reinterpret_cast< char * >( & ( * unsf.visl )[ 0 ][ cId ] ), sizeof( double ) );
    }

       for ( int cId = 0; cId < ug.nTCell; ++ cId )
    {
        file.read( reinterpret_cast< char * >( & ( * unsf.vist )[ 0 ][ cId ] ), sizeof( double ) );
    }

    vector< Real > tmp1( ug.nTCell ), tmp2( ug.nTCell );

       for ( int cId = 0; cId < ug.nTCell; ++ cId )
    {
        tmp1[ cId ] = ( * unsf.timestep )[ 0 ][ cId ];
    }

       for ( int cId = 0; cId < ug.nTCell; ++ cId )
    {
        file.read( reinterpret_cast< char * >( & ( * unsf.timestep )[ 0 ][ cId ] ), sizeof( double ) );
    }

       for ( int cId = 0; cId < ug.nTCell; ++ cId )
    {
        tmp2[ cId ] = ( * unsf.timestep )[ 0 ][ cId ];
    }

    for ( int iCell = 0; iCell < ug.nTCell; ++ iCell )
    {
        for ( int iEqu = 0; iEqu < nscom.nTModel; ++ iEqu )
        {
            file.read( reinterpret_cast< char * >( & ( * unsf.tempr )[ iEqu ][ iCell ] ), sizeof( double ) );
        }
    }

    turbcom.Init();
    uturbf.Init();
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


void UTimestep::CmpTimestep()
{
    this->Init();

    //ReadTmp();

    this->CmpCfl();

    this->CmpSpectrumField();

    if ( nscom.timestepModel == 0 )
    {
        this->CmpLocalTimestep();
    }
    else if ( nscom.timestepModel == 1 )
    {
        this->CmpGlobalTimestep();
    }
    else
    {
        this->CmpLgTimestep();
    }
}

void UTimestep::CmpLocalTimestep()
{
    this->CmpInvTimestep();

    this->CmpVisTimestep();

    this->ModifyTimestep();
}

void UTimestep::CmpInvTimestep()
{
    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        ug.cId  = cId;
        gcom.cvol  = ( * ug.cvol )[ cId ];
        nscom.invsr = ( * unsf.invsr )[ 0 ][ cId ];
        this->CmpCellInvTimestep();
        ( * unsf.timestep )[ 0 ][ cId ] = nscom.timestep;
    }
}

void UTimestep::CmpVisTimestep()
{
    bool flag = vis_model.vismodel > 0 && nscom.visTimestepModel > 0;
    if ( ! flag ) return;

    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        ug.cId  = cId;
        gcom.cvol  = ( * ug.cvol )[ cId ];
        nscom.vissr = ( * unsf.vissr )[ 0 ][ cId ];
        nscom.timestep = ( * unsf.timestep )[ 0 ][ cId ];
        this->CmpCellVisTimestep();
        ( * unsf.timestep )[ 0 ][ cId ] = nscom.timestep;
    }
}

void UTimestep::CmpSpectrumField()
{
    this->CmpInvSpectrumField();
    this->CmpVisSpectrumField();
}

void UTimestep::CmpInvSpectrumField()
{
    Grid * grid = Zone::GetGrid();

    MRField * invsr = ONEFLOW::GetFieldPointer< MRField >( grid, "invsr" );

    ONEFLOW::ZeroField( invsr, 1, grid->nCell );

    for ( int iFace = 0; iFace < grid->nFace; ++ iFace )
    {
        this->SetId( iFace );

        this->PrepareData();

        this->CmpFaceInvSpec();

        this->UpdateInvSpectrumField();
    }
}

void UTimestep::CmpVisSpectrumField()
{
    Grid * grid = Zone::GetGrid();

    MRField * vissr = ONEFLOW::GetFieldPointer< MRField >( grid, "vissr" );

    ONEFLOW::ZeroField( vissr, 1, grid->nCell );

    if ( vis_model.vismodel <= 0 ) return;

    for ( int iFace = 0; iFace < grid->nFace; ++ iFace )
    {
        this->SetId( iFace );

        this->PrepareVisData();

        this->CmpFaceVisSpec();

        this->UpdateVisSpectrumField();
    }
}

void UTimestep::SetId( int fId )
{
    ug.fId = fId;
    ug.lc = ( * ug.lcf )[ fId ];
    ug.rc = ( * ug.rcf )[ fId ];
}

void UTimestep::PrepareData()
{
    gcom.fnx   = ( * ug.fnx   )[ ug.fId ];
    gcom.fny   = ( * ug.fny   )[ ug.fId ];
    gcom.fnz   = ( * ug.fnz   )[ ug.fId ];
    gcom.vfn   = ( * ug.vfn   )[ ug.fId ];
    gcom.farea = ( * ug.farea )[ ug.fId ];

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.q1[ iEqu ] = ( * unsf.q )[ iEqu ][ ug.lc ];
        nscom.q2[ iEqu ] = ( * unsf.q )[ iEqu ][ ug.rc ];
    }

    nscom.gama1 = ( * unsf.gama  )[ 0 ][ ug.lc ];
    nscom.gama2 = ( * unsf.gama  )[ 0 ][ ug.rc ];
}

void UTimestep::PrepareVisData()
{
    gcom.fnx   = ( * ug.fnx   )[ ug.fId ];
    gcom.fny   = ( * ug.fny   )[ ug.fId ];
    gcom.fnz   = ( * ug.fnz   )[ ug.fId ];
    gcom.vfn   = ( * ug.vfn   )[ ug.fId ];
    gcom.farea = ( * ug.farea )[ ug.fId ];

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.q1[ iEqu ] = ( * unsf.q )[ iEqu ][ ug.lc ];
        nscom.q2[ iEqu ] = ( * unsf.q )[ iEqu ][ ug.rc ];
    }

    nscom.gama1 = ( * unsf.gama )[ 0 ][ ug.lc ];
    nscom.gama2 = ( * unsf.gama )[ 0 ][ ug.rc ];

    gcom.cvol1 = ( * ug.cvol )[ ug.lc ];
    gcom.cvol2 = ( * ug.cvol )[ ug.rc ];

    nscom.visl1   = ( * unsf.visl  )[ 0 ][ ug.lc ];
    nscom.visl2   = ( * unsf.visl  )[ 0 ][ ug.rc ];

    nscom.visl = half * ( nscom.visl1 + nscom.visl2 );

    nscom.vist1 = ( * unsf.vist )[ 0 ][ ug.lc ];
    nscom.vist2 = ( * unsf.vist )[ 0 ][ ug.rc ];

    nscom.vist = half * ( nscom.vist1 + nscom.vist2 );

    gcom.ccx1 = ( * ug.ccx )[ ug.lc ];
    gcom.ccx2 = ( * ug.ccx )[ ug.rc ];

    gcom.ccy1 = ( * ug.ccy )[ ug.lc ];
    gcom.ccy2 = ( * ug.ccy )[ ug.rc ];

    gcom.ccz1 = ( * ug.ccz )[ ug.lc ];
    gcom.ccz2 = ( * ug.ccz )[ ug.rc ];
}

void UTimestep::UpdateInvSpectrumField()
{
    ( * unsf.invsr )[ 0 ][ ug.lc ] += nscom.invsr;
    ( * unsf.invsr )[ 0 ][ ug.rc ] += nscom.invsr;
}

void UTimestep::UpdateVisSpectrumField()
{
    ( * unsf.vissr )[ 0 ][ ug.lc ] += nscom.vissr;
    ( * unsf.vissr )[ 0 ][ ug.rc ] += nscom.vissr;
}

void UTimestep::ModifyTimestep()
{
    this->CmpMinTimestep();
    if ( nscom.max_time_ratio <= 0 ) return;

    Real maxPermittedTimestep = nscom.max_time_ratio * nscom.minTimestep;
    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        ( * unsf.timestep )[ 0 ][ cId ] = MIN( ( * unsf.timestep )[ 0 ][ cId ], maxPermittedTimestep );
    }
}

void UTimestep::CmpGlobalTimestep()
{
    this->SetTimestep( ctrl.pdt );
}

void UTimestep::CmpMinTimestep()
{
    nscom.minTimestep = LARGE;
    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        nscom.minTimestep = MIN( ( * unsf.timestep )[ 0 ][ cId ], nscom.minTimestep );
    }
}

void UTimestep::SetTimestep( Real timestep )
{
    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        ( * unsf.timestep )[ 0 ][ cId ] = timestep;
    }
}

void UTimestep::CmpLgTimestep()
{
    this->CmpLocalTimestep();
    this->SetTimestep( nscom.minTimestep );
    ctrl.pdt = nscom.minTimestep;
}


EndNameSpace