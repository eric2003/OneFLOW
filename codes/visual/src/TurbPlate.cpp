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

#include "TurbPlate.h"
#include "NodeField.h"
#include "Zone.h"
#include "ZoneState.h"
#include "TaskState.h"
#include "UnsGrid.h"
#include "NodeMesh.h"
#include "FaceTopo.h"
#include "BcRecord.h"
#include "PointSearch.h"
#include "HXMath.h"
#include "Parallel.h"
#include "PIO.h"
#include "DataBook.h"
#include "DataBase.h"
#include "NsCom.h"
#include "Boundary.h"
#include <algorithm>
#include <iomanip>
using namespace std;

BeginNameSpace( ONEFLOW )

void SetTurbPlateTask()
{
    REGISTER_DATA_CLASS( CreateTurbPlateTask  );
}

void CreateTurbPlateTask( StringField & data )
{
    TurbFlatPlateTask * task = new TurbFlatPlateTask();
    TaskState::task = task;
}

TurbVelCut::TurbVelCut()
{
    sliceInfo.AddSlice( 0.19071e+06 / nscom.reynolds, X_DIR, Y_DIR );
    sliceInfo.AddSlice( 0.10643e+07 / nscom.reynolds, X_DIR, Y_DIR );
    sliceInfo.AddSlice( 0.27034e+07 / nscom.reynolds, X_DIR, Y_DIR );
    sliceInfo.AddSlice( 0.49981e+07 / nscom.reynolds, X_DIR, Y_DIR );
    sliceInfo.AddSlice( 0.76206e+07 / nscom.reynolds, X_DIR, Y_DIR );
    sliceInfo.AddSlice( 0.10274e+08 / nscom.reynolds, X_DIR, Y_DIR );
    nameList.push_back( "q" );
    nameList.push_back( "visl" );
    this->Init();
}

TurbVelCut::~TurbVelCut()
{
}


void TurbVelCut::Dump()
{
    string velocityFile = "results/turbplateflow.dat";

    fstream file;
    PIO::ParallelOpenPrj( file, velocityFile, ios_base::out );
    StringField title;
    title.push_back( "title=\"THE FLOW FIELD OF ONEFLOW\"" );
    title.push_back( "variables=" );
    title.push_back( "\"y+\"" );
    title.push_back( "\"u+\"" );

    for ( UInt i = 0; i < title.size(); ++ i )
    {
        file << title[ i ] << endl;
    }

    size_t nSlice = sliceData.size();
    for ( int i = 0; i < nSlice; ++ i )
    {
        LamData * lamData = sliceData[ i ];
        this->Dump( lamData, file, sliceInfo.dir2[ i ] );
    }

    PIO::Close( file );
}

void TurbVelCut::Dump( LamData * lamData, fstream & file, int axis )
{
    int nNode = lamData->GetNNode();

    file << "zone  i = " << nNode << " \n";

    int wordWidth = 20;

    RealField & x = lamData->data[ 0 ]->x;
    RealField & y = lamData->data[ 0 ]->y;
    RealField & z = lamData->data[ 0 ]->z;

    lamData->SortDataByAxis( axis );

    RealField2D & qdata   = lamData->data[ 0 ]->slicedata;
    RealField2D & visdata = lamData->data[ 1 ]->slicedata;

    int ywId = lamData->FindYIndex();
    Real xw = x[ ywId ];
    Real yw = y[ ywId ];
    Real zw = z[ ywId ];

    Real rw = qdata[ 0 ][ ywId ];
    Real uw = qdata[ 1 ][ ywId ];
    Real vw = qdata[ 2 ][ ywId ];
    Real ww = qdata[ 3 ][ ywId ];
    Real pw = qdata[ 4 ][ ywId ];

    Real visw = visdata[ 0 ][ ywId ];

    for ( int iNode = 0; iNode < nNode; ++ iNode )
    {
        Real xm = x[ iNode ];
        Real ym = y[ iNode ];
        Real zm = z[ iNode ];

        Real rm = qdata[ 0 ][ iNode ];
        Real um = qdata[ 1 ][ iNode ];
        Real vm = qdata[ 2 ][ iNode ];
        Real wm = qdata[ 3 ][ iNode ];
        Real pm = qdata[ 4 ][ iNode ];

        Real vis = visdata[ 0 ][ iNode ];

        Real dudy = uw / yw;
        Real tauw = visw * dudy;

        Real utau = sqrt( tauw / ( rw * nscom.reynolds ) );

        //注意这里的定义
        Real up = um / utau;
        Real yp = utau * ym * nscom.reynolds / ( vis / rm );

        file << setiosflags( ios::left );
        file << setiosflags( ios::scientific );
        file << setprecision( 10 );
        file << setw( wordWidth ) << yp << setw( wordWidth ) << up << endl;
    }
}

TurbFriCut::TurbFriCut()
{
    sliceInfo.AddSlice( 5.0e-7, Y_DIR, X_DIR );
    nameList.push_back( "q" );
    nameList.push_back( "visl" );
    this->Init();
}

TurbFriCut::~TurbFriCut()
{
}

void TurbFriCut::Dump()
{
    string frictionFile = "results/turbplate_cf.dat";

    fstream file;
    PIO::ParallelOpenPrj( file, frictionFile, ios_base::out );
    StringField title;
    title.push_back( "title=\"THE FLOW FIELD OF ONEFLOW\"" );
    title.push_back( "variables=" );
    title.push_back( "\"x\"" );
    title.push_back( "\"cf\"" );

    for ( UInt i = 0; i < title.size(); ++ i )
    {
        file << title[ i ] << endl;
    }

    size_t nSlice = sliceData.size();
    for ( int i = 0; i < nSlice; ++ i )
    {
        LamData * lamData = sliceData[ i ];
        this->Dump( lamData, file, sliceInfo.dir2[ i ] );
    }

    PIO::Close( file );
}

void TurbFriCut::Dump( LamData * lamData, fstream & file, int axis )
{
    int nNode = lamData->GetNNode();

    file << "zone  i = " << nNode << " \n";

    int wordWidth = 20;

    RealField & x = lamData->data[ 0 ]->x;
    RealField & y = lamData->data[ 0 ]->y;
    RealField & z = lamData->data[ 0 ]->z;

    lamData->SortDataByAxis( axis );

    RealField2D & qdata   = lamData->data[ 0 ]->slicedata;
    RealField2D & visdata = lamData->data[ 1 ]->slicedata;

    Real vel_inf = 1.0;

    for ( int iNode = 0; iNode < nNode; ++ iNode )
    {
        Real xm = x[ iNode ];
        Real ym = y[ iNode ];
        Real zm = z[ iNode ];

        Real rm = qdata[ 0 ][ iNode ];
        Real um = qdata[ 1 ][ iNode ];
        Real vm = qdata[ 2 ][ iNode ];
        Real wm = qdata[ 3 ][ iNode ];
        Real pm = qdata[ 4 ][ iNode ];

        Real vis = visdata[ 0 ][ iNode ];

        Real dudy = um / ym;
        Real xf = rm * vel_inf * xm / vis * nscom.reynolds;
        Real cf = 2 * vis * dudy / nscom.reynolds;

        file << setiosflags( ios::left );
        file << setiosflags( ios::scientific );
        file << setprecision( 10 );
        file << setw( wordWidth ) << xf << setw( wordWidth ) << cf << endl;
    }
}

TurbFlatPlateTask::TurbFlatPlateTask()
{
    velCut = new TurbVelCut();
    friCut = new TurbFriCut();
}

TurbFlatPlateTask::~TurbFlatPlateTask()
{
    delete velCut;
    delete friCut;
}

void TurbFlatPlateTask::Run()
{
    this->OutProfile( velCut );
    this->OutProfile( friCut );
}

void TurbFlatPlateTask::OutProfile( CuttingClass * cut )
{
    for ( int zId = 0; zId < ZoneState::nZones; ++ zId )
    {
        ZoneState::zid = zId;
        if ( ! ZoneState::IsValidZone( zId ) ) continue;
        cut->Slice();
    }

    cut->Swap();

    if ( Parallel::pid != Parallel::serverid ) return;

    cut->Dump();
}

EndNameSpace