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

#include "CgnsZbcConn.h"
#include "CgnsBcConn.h"
#include "CgnsBcBoco.h"
#include "CgnsZone.h"
#include "CgnsBase.h"
#include "Boundary.h"
#include "StrUtil.h"
#include "Dimension.h"
#include "HXMath.h"
#include "HXStd.h"
#include "StrRegion.h"
#include "StrGrid.h"
#include "GridMediator.h"
#include "FaceSolver.h"
#include "BcRecord.h"
#include <iostream>

using namespace std;

BeginNameSpace( ONEFLOW )
#ifdef ENABLE_CGNS

CgnsZbcConn::CgnsZbcConn( CgnsZone * cgnsZone )
{
    this->cgnsZone = cgnsZone;
    this->nConn = 0;
}

CgnsZbcConn::~CgnsZbcConn()
{
    for ( int iConn = 0; iConn < this->nConn; ++ iConn )
    {
        delete this->cgnsBcConns[ iConn ];
    }
}

void CgnsZbcConn::AddCgnsConnBcRegion( CgnsBcConn * cgnsBcConn )
{
    this->cgnsBcConns.push_back( cgnsBcConn );
    int id = this->cgnsBcConns.size();
    cgnsBcConn->bcId = id;
}

CgnsBcConn * CgnsZbcConn::GetCgnsBc( int iConn )
{
    return this->cgnsBcConns[ iConn ];
}

void CgnsZbcConn::CreateCgnsZbc()
{
    for ( int iConn = 0; iConn < this->nConn; ++ iConn )
    {
        CgnsBcConn * cgnsBcConn = new CgnsBcConn( this->cgnsZone );
        this->AddCgnsConnBcRegion( cgnsBcConn );
    }
}

void CgnsZbcConn::PrintZnconn()
{
    cout << "   nConn        = " << this->nConn << endl;
}

void CgnsZbcConn::ReadZnconn( int nConn )
{
    this->nConn = nConn;
    this->PrintZnconn();
}

void CgnsZbcConn::ReadZnconn()
{
    int fileId = cgnsZone->cgnsBase->fileId;
    int baseId = cgnsZone->cgnsBase->baseId;
    int zId = cgnsZone->zId;

    cg_nconns( fileId, baseId, zId, & this->nConn );
    this->PrintZnconn();
}

void CgnsZbcConn::ReadCgnsZbcConn()
{
    this->ReadZnconn();
    this->CreateCgnsZbc();
    for ( int iConn = 0; iConn < this->nConn; ++ iConn )
    {
        CgnsBcConn * cgnsBcConn = this->GetCgnsBc( iConn );
        cgnsBcConn->ReadCgnsBcConn();
    }
}

void CgnsZbcConn::SetPeriodicBc()
{
    for ( int iConn = 0; iConn < this->nConn; ++ iConn )
    {
        CgnsBcConn * cgnsBcConn = this->GetCgnsBc( iConn );
        cgnsBcConn->SetPeriodicBc();
    }
}

void CgnsZbcConn::ConvertToInnerDataStandard()
{
    for ( int iConn = 0; iConn < this->nConn; ++ iConn )
    {
        CgnsBcConn * cgnsBcConn = this->GetCgnsBc( iConn );
        cgnsBcConn->ConvertToInnerDataStandard();
    }
}


#endif
EndNameSpace