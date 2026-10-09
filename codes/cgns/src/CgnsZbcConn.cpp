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

    OneFLOW is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "CgnsZbcConn.h"
#include "CgnsBcConn.h"
#include "CgnsBcBoco.h"
#include "CgnsZone.h"
#include "CgnsFile.h"
#include "CgnsBase.h"
#include "Boundary.h"
#include "StringUtils.h"
#include "Dimension.h"
#include "HXMath.h"
#include "HXStd.h"
#include "StrRegion.h"
#include "StrGrid.h"
#include "GridMediator.h"
#include "FaceSolver.h"
#include "BcRecord.h"
#include <iostream>
#include <utility>
#include <stdexcept>



BeginNameSpace( ONEFLOW )
#ifdef ENABLE_CGNS

CgnsZbcConn::CgnsZbcConn( CgnsZone & cgnsZone )
    : cgnsZone( cgnsZone )
{
    this->nConnToCreate = 0;
}

CgnsZbcConn::~CgnsZbcConn() = default;

int CgnsZbcConn::GetNConn() const
{
    return static_cast< int >( this->cgnsBcConns.size() );
}

void CgnsZbcConn::AddCgnsConnBcRegion( std::unique_ptr< CgnsBcConn > cgnsBcConn )
{
    if ( cgnsBcConn == nullptr )
    {
        throw std::invalid_argument( "CgnsZbcConn: cannot add a null connection" );
    }

    CgnsBcConn * bcConn = cgnsBcConn.get();
    this->cgnsBcConns.push_back( std::move( cgnsBcConn ) );
    int id = this->cgnsBcConns.size();
    bcConn->bcId = id;
}

CgnsBcConn & CgnsZbcConn::GetCgnsBc( int iConn )
{
    return *this->cgnsBcConns.at( iConn );
}

void CgnsZbcConn::CreateCgnsZbc()
{
    for ( int iConn = 0; iConn < this->nConnToCreate; ++ iConn )
    {
this->AddCgnsConnBcRegion( std::make_unique< CgnsBcConn >( &this->cgnsZone ) );
    }
}

void CgnsZbcConn::PrintZnconn()
{
    std::cout << "   nConn        = " << this->nConnToCreate << std::endl;
}

void CgnsZbcConn::ReadZnconn( int nConn )
{
    this->nConnToCreate = nConn;
    this->PrintZnconn();
}

void CgnsZbcConn::ReadZnconn()
{
    int fileId = cgnsZone.cgnsBase.cgnsFile->fileId;
    int baseId = cgnsZone.cgnsBase.baseId;
    int zId = cgnsZone.zId;

    cg_nconns( fileId, baseId, zId, & this->nConnToCreate );
    this->PrintZnconn();
}

void CgnsZbcConn::ReadCgnsZbcConn()
{
    this->ReadZnconn();
    this->CreateCgnsZbc();
    for ( int iConn = 0; iConn < this->GetNConn(); ++ iConn )
    {
        CgnsBcConn & cgnsBcConn = this->GetCgnsBc( iConn );
        cgnsBcConn.ReadCgnsBcConn();
    }
}

void CgnsZbcConn::DumpCgnsZbcConn()
{
    this->PrintZnconn();
    for ( int iConn = 0; iConn < this->GetNConn(); ++ iConn )
    {
        CgnsBcConn & cgnsBcConn = this->GetCgnsBc( iConn );
        cgnsBcConn.DumpCgnsBcConn();
    }
}

void CgnsZbcConn::SetPeriodicBc()
{
    for ( int iConn = 0; iConn < this->GetNConn(); ++ iConn )
    {
        CgnsBcConn & cgnsBcConn = this->GetCgnsBc( iConn );
        cgnsBcConn.SetPeriodicBc();
    }
}

void CgnsZbcConn::ConvertToInnerDataStandard()
{
    for ( int iConn = 0; iConn < this->GetNConn(); ++ iConn )
    {
        CgnsBcConn & cgnsBcConn = this->GetCgnsBc( iConn );
        cgnsBcConn.ConvertToInnerDataStandard();
    }
}


#endif
EndNameSpace
