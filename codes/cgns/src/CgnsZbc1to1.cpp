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

#include "CgnsZbc1to1.h"
#include "CgnsBc1to1.h"
#include "CgnsBcBoco.h"
#include "CgnsZone.h"
#include "CgnsBase.h"
#include "CgnsFile.h"
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

CgnsZbc1to1::CgnsZbc1to1( CgnsZone & cgnsZone )
    : cgnsZone( cgnsZone )
{
    this->n1to1ToCreate = 0;
}

CgnsZbc1to1::~CgnsZbc1to1() = default;

int CgnsZbc1to1::GetN1to1() const
{
    return static_cast< int >( this->cgnsBc1to1s.size() );
}

void CgnsZbc1to1::AddCgns1To1BcRegion( std::unique_ptr< CgnsBc1to1 > cgnsBc1to1 )
{
    if ( cgnsBc1to1 == nullptr )
    {
        throw std::invalid_argument( "CgnsZbc1to1: cannot add a null interface" );
    }

    CgnsBc1to1 * bc1to1 = cgnsBc1to1.get();
    this->cgnsBc1to1s.push_back( std::move( cgnsBc1to1 ) );
    int id = this->cgnsBc1to1s.size();
    bc1to1->bcId = id;
}

CgnsBc1to1 & CgnsZbc1to1::GetCgnsBcRegion1to1( int i1to1 )
{
    return *this->cgnsBc1to1s.at( i1to1 );
}

const CgnsBc1to1 & CgnsZbc1to1::GetCgnsBcRegion1to1( int i1to1 ) const
{
    return *this->cgnsBc1to1s.at( i1to1 );
}

void CgnsZbc1to1::PrintZn1to1()
{
    std::cout << "   n1to1        = " << this->n1to1ToCreate << std::endl;
}

void CgnsZbc1to1::CreateCgnsZbc()
{
    for ( int i1to1 = this->GetN1to1(); i1to1 < this->n1to1ToCreate; ++ i1to1 )
    {
        this->AddCgns1To1BcRegion( std::make_unique< CgnsBc1to1 >( &this->cgnsZone ) );
    }
}

void CgnsZbc1to1::ConvertToInnerDataStandard()
{
    for ( int i1to1 = 0; i1to1 < this->GetN1to1(); ++ i1to1 )
    {
        CgnsBc1to1 & cgnsBc1to1 = this->GetCgnsBcRegion1to1( i1to1 );
        cgnsBc1to1.ConvertToInnerDataStandard();
    }
}

void CgnsZbc1to1::ReadZn1to1( int n1to1 )
{
    this->n1to1ToCreate = n1to1;
    this->PrintZn1to1();
}

void CgnsZbc1to1::ReadZn1to1()
{
    int fileId = cgnsZone.cgnsBase.cgnsFile->fileId;
    int baseId = cgnsZone.cgnsBase.baseId;
    int zId = cgnsZone.zId;

    // find out how many general interfaces there are in this zone
    // the following is the number of structured grid interface
    cg_n1to1( fileId, baseId, zId, & this->n1to1ToCreate );

    this->PrintZn1to1();
}

void CgnsZbc1to1::ReadCgnsZbc1to1()
{
    this->ReadZn1to1();
    this->CreateCgnsZbc();

    for ( int i1to1 = 0; i1to1 < this->GetN1to1(); ++ i1to1 )
    {
        CgnsBc1to1 & cgnsBc1to1 = this->GetCgnsBcRegion1to1( i1to1 );
        cgnsBc1to1.ReadCgnsBc1To1();
    }
}
void CgnsZbc1to1::DumpCgnsZbc1to1()
{
    this->PrintZn1to1();

    for ( int i1to1 = 0; i1to1 < this->GetN1to1(); ++ i1to1 )
    {
        CgnsBc1to1 & cgnsBc1to1 = this->GetCgnsBcRegion1to1( i1to1 );
        cgnsBc1to1.DumpCgnsBc1To1();
    }
}


#endif
EndNameSpace
