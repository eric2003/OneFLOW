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

#include "CgnsBase.h"
#include "CgnsFile.h"
#include "CgnsZone.h"
#include "CgnsZoneUtils.h"
#include "StringUtils.h"
#include "Dimension.h"
#include "CgnsFamilyBc.h"
#include "CgnsVariable.h"

#include <iostream>
#include <memory>
#include <stdexcept>
#include <utility>

BeginNameSpace( ONEFLOW )

#ifdef ENABLE_CGNS

CgnsBase::CgnsBase()
    : cgnsFile( nullptr ),
      baseId( 0 ),
      nZones( 0 ),
      celldim( 0 ),
      phydim( 0 )
{
}

CgnsBase::CgnsBase( CgnsFile * cgnsFile )
    : cgnsFile( cgnsFile ),
      baseId( 0 ),
      nZones( 0 ),
      celldim( 0 ),
      phydim( 0 )
{
}

CgnsBase::~CgnsBase() = default;

void CgnsBase::FreeZoneList()
{
    this->cgnsZones.clear();
    this->zoneNameMap.clear();
    this->nZones = 0;
}


CgnsZone * CgnsBase::GetCgnsZone( int iZone )
{
    // Zone indices are zero-based internally.
    if ( iZone < 0 || static_cast< size_t >( iZone ) >= this->cgnsZones.size() )
    {
        throw std::out_of_range( "CgnsBase::GetCgnsZone: zone index is out of range" );
    }

    return this->cgnsZones[ static_cast< size_t >( iZone ) ].get();
}

CgnsZone * CgnsBase::GetCgnsZoneByName( const std::string & zoneName )
{
    const auto iter = this->zoneNameMap.find( zoneName );
    if ( iter == this->zoneNameMap.end() )
    {
        throw std::out_of_range( "CgnsBase::GetCgnsZoneByName: unknown zone name '" + zoneName + "'" );
    }

    return this->GetCgnsZone( iter->second - 1 );
}

int CgnsBase::GetNZones()
{
    return static_cast< int >( this->cgnsZones.size() );
}

void CgnsBase::SetDefaultCgnsBaseBasicInfo()
{
    //this->celldim = Dim::dimension;
    //this->phydim  = Dim::dimension;

    this->celldim = THREE_D;
    this->phydim  = THREE_D;
  
    this->baseName = ONEFLOW::AddString( "Base", this->baseId );
}

void CgnsBase::AddCgnsZone( std::unique_ptr< CgnsZone > cgnsZone )
{
    if ( ! cgnsZone )
    {
        throw std::invalid_argument( "CgnsBase::AddCgnsZone: cannot add a null zone" );
    }

    CgnsZone * zone = cgnsZone.get();
    zone->zId = static_cast< int >( this->cgnsZones.size() + 1 );
    this->cgnsZones.push_back( std::move( cgnsZone ) );
}

void CgnsBase::AllocateAllCgnsZones()
{
    if ( this->nZones < 0 )
    {
        throw std::invalid_argument( "CgnsBase::AllocateAllCgnsZones: zone count cannot be negative" );
    }

    if ( ! this->cgnsZones.empty() )
    {
        throw std::logic_error( "CgnsBase::AllocateAllCgnsZones: zones have already been allocated" );
    }

    for ( int iZone = 0; iZone < this->nZones; ++ iZone )
    {
        auto cgnsZone = std::make_unique< CgnsZone >( *this );
        CgnsZone * zone = cgnsZone.get();
        this->AddCgnsZone( std::move( cgnsZone ) );

        zone->Create();
    }
}

void CgnsBase::ReadCgnsBaseBasicInfo()
{
    CgnsTraits::char33 cgnsBaseName;

    double double_base_id;
    cg_base_id( this->cgnsFile->fileId, this->baseId, & double_base_id );
    std::cout << "   double_base_id = " << double_base_id << "\n";
    //Check the cell and physical dimensions of the bases.
    cg_base_read( this->cgnsFile->fileId, this->baseId, cgnsBaseName, & this->celldim, & this->phydim );
    this->baseName = cgnsBaseName;
    std::cout << "   baseId = " << this->baseId << " baseName = " << cgnsBaseName << "\n";
    std::cout << "   cell dim = " << this->celldim << " physical dim = " << this->phydim << "\n";
}

void CgnsBase::DumpCgnsBaseBasicInfo()
{
    cg_base_write( this->cgnsFile->fileId, this->baseName.c_str(), this->celldim, this->phydim, &this->baseId );
    std::cout << " baseId = " << this->baseId << " baseName = " << this->baseName << "\n";
}

void CgnsBase::ReadNumberOfCgnsZones()
{
    int zoneCount = -1;
    const int status = cg_nzones( this->cgnsFile->fileId, this->baseId, & zoneCount );
    if ( status != CG_OK )
    {
        throw std::runtime_error( "CgnsBase::ReadNumberOfCgnsZones: " + std::string( cg_get_error() ) );
    }
    if ( zoneCount < 0 )
    {
        throw std::runtime_error( "CgnsBase::ReadNumberOfCgnsZones: CGNS returned a negative zone count" );
    }

    this->nZones = zoneCount;
}

CgnsZone * CgnsBase::CreateCgnsZone()
{
    auto cgnsZone = std::make_unique< CgnsZone >( *this );
    CgnsZone * zone = cgnsZone.get();
    this->AddCgnsZone( std::move( cgnsZone ) );
    zone->Create();
    return zone;
}

void CgnsBase::CreateCgnsZones( int nZones )
{
    if ( nZones < 0 )
    {
        throw std::invalid_argument( "CgnsBase::CreateCgnsZones: zone count cannot be negative" );
    }

    if ( ! this->cgnsZones.empty() )
    {
        throw std::logic_error( "CgnsBase::CreateCgnsZones: zones have already been created" );
    }

    this->nZones = nZones;
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        this->CreateCgnsZone();
    }
}

void CgnsBase::ConstructZoneNameMap()
{
    std::map< std::string, int > stagedZoneNameMap;
    for ( size_t iZone = 0; iZone < this->cgnsZones.size(); ++ iZone )
    {
        const CgnsZone * cgnsZone = this->cgnsZones[ iZone ].get();
        const auto inserted = stagedZoneNameMap.emplace( cgnsZone->zoneName, cgnsZone->zId );
        if ( ! inserted.second )
        {
            throw std::runtime_error( "CgnsBase::ConstructZoneNameMap: duplicate zone name '" + cgnsZone->zoneName + "'" );
        }
    }

    this->zoneNameMap.swap( stagedZoneNameMap );
}

void CgnsBase::ReadAllCgnsZones()
{
    std::cout << "** Reading CGNS Grid In Base " << this->baseId << "\n";
    std::cout << "   Reading CGNS Family Specified BC \n";
    this->ReadFamilySpecifiedBc();
    std::cout << "   numberOfCgnsZones       = " << this->nZones << "\n\n";

    for ( size_t iZone = 0; iZone < this->cgnsZones.size(); ++ iZone )
    {
        std::cout << "==>iZone = " << iZone << " numberOfCgnsZones = " << this->cgnsZones.size() << "\n";
        CgnsZone * cgnsZone = this->GetCgnsZone( static_cast< int >( iZone ) );
        cgnsZone->ReadCgnsGrid();
    }
}

void CgnsBase::DumpAllCgnsZones()
{
    std::cout << "** Dumping CGNS Grid In Base " << this->baseId << "\n";
    std::cout << "   Dumping CGNS Family Specified BC \n";
    //this->ReadFamilySpecifiedBc();
    std::cout << "   numberOfCgnsZones       = " << this->nZones << "\n\n";

    for ( size_t iZone = 0; iZone < this->cgnsZones.size(); ++ iZone )
    {
        std::cout << "==>iZone = " << iZone << " numberOfCgnsZones = " << this->cgnsZones.size() << "\n";
        CgnsZone * cgnsZone = this->GetCgnsZone( static_cast< int >( iZone ) );
        cgnsZone->DumpCgnsGrid();
    }
}

void CgnsBase::ProcessCgnsZones()
{
    this->ConvertToInnerDataStandard();

    this->ConstructZoneNameMap();

    for ( size_t iZone = 0; iZone < this->cgnsZones.size(); ++ iZone )
    {
        std::cout << "==>iZone = " << iZone << " numberOfCgnsZones = " << this->cgnsZones.size() << "\n";
        std::cout << "cgnsZone->SetPeriodicBc\n";
        CgnsZone * cgnsZone = this->GetCgnsZone( static_cast< int >( iZone ) );
        cgnsZone->SetPeriodicBc();
    }
}

void CgnsBase::ConvertToInnerDataStandard()
{
    std::cout << "   ConvertToInnerDataStandard \n";

    for ( size_t iZone = 0; iZone < this->cgnsZones.size(); ++ iZone )
    {
        std::cout << "==>iZone = " << iZone << " numberOfCgnsZones = " << this->cgnsZones.size() << "\n";
        CgnsZone * cgnsZone = this->GetCgnsZone( static_cast< int >( iZone ) );
        cgnsZone->ConvertToInnerDataStandard();
    }
}

void CgnsBase::SetFamilyBc( BCType_t & bcType, const std::string & bcRegionName )
{
    this->familyBc->SetFamilyBc( bcType, bcRegionName );
}

BCType_t CgnsBase::GetFamilyBcType( const std::string & bcFamilyName )
{
    return this->familyBc->GetFamilyBcType( bcFamilyName );
}

void CgnsBase::ReadFamilySpecifiedBc()
{
    this->familyBc = std::make_unique< CgnsFamilyBc >( this );
    this->familyBc->ReadFamilySpecifiedBc();
}

CgnsZone * CgnsBase::WriteZoneInfo( const std::string & zoneName, ZoneType_t zoneType, cgsize_t * isize )
{
    auto cgnsZone = std::make_unique< CgnsZone >( *this );
    CgnsZone * zone = cgnsZone.get();
    zone->WriteZoneInfo( zoneName, zoneType, isize );

    this->AddCgnsZone( std::move( cgnsZone ) );
    return zone;
}

CgnsZone * CgnsBase::WriteZone( const std::string & zoneName )
{
    cgsize_t isize[ 9 ];
    this->SetTestISize( isize );

    return this->WriteZoneInfo( zoneName, CGNS_ENUMV( Structured ), isize );
}

void CgnsBase::SetTestISize( cgsize_t * isize )
{
    int nijk = 5;
    for ( int n = 0; n < 3; n ++ )
    {
        isize[ n     ] = nijk;
        isize[ n + 3 ] = nijk - 1;
        isize[ n + 6 ] = 0;
    }
}

void CgnsBase::GoToBase()
{
    cg_goto( this->cgnsFile->fileId, this->baseId, "end" );
}

void CgnsBase::GoToNode( const std::string & nodeName, int ith )
{
    cg_goto( this->cgnsFile->fileId, this->baseId, nodeName.c_str(), ith, NULL );
}

void CgnsBase::ReadArray()
{
    CgnsUserData cgnsUserData( this );
    cgnsUserData.ReadUserData();
}

void CgnsBase::ReadReferenceState()
{
    this->GoToBase();

    CGNS_ENUMT(DataClass_t) id;
    cg_dataclass_read( & id );
    std::cout << "DataClass id = " << id << "\n";
    std::cout << "DataClass = " << DataClassName[ id ] << "\n";

    char * state;
    cg_state_read( & state );
    std::cout << "ReferenceState = " << state << "\n";

    this->GoToNode( "ReferenceState_t", 1 );
    int narrays = -1;
    cg_narrays( & narrays );
    std::cout << " narrays = " << narrays << "\n";

    for ( int n = 1; n <= narrays; ++ n )
    {
        CGNS_ENUMT(DataType_t) idata;
        int idim;
        cgsize_t idimvec;
        char arrayname[33];
        cg_array_info( n, arrayname, & idata, & idim, & idimvec );
        std::cout << " DataTypeName = " << DataTypeName[ idata ] << "\n";
        double data;
        cg_array_read_as( n, CGNS_ENUMV(RealDouble), & data );
        std::cout << "Variable = " << arrayname << "\n";
        std::cout << "   data = " << data << "\n";
    }
}

void CgnsBase::ReadBaseDescriptor()
{
    this->GoToBase();

    //find out how many descriptors are here:
    int ndescriptors = -1;
    cg_ndescriptors( & ndescriptors );
    std::cout << " ndescriptors = " << ndescriptors << "\n";
    for ( int n = 1; n <= ndescriptors; ++ n )
    {
        //read descriptor
        char * text = nullptr;
        char name[ 33 ];
        cg_descriptor_read( n, name, & text );
        std::unique_ptr< char[] > descriptorText( text );
        std::cout << "The descriptor is : " << name << "," << descriptorText.get() << "\n";
    }
}

void CgnsBase::ReadConvergence()
{
    this->GoToBase();

    int nIterations;
    char * text = nullptr;
    cg_convergence_read( &nIterations, & text );
    std::unique_ptr< char[] > convergenceText( text );
    std::cout << "nIterations = " << nIterations << " text = " << convergenceText.get() << "\n";

    this->GoToNode( "ConvergenceHistory_t", 1 );
    int narrays = -1;
    cg_narrays( & narrays );
    std::cout << " narrays = " << narrays << "\n";

    for ( int n = 1; n <= narrays; ++ n )
    {
        CGNS_ENUMT( DataType_t ) itype;

        int idim;
        cgsize_t idimvec;
        char arrayname[ 33 ];
        cg_array_info( n, arrayname, & itype, & idim, & idimvec );
        std::vector< double > varArray( idimvec );
        std::cout << "Datatype = " << itype << " DataTypeName = " << DataTypeName[ itype ] << "\n";
        cg_array_read_as( n, itype, &varArray[ 0 ] );
        std::cout << " VarArrayName = " << arrayname << "\n";
        for ( int i = 0; i < idimvec; ++ i )
        {
            std::cout << varArray[ i ] << " ";
        }
        std::cout << "\n";
    }
}

void CgnsBase::ReadCgnsZones()
{
    if ( ! this->cgnsZones.empty() )
    {
        throw std::logic_error( "CgnsBase::ReadCgnsZones: zones have already been read or allocated" );
    }

    this->ReadNumberOfCgnsZones();
    if ( this->nZones < 0 )
    {
        throw std::runtime_error( "CgnsBase::ReadCgnsZones: CGNS returned a negative zone count" );
    }

    for ( int iZone = 0; iZone < this->nZones; ++ iZone )
    {
        int zoneId = iZone + 1;

        auto cgnsZone = std::make_unique< CgnsZone >( *this );
        CgnsZone * zone = cgnsZone.get();
        zone->zId = zoneId;
        this->AddCgnsZone( std::move( cgnsZone ) );
        zone->ReadCgnsZoneBasicInfo();
    }
}

void CgnsBase::ReadFlowEqn()
{
    this->ReadCgnsZones();

    for ( size_t iZone = 0; iZone < this->cgnsZones.size(); ++ iZone )
    {
        CgnsZone * cgnsZone = this->GetCgnsZone( static_cast< int >( iZone ) );
        cgnsZone->ReadFlowEqn();
    }
}

#endif
EndNameSpace
