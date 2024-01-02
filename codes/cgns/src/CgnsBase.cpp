/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2024 He Xin and the OneFLOW contributors.
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
#include "CgnsZoneUtil.h"
#include "StrUtil.h"
#include "Dimension.h"
#include "CgnsFamilyBc.h"
#include "CgnsVariable.h"

#include <iostream>

BeginNameSpace( ONEFLOW )

#ifdef ENABLE_CGNS

CgnsBase::CgnsBase()
{
    this->cgnsFile = 0;
    this->familyBc = 0;
    this->freeFlag = false;
}

CgnsBase::CgnsBase( CgnsFile * cgnsFile )
{
    this->cgnsFile = cgnsFile;
    this->familyBc = 0;
    this->freeFlag = false;
}

CgnsBase::~CgnsBase()
{
    delete this->familyBc;
    if ( this->freeFlag )
    {
        this->FreeZoneList();
    }
}

void CgnsBase::FreeZoneList()
{
    for ( int i = 0; i < cgnsZones.size(); ++ i )
    {
        delete cgnsZones[ i ];
    }
}


CgnsZone * CgnsBase::GetCgnsZone( int iZone )
{
    //iZone base on 0
    return this->cgnsZones[ iZone ];
}

CgnsZone * CgnsBase::GetCgnsZoneByName( const std::string & zoneName )
{
    std::map< std::string, int >::iterator iter;
    iter = zoneNameMap.find( zoneName );
    int iZone = iter->second - 1;
    return this->GetCgnsZone( iZone );
}

int CgnsBase::GetNZones()
{
    return this->cgnsZones.size();
}

void CgnsBase::SetDefaultCgnsBaseBasicInfo()
{
    //this->celldim = Dim::dimension;
    //this->phydim  = Dim::dimension;

    this->celldim = THREE_D;
    this->phydim  = THREE_D;
  
    this->baseName = ONEFLOW::AddString( "Base", this->baseId );
}

void CgnsBase::AddCgnsZone( CgnsZone * cgnsZone )
{
    cgnsZones.push_back( cgnsZone );
    int zId = cgnsZones.size();
    cgnsZone->zId = zId;
}

void CgnsBase::AllocateAllCgnsZones()
{
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        CgnsZone * cgnsZone = new CgnsZone( this );

        this->AddCgnsZone( cgnsZone );

        cgnsZone->Create();
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
    //Read the number of zones in the grid.
    cg_nzones( this->cgnsFile->fileId, this->baseId, & this->nZones );
}

CgnsZone * CgnsBase::CreateCgnsZone()
{
    CgnsZone * cgnsZone = new CgnsZone( this );

    this->AddCgnsZone( cgnsZone );

    cgnsZone->Create();

    return cgnsZone;
}

void CgnsBase::CreateCgnsZones( int nZones )
{
    this->nZones = nZones;
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        CgnsZone * cgnsZone = this->CreateCgnsZone();
    }
}

void CgnsBase::ConstructZoneNameMap()
{
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        CgnsZone * cgnsZone = this->GetCgnsZone( iZone );
        zoneNameMap[ cgnsZone->zoneName ] = cgnsZone->zId;
    }
}

void CgnsBase::ReadAllCgnsZones()
{
    std::cout << "** Reading CGNS Grid In Base " << this->baseId << "\n";
    std::cout << "   Reading CGNS Family Specified BC \n";
    this->ReadFamilySpecifiedBc();
    std::cout << "   numberOfCgnsZones       = " << this->nZones << "\n\n";

    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        std::cout << "==>iZone = " << iZone << " numberOfCgnsZones = " << this->nZones << "\n";
        CgnsZone * cgnsZone = this->GetCgnsZone( iZone );
        cgnsZone->ReadCgnsGrid();
    }
}

void CgnsBase::DumpAllCgnsZones()
{
    std::cout << "** Dumping CGNS Grid In Base " << this->baseId << "\n";
    std::cout << "   Dumping CGNS Family Specified BC \n";
    //this->ReadFamilySpecifiedBc();
    std::cout << "   numberOfCgnsZones       = " << this->nZones << "\n\n";

    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        std::cout << "==>iZone = " << iZone << " numberOfCgnsZones = " << this->nZones << "\n";
        CgnsZone * cgnsZone = this->GetCgnsZone( iZone );
        cgnsZone->DumpCgnsGrid();
    }
}

void CgnsBase::ProcessCgnsZones()
{
    this->ConvertToInnerDataStandard();

    this->ConstructZoneNameMap();

    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        std::cout << "==>iZone = " << iZone << " numberOfCgnsZones = " << this->nZones << "\n";
        std::cout << "cgnsZone->SetPeriodicBc\n";
        CgnsZone * cgnsZone = this->GetCgnsZone( iZone );
        cgnsZone->SetPeriodicBc();
    }
}

void CgnsBase::ConvertToInnerDataStandard()
{
    std::cout << "   ConvertToInnerDataStandard \n";

    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        std::cout << "==>iZone = " << iZone << " numberOfCgnsZones = " << this->nZones << "\n";
        CgnsZone * cgnsZone = this->GetCgnsZone( iZone );
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
    this->familyBc = new CgnsFamilyBc( this );
    this->familyBc->ReadFamilySpecifiedBc();
}

CgnsZone * CgnsBase::WriteZoneInfo( const std::string & zoneName, ZoneType_t zoneType, cgsize_t * isize )
{
    int cgzone = -1;
    cg_zone_write( this->cgnsFile->fileId, this->baseId, zoneName.c_str(), isize, zoneType, & cgzone );
    this->freeFlag = true;

    CgnsZone * cgnsZone = new CgnsZone( this );

    this->AddCgnsZone( cgnsZone );

    cgnsZone->WriteZoneInfo( zoneName, zoneType, isize );

    return cgnsZone;
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
        char *text, name[33];
        cg_descriptor_read( n, name, &text );
        std::cout << "The descriptor is : " << name << "," << text << "\n";
        delete[ ] text;
    }
}

void CgnsBase::ReadConvergence()
{
    this->GoToBase();

    int nIterations;
    char *text;
    cg_convergence_read( &nIterations, &text );
    std::cout << "nIterations = " << nIterations << " text = " << text << "\n";
    delete[ ] text;

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
    this->ReadNumberOfCgnsZones();

    for ( int iZone = 0; iZone < this->nZones; ++ iZone )
    {
        int zoneId = iZone + 1;

        CgnsZone * cgnsZone = new CgnsZone( this );
        cgnsZone->zId = zoneId;
        this->AddCgnsZone( cgnsZone );
        cgnsZone->ReadCgnsZoneBasicInfo();
    }
}

void CgnsBase::ReadFlowEqn()
{
    this->ReadCgnsZones();

    for ( int iZone = 0; iZone < this->nZones; ++ iZone )
    {
        CgnsZone * cgnsZone = this->GetCgnsZone( iZone );
        cgnsZone->ReadFlowEqn();
    }
}

#endif
EndNameSpace
