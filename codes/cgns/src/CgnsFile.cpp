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

#include "CgnsFile.h"
#include "CgnsBase.h"
#include "CgnsFactory.h"
#include "Prj.h"
#include "StrUtil.h"
#include "CgnsZone.h"
#include "Stop.h"
#include <iostream>
#include <iomanip>

BeginNameSpace( ONEFLOW )

std::string GetCgnsFileTypeName( int file_type )
{
    std::string fileTypeName;
    if ( file_type == CG_FILE_ADF )
    {
        fileTypeName = "CG_FILE_ADF";
    }
    else if ( file_type == CG_FILE_HDF5 )
    {
        fileTypeName = "CG_FILE_HDF5";
    }
    else if ( file_type == CG_FILE_ADF2 )
    {
        fileTypeName = "CG_FILE_ADF2";
    }
    else
    {
        fileTypeName = "CG_FILE_NONE";
    }
    return fileTypeName;
}

CgnsFile::CgnsFile()
{
    this->openStatus = CG_ERROR;
}

CgnsFile::CgnsFile( const std::string & fileName, int openMode )
{
    this->OpenCgnsFile( fileName, openMode );
}

CgnsFile::~CgnsFile()
{
    FreeBaseList();

    if ( this->openStatus != CG_OK ) return;
    this->CloseCgnsFile();
}

void CgnsFile::OpenCgnsFile( const std::string & fileName, int cgnsOpenMode )
{
    this->fileName = fileName;
    this->openMode = openMode;

    this->openStatus = cg_open( fileName.c_str(), cgnsOpenMode, & this->fileId );
    std::string stars("**************************************************************");
    std::cout << stars << "\n";
    std::cout << "   CGNS File Index = " << this->fileId << "\n";

    if ( this->openStatus != CG_OK )
    {
        //cg_error_exit();
        Stop( cg_get_error() );
    }

    float fileVersion = -1;
    cg_version( this->fileId, & fileVersion );

    std::cout << "   CGNS File Version = " << std::setiosflags( std::ios::fixed ) << std::setprecision( 4 ) << fileVersion << "\n";

    int precision = -1;
    cg_precision( this->fileId, & precision );

    std::cout << "   CGNS Precision = " << precision << "\n";

    int file_type = -1;
    cg_get_file_type( this->fileId, & file_type );

    std::cout << "   CGNS File Type = " << file_type << " FileTypeName = " << GetCgnsFileTypeName( file_type ) << "\n";
    std::cout << stars << "\n";
}

void CgnsFile::CloseCgnsFile()
{
    cg_close( this->fileId );
}

CgnsBase * CgnsFile::WriteBase( const std::string & baseName )
{
    int celldim = 3;
    int physdim = 3;
    return this->WriteBase( baseName, celldim, physdim );
}

CgnsBase * CgnsFile::WriteBase( const std::string & baseName, int celldim, int physdim )
{
    int baseId = -1;
    cg_base_write( fileId, baseName.c_str(), celldim, physdim, & baseId );
    std::cout << " CGNS Base index = " << baseId << "\n";
    this->currBaseId = baseId;
    return this->AddBase( fileId, baseName, celldim, physdim, baseId );
}

void CgnsFile::FreeBaseList()
{
    for ( int i = 0; i < baseList.size(); ++ i )
    {
        delete baseList[ i ];
    }
}

CgnsBase * CgnsFile::AddBase( int fileId, const std::string & baseName, int celldim, int physdim, int baseId )
{
    CgnsBase * base = new CgnsBase( this );
    base->baseName = baseName;
    base->celldim = celldim;
    base->phydim = physdim;
    base->baseId = baseId;
    baseList.push_back( base );
    return base;
}

void CgnsFile::GoPath( const std::string & path )
{
    cg_gopath( this->fileId, path.c_str() );
}

void CgnsFile::ReadNumberOfBases()
{
    cg_nbases( this->fileId, & this->nBases );
    std::cout << " Total number of CGNS Base = " << this->nBases << "\n";
}

CgnsBase * CgnsFile::CreateCgnsBase()
{
    int iBase = baseList.size();
    int baseId = iBase + 1;

    CgnsBase * cgnsBase = new CgnsBase( this );
    this->baseList.push_back( cgnsBase );
    cgnsBase->baseId = baseId;

    return cgnsBase;
}

void CgnsFile::ReadBases()
{
    this->ReadNumberOfBases();
    for ( int iBase = 0; iBase < this->nBases; ++ iBase )
    {
        int baseId = iBase + 1;
        CgnsBase * cgnsBase = new CgnsBase( this );
        this->baseList.push_back( cgnsBase );
        cgnsBase->baseId = baseId;
        cgnsBase->ReadCgnsBaseBasicInfo();
    }
}

void CgnsFile::ReadArray()
{
    this->ReadBases();
    for ( int iBase = 0; iBase < this->nBases; ++ iBase )
    {
        CgnsBase * cgnsBase = this->baseList[ iBase ];
        cgnsBase->ReadArray();
    }
}

void CgnsFile::ReadReferenceState()
{
    this->ReadBases();
    for ( int iBase = 0; iBase < this->nBases; ++ iBase )
    {
        CgnsBase * cgnsBase = this->baseList[ iBase ];
        cgnsBase->ReadReferenceState();
    }
}

void CgnsFile::ReadBaseDescriptor()
{
    this->ReadBases();
    for ( int iBase = 0; iBase < this->nBases; ++ iBase )
    {
        CgnsBase * cgnsBase = this->baseList[ iBase ];
        cgnsBase->ReadBaseDescriptor();
    }
}

void CgnsFile::WriteBaseDescriptor()
{
    std::cout << "Program write_descriptor\n";
    CgnsBase * cgnsBase = this->WriteBase( "base1" );

    //cg_goto must be called or an error will occur
    cg_goto( this->fileId, cgnsBase->baseId, "end" );
    cg_descriptor_write("Information","info1");
    cg_descriptor_write("hello world","haha ! hello world!");

    cgnsBase = this->WriteBase( "base2" );
    cg_goto( this->fileId, cgnsBase->baseId, "end" );
    cg_descriptor_write("descript1","des1");
    cg_descriptor_write("descript2","des1");
    cg_descriptor_write("descript3","des1");

    cgnsBase = this->WriteBase( "base3" );
    cg_goto( this->fileId, cgnsBase->baseId, "end" );
    cg_descriptor_write("mydes","mydes1");
    cg_descriptor_write("mydes","mydes2");
    cg_descriptor_write("mydes","mydes3");
}

void CgnsFile::ReadConvergence()
{
    this->ReadBases();
    for ( int iBase = 0; iBase < this->nBases; ++ iBase )
    {
        CgnsBase * cgnsBase = this->baseList[ iBase ];
        cgnsBase->ReadConvergence();
    }
}

void CgnsFile::ReadFlowEqn()
{
    this->ReadBases();
    for ( int iBase = 0; iBase < this->nBases; ++ iBase )
    {
        CgnsBase * cgnsBase = this->baseList[ iBase ];
        cgnsBase->ReadFlowEqn();
    }
}

EndNameSpace
