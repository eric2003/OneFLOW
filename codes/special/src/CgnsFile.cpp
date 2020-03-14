/*-----------------------this->----------------------------------------------------*\
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

#include "CgnsFile.h"
#include "CgnsFactory.h"
#include "Prj.h"
#include "StrUtil.h"
#include "CgnsZone.h"
#include <iostream>
#include <iomanip>
using namespace std;

BeginNameSpace( ONEFLOW )

CgnsFile::CgnsFile()
{
    this->openStatus = CG_ERROR;
}

CgnsFile::CgnsFile( const string & fileName, int openMode )
{
    this->fileName = fileName;
    this->openMode = openMode;
    this->OpenCgnsFile( fileName, openMode );
}

CgnsFile::~CgnsFile()
{
    if ( this->openStatus != CG_OK ) return;
    this->CloseCgnsFile();
}

void CgnsFile::OpenCgnsFile( const string & fileName, int cgnsOpenMode )
{
    this->openStatus = cg_open( fileName.c_str(), cgnsOpenMode, & this->fileId );
    cout << " Current CGNS File Index = " << this->fileId << "\n";
    if ( this->openStatus != CG_OK )
    {
        cg_error_exit();
    }
}

void CgnsFile::CloseCgnsFile()
{
    cg_close( this->fileId );
}

void CgnsFile::WriteBase( const string & baseName )
{
    int celldim = 3;
    int physdim = 3;
    this->WriteBase( baseName, celldim, physdim );
}

void CgnsFile::WriteBase( const string & baseName, int celldim, int physdim )
{
    int baseId = -1;
    cg_base_write( fileId, baseName.c_str(), celldim, physdim, & baseId );
    cout << " CGNS Base index = " << baseId << "\n";
    this->currBaseId = baseId;
}

void CgnsFile::WriteCoor()
{
    #define NUM_SIDE 5
    float coord[ NUM_SIDE * NUM_SIDE * NUM_SIDE ];

    cgsize_t size[ 9 ];

    for ( int n = 0; n < 3; n ++ )
    {
        size[ n     ] = NUM_SIDE;
        size[ n + 3 ] = NUM_SIDE - 1;
        size[ n + 6 ] = 0;
    }

    int nzones = 250;

    for ( int nz = 1; nz <= nzones; nz ++ )
    {
        string name = AddString( "Zone", nz );
        int cgzone = -1;
        int cgcoord = -1;
        cg_zone_write( this->fileId, this->currBaseId, name.c_str(), size, CGNS_ENUMV( Structured ), &cgzone );
        //cg_coord_write( this->fileId, this->currBaseId, cgzone, CGNS_ENUMV( RealSingle ), "CoordinateX", coord, &cgcoord );
        //cg_coord_write( this->fileId, this->currBaseId, cgzone, CGNS_ENUMV( RealSingle ), "CoordinateY", coord, &cgcoord );
        //cg_coord_write( this->fileId, this->currBaseId, cgzone, CGNS_ENUMV( RealSingle ), "CoordinateZ", coord, &cgcoord );
    }
}

void CgnsFile::WriteLink( const string & linkFileName )
{
    int nzones = 250;
    int curr_base_id = 1;
    for ( int nz = 1; nz <= nzones; nz ++ )
    {
        string name     = AddString( "Link to Zone", nz );
        string linkpath = AddString( "/Base/Zone", nz );

        cg_goto( this->fileId, curr_base_id, "end" );
        cg_link_write( name.c_str(), linkFileName.c_str(), linkpath.c_str() );
    }
}


EndNameSpace