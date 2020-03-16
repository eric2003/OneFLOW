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
#include "CgnsBase.h"
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
    FreeBaseList();

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

CgnsBase * CgnsFile::WriteBase( const string & baseName )
{
    int celldim = 3;
    int physdim = 3;
    return this->WriteBase( baseName, celldim, physdim );
}

CgnsBase * CgnsFile::WriteBase( const string & baseName, int celldim, int physdim )
{
    int baseId = -1;
    cg_base_write( fileId, baseName.c_str(), celldim, physdim, & baseId );
    cout << " CGNS Base index = " << baseId << "\n";
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

CgnsBase * CgnsFile::AddBase( int fileId, const string & baseName, int celldim, int physdim, int baseId )
{
    CgnsBase * base = new CgnsBase();
    base->fileId = fileId;
    base->baseName = baseName;
    base->celldim = celldim;
    base->phydim = physdim;
    base->baseId = baseId;
    baseList.push_back( base );
    return base;
}


EndNameSpace