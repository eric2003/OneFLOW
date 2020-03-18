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

#include "FileO.h"
#include "FileUtil.h"
#include "Prj.h"

BeginNameSpace( ONEFLOW )

FileO::FileO()
{
    file = new fstream();
    sep = " ";
    nWord = 5;
    nWidth = 5;
    nCount = 0;
}

FileO::~FileO()
{
    delete file;
}

void FileO::OpenPrjFile( const string & fileName, const ios_base::openmode & fileOpenMode )
{
    this->fileName     = fileName;
    this->fileOpenMode = fileOpenMode;
    ONEFLOW::OpenPrjFile( * file, fileName, fileOpenMode );
}

void FileO::CloseFile()
{
    ONEFLOW::CloseFile( * file );
}

void FileO::DumpCoorAscii( RealField & coor )
{
    int nCountMax = 10000;
    int nPoint = coor.size();
    nWidth = 15;
    for ( int i = 0; i < nPoint; ++ i )
    {
        ( * file ) << setw( nWidth ) << coor[ i ];
        nCount ++;
        if ( nCount % nWord == 0 )
        {
            if ( nCount >= nCountMax ) nCount = 0;
            ( * file ) << "\n";
        }
    }
}


EndNameSpace