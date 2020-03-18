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

#pragma once
#include "Configure.h"
#include "HXDefine.h"
#include <fstream>
#include <iomanip>
using namespace std;

BeginNameSpace( ONEFLOW )

class FileO
{
public:
    FileO();
    ~FileO();
public:
    int nWord;
    int nWidth;
    int nCount;
    std::string fileName;
    ios_base::openmode fileOpenMode;
    fstream * file;
    string sep;
public:
    void OpenPrjFile( const string & fileName, const ios_base::openmode & fileOpenMode );
    void CloseFile();

    void WriteEndLine()
    {
        ( * file ) << "\n";
    }

    template< typename T >
    void Write( const T & value )
    {
        ( * file ) << value << sep;
    }

    template< typename T >
    void WriteFormat( const T & value )
    {
        ( * file ) << setw( nWidth ) << value;
    }

    template< typename T >
    void WriteLine( const T & value )
    {
        ( * file ) << value << "\n";
    }

    void DumpCoorAscii( RealField & coor );
};

EndNameSpace