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
#include <string>
#include <fstream>
using namespace std;

BeginNameSpace( ONEFLOW )

void OpenLogFile( int logFileIndex, fstream & file );
void CloseLogFile( fstream & file );
class LogFile;
extern LogFile logFile;

class LogFile
{
public:
    LogFile();
    ~LogFile();
    fstream my_fstream;
public:
    void Open();
    void Close();
};

template< typename T >
LogFile & operator << ( LogFile & f, const T & value )
{
    f.Open();
    f.my_fstream << value;
    f.Close();
    return f;
}

EndNameSpace
