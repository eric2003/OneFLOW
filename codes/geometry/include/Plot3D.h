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
#include "HXDefine.h"
#include <vector>
#include <string>
#include <fstream>
using namespace std;

BeginNameSpace( ONEFLOW )

const int ASCII  = 0;
const int BINARY = 1;
class GridMediator;
class FileIO;
class ZgridMediator;

class Plot3D
{
public:
    Plot3D();
    ~Plot3D();
public:
    static void ReadPlot3D( GridMediator * gridMediator );
    static void ReadCoor( GridMediator * gridMediator );
    static void ReadCoorBinary( GridMediator * gridMediator );
    static void ReadCoorAscii ( GridMediator * gridMediator );
    static void ReadCoor( FileIO * ioFile, RealField & coordinate );
    static void ReadCoor( FileIO * ioFile, RealField & coor, int total_size );
    static void ReadBc( GridMediator * gridMediator );
public:
    static void DumpCoor( GridMediator * gridMediator );
    static void DumpCoorBinary( GridMediator * gridMediator );
    static void DumpCoorAscii( GridMediator * gridMediator );
    static void DumpCoorAscii( fstream & file, RealField & coor );
    static void DumpBc( GridMediator * gridMediator );
    static void Plot3DToCgns( ZgridMediator * zgridMediator );

};

bool GetPlot3D_NKFlag();

EndNameSpace