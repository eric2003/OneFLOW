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
#include <list>
#include <string>
#include <map>
using namespace std;

BeginNameSpace( ONEFLOW )

class StrGrid;
class StrBcSetting
{
public:
    StrBcSetting();
    ~StrBcSetting();
public:
    IntField iminList, imaxList, jminList, jmaxList, kminList, kmaxList, bcTypeList;
public:
    map< string, int > boundaryMap;
public:
    void ConstructBcMap();
    int GetBcType( const string & bcTypeName );
public:
    void PushBc( int imin, int imax, int jmin, int jmax, int kmin, int kmax, int bcType );
    void SetBcRegion( StrGrid * grid );
    int CalcNBcRegion();
};


EndNameSpace