/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2022 He Xin and the OneFLOW contributors.
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
#include "HXCgns.h"
#include <string>
#include <vector>


BeginNameSpace( ONEFLOW )

class CgnsBase;
class CgnsFile;

class CgnsTest
{
public:
    CgnsTest();
    ~CgnsTest();
public:
    std::string fileName;
public:
    void Init();
    void Run();
    void Test();
public:
    void SetDefaultGridName();
    void WriteReferenceState();
    void ReadReferenceState();
    void WriteConvergence();
    void ReadConvergence();
    void WriteDescriptor();
    void ReadDescriptor();
    void WriteSimpleMultiBaseTest();
    void ReadSimpleMultiBaseTest();
    void TestCgnsLink();
    void WriteFlowEqn();
    void ReadFlowEqn();
private:
    void SetISize( cgsize_t * isize );
public:
    void WriteDouble( const std::string & varName, const double & varValue );
public:
    void ReadEmptyCgnsFile();
    void WriteEmptyCgnsFile();
public:
    void WriteArray();
    void WriteArray( CgnsFile * cgnsFile, CgnsBase * cgnsBase );
    void ReadArray();
    void GetArray( std::vector< std::vector< float > > & myfloat2d );
    void WriteTest();
public:
    int read_bcpnts_unst();
    int write_bcpnts_unst();
    int write_grid_unst();
    int read_grid_unst();
    void mytest_read();
    void mytest_write();
};


EndNameSpace
