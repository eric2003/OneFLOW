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
#include "Test.h"
#include "FileIO.h"
#include "FileUtil.h"
#include "Prj.h"
#include <iostream>
#include <fstream>
using namespace std;

BeginNameSpace( ONEFLOW )

Test::Test()
{
    ;
}

Test::~Test()
{
    ;
}

void Test::Run()
{
    int ni = 101;
    Real xlen = 2.0;
    Real dx = xlen / ( ni - 1 );
    fstream file;
    OpenPrjFile( file, "/test/vencat.dat", ios_base::out );
    StringField title;
    title.push_back( "title=\"THE FLOW FIELD OF ONEFLOW\"" );
    title.push_back( "variables=" );
    title.push_back( "\"x\"" );
    title.push_back( "\"vencat\"" );
    title.push_back( "\"gvencat\"" );

    for ( UInt i = 0; i < title.size(); ++ i )
    {
        file << title[ i ] << endl;
    }
    int nj = 10;
    for ( int j = 0; j < nj; ++ j )
    {
        Real c = 2.0 + j;
        file << " zone  i = " << ni << endl;
        for ( int i = 0; i < ni; ++ i )
        {
            Real x = i * dx;
            Real f = VencatC( x, c );
            Real g = f * x;
            //cout << x << " " << f << endl;
            file << x << " " << f << " " << g << endl;
        }
    }

    CloseFile( file );
}

Real Test::Vencat( Real x )
{
    //Real c = 2.0;
    Real c = 9.0;
    Real v1 = 1 + 2 * x;
    Real v2 = 1 + x + 2 * x * x;
    return v1 / v2;
}

Real Test::VencatC( Real x, Real c )
{
    Real v1 = 1 + c * x;
    Real v2 = 1 + x + c * x * x;
    return v1 / v2;
}


Real Test::VencatEric( Real x )
{
    Real v1 = 1 + 3 * x;
    Real v2 = 1 + x + 2 * x * x;
    return v1 / v2;
}


void FunTest()
{
    Test test;
    test.Run();
}


EndNameSpace
