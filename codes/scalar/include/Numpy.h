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
#include <vector>
#include <string>
using namespace std;

BeginNameSpace( ONEFLOW )

class Numpy
{
public:
    Numpy();
    ~Numpy();
public:
    static void Ones( vector< double > & var );
    static void Set( vector< double > & var, int st, int ed, double v );
    static void Linspace( vector< double > & var, double st, double ed );
    static void Plot( vector< double > & x, vector< double > & f );
    static void Plot( const string & fileName, vector< double > & x, vector< double > & f );
    static void Copy( vector< double > & a, vector< double > & b );
    static void ToTecplot( const string & fileName, vector< double > & x, vector< double > & f );
    static void ToTecplot( const string & fileName, vector< double > & x, vector< double > & u, vector< double > & v );
    static void Analysis( const string & fileName, vector< double > & x, vector< vector< double > > & du );
    static void AnalysisNew( const string & fileName, vector< vector< double > > & x, vector< vector< double > > & du );
    static void DrawL1Norm( const string & fileName, vector< double > & dxList, vector< double > & l1NormList );
    static void DrawNorms( const string & fileName, vector< double > & dxList, vector< double > & l1NormList, vector< double > & l2NormList );

};

EndNameSpace