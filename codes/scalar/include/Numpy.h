/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2023 He Xin and the OneFLOW contributors.
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
#include <fstream>


BeginNameSpace( ONEFLOW )

class Numpy
{
public:
    Numpy();
    ~Numpy();
public:
    static void OpenPrjFile( std::fstream & file, const std::string & fileName, const std::ios_base::openmode & openMode );
    static std::string AddFileName( const std::string & prefix, const std::string & fileName );
    static void Ones( std::vector< double > & var );
    static void Set( std::vector< double > & var, int st, int ed, double v );
    static void Linspace( std::vector< double > & var, double st, double ed );
    static void Plot( std::vector< double > & x, std::vector< double > & f );
    static void Plot( const std::string & fileName, std::vector< double > & x, std::vector< double > & f );
    static void Copy( std::vector< double > & a, std::vector< double > & b );
    static void ToTecplot( const std::string & fileName, std::vector< double > & x, std::vector< double > & f );
    static void ToTecplot( const std::string & fileName, std::vector< double > & x, std::vector< double > & u, std::vector< double > & v );
    static void Analysis( const std::string & fileName, std::vector< double > & x, std::vector< std::vector< double > > & du );
    static void AnalysisNew( const std::string & fileName, std::vector< std::vector< double > > & x, std::vector< std::vector< double > > & du );
    static void DrawL1Norm( const std::string & fileName, std::vector< double > & dxList, std::vector< double > & l1NormList );
    static void DrawNorms( const std::string & fileName, std::vector< double > & dxList, std::vector< double > & l1NormList, std::vector< double > & l2NormList );

};

EndNameSpace
