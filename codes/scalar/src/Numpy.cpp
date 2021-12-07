/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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

#include "Numpy.h"
#include "Prj.h"
#include "StrUtil.h"
#include "FileUtil.h"
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

BeginNameSpace( ONEFLOW )

Numpy::Numpy()
{
}

Numpy::~Numpy()
{
}

string Numpy::AddFileName( const std::string & prefix, const std::string & fileName )
{
    std::string fileNameNew = AddString(prefix, "/", fileName );
    return fileNameNew;
}

void Numpy::OpenPrjFile( std::fstream & file, const std::string & fileName, const std::ios_base::openmode & openMode )
{
    std::string fileNameNew = Numpy::AddFileName( "results", fileName );
    ONEFLOW::OpenPrjFile( file, fileNameNew, std::ios_base::out );
}

void Numpy::Ones( std::vector< double > & var )
{
    int nSize = var.size();
    for ( int i = 0; i < nSize; ++ i )
    {
        var[ i ] = 1.0;
    }
}

void Numpy::Set( std::vector< double > & var, int st, int ed, double v )
{
    for ( int i = st; i < ed; ++ i )
    {
        var[ i ] = v;
    }
}

void Numpy::Linspace( std::vector< double > & var, double st, double ed )
{
    int nSize = var.size();
    double ds = ed - st;
    double dx = ds / ( nSize - 1.0 );
    for ( int i = 0; i < nSize; ++ i )
    {
        var[ i ] = st + i * dx;
    }
}

void Numpy::Copy( std::vector< double > & a, std::vector< double > & b )
{
    int nSize = a.size();
    for ( int i = 0; i < nSize; ++ i )
    {
        b[ i ] = a[ i ];
    }
}

void Numpy::Plot( std::vector< double > & x, std::vector< double > & f )
{
    Numpy::Plot( "plot.plt", x, f );
}

void Numpy::Plot( const std::string & fileName, std::vector< double > & x, std::vector< double > & f )
{
    std::fstream file;
    Numpy::OpenPrjFile( file, fileName, std::ios_base::out );

    int nSize = x.size();
    file << nSize << " ";
    for ( int i = 0; i < nSize; ++ i )
    {
        file << x[ i ] << " " << f[ i ] << " ";
    }

    CloseFile( file );
}

void Numpy::ToTecplot( const std::string & fileName, std::vector< double > & x, std::vector< double > & f )
{
    std::fstream file;
    Numpy::OpenPrjFile( file, fileName, std::ios_base::out );
    int nSize = x.size();

    file << "TITLE = " << "\"OneFLOW X-Y Plot\"" << endl;
    file << "VARIABLES = " << "\"x\", " << "\"u\" " << "\n";
    file << "ZONE T = " << "\"scalar results\"," << " I = " << nSize << "\n";
    for ( int i = 0; i < nSize; ++ i )
    {
        file << x[ i ] << " " << f[ i ] << "\n";
    }

    CloseFile( file );
}

void Numpy::ToTecplot( const std::string & fileName, std::vector< double > & x, std::vector< double > & u, std::vector< double > & v )
{
    std::fstream file;
    Numpy::OpenPrjFile( file, fileName, std::ios_base::out );
    int nSize = x.size();

    file << "TITLE = " << "\"OneFLOW X-Y Plot\"" << endl;
    file << "VARIABLES = " << "\"x\", " << "\"u\", " << "\"v\" " << "\n";
    file << "ZONE T = " << "\"scalar results\"," << " I = " << nSize << "\n";
    for ( int i = 0; i < nSize; ++ i )
    {
        file << x[ i ] << " " << u[ i ] << " " << v[ i ]  << "\n";
    }

    CloseFile( file );
}

void Numpy::Analysis( const std::string & fileName, std::vector< double > & x, std::vector< std::vector< double > > & du )
{
    std::fstream file;
    Numpy::OpenPrjFile( file, fileName, std::ios_base::out );

    file << "TITLE = " << "\"OneFLOW X-Y Plot\"" << endl;

    int nSize = x.size();
    file << "ZONE T = " << "\"scalar analysis results\"," << " I = " << nSize << "\n";
    for ( int i = 0; i < nSize; ++ i )
    {
        file << x[ i ];
        for ( int j = 0; j < du.size(); ++ j )
        {
            file << " " << du[ j ][ i ];
        }
        file << endl;
    }

    CloseFile( file );
}

void Numpy::AnalysisNew( const std::string & fileName, std::vector< std::vector< double > > & x, std::vector< std::vector< double > > & du )
{
    std::fstream file;
    Numpy::OpenPrjFile( file, fileName, std::ios_base::out );

    file << "TITLE = " << "\"OneFLOW X-Y Plot\"" << endl;

    for ( int j = 0; j < du.size(); ++ j )
    {
        int nSize = x[ j ].size();
        file << "ZONE T = " << "\"scalar analysis results\"," << " I = " << nSize << "\n";
        for ( int i = 0; i < nSize; ++ i )
        {
            file << x[ j ][ i ] << " " << du[ j ][ i ] << endl;
        }
    }

    CloseFile( file );
}

void Numpy::DrawL1Norm( const std::string & fileName, std::vector< double > & dxList, std::vector< double > & l1NormList )
{
    std::fstream file;
    Numpy::OpenPrjFile( file, fileName, std::ios_base::out );

    file << "TITLE = " << "\"OneFLOW X-Y Plot\"" << endl;

    int nSize = dxList.size();
    file << "VARIABLES = " << "\"dx\", " << "\"norm\" " << "\n";
    file << "ZONE T = " << "\"scalar Norm analysis results\"," << " I = " << nSize << "\n";
    for ( int i = 0; i < nSize; ++ i )
    {
        file << dxList[ i ] << " " << l1NormList[ i ] << endl;
    }

    CloseFile( file );
}

void Numpy::DrawNorms( const std::string & fileName, std::vector< double > & dxList, std::vector< double > & l1NormList, std::vector< double > & l2NormList )
{
    std::fstream file;
    Numpy::OpenPrjFile( file, fileName, std::ios_base::out );

    file << "TITLE = " << "\"OneFLOW X-Y Plot\"" << endl;

    int nSize = dxList.size();
    file << "VARIABLES = " << "\"dx\" " << "\"norm1\" " << "\"norm2\" " << "\n";
    file << "ZONE T = " << "\"scalar Norm analysis results\"," << " I = " << nSize << "\n";
    for ( int i = 0; i < nSize; ++ i )
    {
        file << dxList[ i ] << " " << l1NormList[ i ] << " " << l2NormList[ i ] << endl;
    }
    CloseFile( file );
}



EndNameSpace
