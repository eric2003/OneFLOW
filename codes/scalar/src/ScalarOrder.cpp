/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2024 He Xin and the OneFLOW contributors.
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

#include "ScalarOrder.h"
#include "ScalarSolver.h"
#include "Numpy.h"
#include "HXMath.h"
#include <iostream>
#include <iomanip>
#include <vector>


BeginNameSpace( ONEFLOW )

ScalarPara::ScalarPara()
{
}

ScalarPara::~ScalarPara()
{
}

ScalarOrder::ScalarOrder()
{
}

ScalarOrder::~ScalarOrder()
{
}

void ScalarOrder::Run1()
{
    int nx0 = 40;
    double len = 2.0;
    double dx = len / ( nx0 );
    int nt = 25;
    double dt = 0.025;
    double c  = 1;
    double timeN = 0.625;

    ScalarPara para;
    para.c = 1.0;
    para.nx = nx0 + 1;
    para.dx = dx;
    para.timeN = timeN;
    para.len = len;
    para.dt = dt;
    para.nt = nt;

    std::vector< std::vector< double > > duList;
    std::vector< std::vector< double > > xList;

    int nTest = 5;

    ScalarSolver * scalarSolver = new ScalarSolver();
    int icoef = 1;
    for ( int i = 0; i < nTest; ++ i )
    {
        icoef *= 2;
        para.nt = nt * icoef;
        para.dt = para.timeN / para.nt;
        para.nx = nx0 * icoef + 1;
        para.dx = len / ( para.nx - 1.0 );
        std::cout << " i = " << i << " dt = " << para.dt << " para.nx=" << para.nx << "  para.dx = " << para.dx << "\n";
        scalarSolver->RunTest( &para );
        duList.push_back( para.du );
        xList.push_back( para.x );
    }
   
    delete scalarSolver;

    Numpy::AnalysisNew( "analysisNew.plt", xList, duList );
   
}

void ScalarOrder::Run()
{
    int nx0 = 40;
    double len = 2.0;
    double dx = len / ( nx0 );
    int nt = 25;
    double dt = 0.025;
    double c  = 1;
    double timeN = 0.625;
    //double timeN = 0.25;

    ScalarPara para;
    para.c = 1.0;
    para.nx = nx0 + 1;
    para.dx = dx;
    para.timeN = timeN;
    para.len = len;
    para.dt = dt;
    para.nt = nt;

    std::vector< std::vector< double > > duList;
    std::vector< std::vector< double > > xList;

    std::vector< double > l1NormList, l2NormList;
    std::vector< double > dxList;

    int nTest = 10;

    ScalarSolver * scalarSolver = new ScalarSolver();
    int icoef = 1;
    for ( int i = 0; i < nTest; ++ i )
    {
        para.nx = nx0 * icoef + 1;
        para.dx = len / ( para.nx - 1.0 );

        //para.nt = nt;
        //para.dt = para.dx / 2.0;
        //para.timeN = para.dt * para.nt;
        para.nt = nt * icoef;
        para.dt = para.timeN / para.nt;
        std::cout << " i = " << i << " dt = " << para.dt << " para.nx=" << para.nx << "  para.dx = " << para.dx << "\n";
        scalarSolver->RunTest( &para );
        duList.push_back( para.du );
        xList.push_back( para.x );
        dxList.push_back( para.dx );
        l1NormList.push_back( para.l1Norm );
        l2NormList.push_back( para.l2Norm );
        icoef *= 2;
    }
    delete scalarSolver;

    Numpy::AnalysisNew( "analysisNew.plt", xList, duList );
    Numpy::DrawL1Norm( "l1Norm.plt", dxList, l1NormList );
    Numpy::DrawNorms( "l1l2Norm.plt", dxList, l1NormList, l2NormList );

    std::vector< double > l1p, l2p;
    for ( int i = 1; i < l1NormList.size(); ++ i )
    {
        double r1 = l1NormList[ i - 1 ] / l1NormList[ i ];
        double r2 = l2NormList[ i - 1 ] / l2NormList[ i ];
        double or1 = log( r1 ) / log( 2 );
        double or2 = log( r2 ) / log( 2 );
        l1p.push_back( or1 );
        l2p.push_back( or2 );
        int width = 8;
        std::cout << " i = " << i;
        std::cout << " order1 = " << std::setw( width ) << std::setiosflags( std::ios::fixed ) << std::setprecision( 6 ) << or1;
        std::cout << " order2 = " << std::setw( width ) << std::setiosflags( std::ios::fixed ) << std::setprecision( 6 ) << or2 << "\n";
    }
    

}

EndNameSpace
