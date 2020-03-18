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

#include "GridMachine.h"
#include "PointMachine.h"
#include "LineMachine.h"
#include "DomainMachine.h"
#include "BlockMachine.h"
#include "Dimension.h"
#include "DataBase.h"
#include "FileIO.h"
#include "HXMath.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

GridMachine grid_Machine;

GridMachine::GridMachine()
{
}

GridMachine::~GridMachine()
{
}

void GridMachine::ReadScript()
{
    string fileName = GetDataValue< string >( "gridLayoutFileName" );
    string separator = " =\r\n\t#$,;\"(){}";

    FileIO ioFile;

    ioFile.OpenPrjFile( fileName, ios_base::in );
    ioFile.SetDefaultSeparator( separator );

    while ( ! ioFile.ReachTheEndOfFile() )
    {
        bool resultFlag = ioFile.ReadNextMeaningfulLine();
        if ( ! resultFlag ) break;

        string keyWord = ioFile.ReadNextWord();
        string word;

        if ( keyWord == "Point" )
        {
            int id = ioFile.ReadNextDigit< int >();

            Real x = ioFile.ReadNextDigit< Real >();
            Real y = ioFile.ReadNextDigit< Real >();
            Real z = ioFile.ReadNextDigit< Real >();

            point_Machine.AddPoint( x, y, z, id );
        }
        else if ( keyWord == "Line" )
        {
            int id = ioFile.ReadNextDigit< int >();
            int p1 = ioFile.ReadNextDigit< int >();
            int p2 = ioFile.ReadNextDigit< int >();

            line_Machine.AddLine( p1, p2, id );
        }
        else if ( keyWord == "Dim" )
        {
            line_Machine.AddDimension( & ioFile );
        }
        else if ( keyWord == "Ds" )
        {
            line_Machine.AddDs( & ioFile );
        }
        else if ( keyWord == "Boundary" )
        {
            domain_Machine.AddBcType( & ioFile );
        }
        else if ( keyWord == "Add" )
        {
            block_Machine.AddFaceToBlock( & ioFile );
        }
        
    };

    ioFile.CloseFile();
}

void GridMachine::GeneGrid()
{
    GenerateAllLineMesh();
    GenerateFaceBlockLink();
}

void GridMachine::GenerateFaceBlockLink()
{
    block_Machine.GenerateFaceBlockLink();
}

void GridMachine::GenerateAllLineMesh()
{
    line_Machine.GenerateAllLineMesh();
}

void GridMachine::DumpStandardGrid()
{
    block_Machine.DumpStandardGrid();
}

void GridMachine::ConstructBlockTopo()
{
    block_Machine.ConstructBlockTopo();
}


EndNameSpace