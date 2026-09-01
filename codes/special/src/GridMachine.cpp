/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
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
#include "GridPara.h"
#include "TextFileParser.h"
#include "HXMath.h"
#include <iostream>


BeginNameSpace( ONEFLOW )

GridMachine grid_Machine;

GridMachine::GridMachine()
{
}

GridMachine::~GridMachine()
{
}

void GridMachine::Run()
{
    std::string fileName = GetDataValue< std::string >( "gridLayoutFileName" );
    this->ReadScript( fileName );
    this->GeneGrid();
}

void GridMachine::ReadScript( const std::string & fileName )
{
    std::string separator = " =\r\n\t#$,;\"(){}";

    TextFileParser textFileParser;

    textFileParser.OpenPrjFile( fileName, std::ios_base::in );
    textFileParser.SetDefaultSeparator( separator );

    while ( ! textFileParser.ReachTheEndOfFile() )
    {
        bool resultFlag = textFileParser.ReadNextMeaningfulLine();
        if ( ! resultFlag ) break;

        std::string keyWord = textFileParser.ReadNextWord();
        std::string word;

        if ( keyWord == "Point" )
        {
            int id = textFileParser.ReadNextDigit< int >();

            Real x = textFileParser.ReadNextDigit< Real >();
            Real y = textFileParser.ReadNextDigit< Real >();
            Real z = textFileParser.ReadNextDigit< Real >();

            point_Machine.AddPoint( x, y, z, id );
        }
        else if ( keyWord == "Line" )
        {
            int id = textFileParser.ReadNextDigit< int >();
            int p1 = textFileParser.ReadNextDigit< int >();
            int p2 = textFileParser.ReadNextDigit< int >();

            line_Machine.AddLine( p1, p2, id );
        }
        else if ( keyWord == "Circle" )
        {
            int id = textFileParser.ReadNextDigit< int >();
            int p1 = textFileParser.ReadNextDigit< int >();
            int pc = textFileParser.ReadNextDigit< int >();
            int p2 = textFileParser.ReadNextDigit< int >();

            line_Machine.AddCircle( p1, pc, p2, id );
        }
        else if ( keyWord == "Dim" )
        {
            line_Machine.AddDimension( & textFileParser );
        }
        else if ( keyWord == "Ds" )
        {
            line_Machine.AddDs( & textFileParser );
        }
        else if ( keyWord == "Boundary" )
        {
            domain_Machine.AddBcType( & textFileParser );
        }
        else if ( keyWord == "Add" )
        {
            block_Machine.AddFaceToBlock( & textFileParser );
        }
        
    };

    textFileParser.CloseFile();
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


EndNameSpace
