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

    OneFLOW is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "BlockMachine.h"
#include "TextFileParser.h"
#include "BlockFaceSolver.h"


BeginNameSpace( ONEFLOW )

BlockMachine block_Machine;

void BlockMachine::AddFaceToBlock( TextFileParser & textFileParser )
{
    std::string word = textFileParser.ReadNextWord();
    if ( word == "L2F" )
    {
        int faceid = textFileParser.ReadNextDigit< int >();
        int pos = textFileParser.ReadNextDigit< int >();
        int lineid = textFileParser.ReadNextDigit< int >();
        this->AddLineToFace( faceid, pos, lineid );
    }
    else if ( word == "F2B" )
    {
        int blockid = textFileParser.ReadNextDigit< int >();
        int pos = textFileParser.ReadNextDigit< int >();
        int faceid = textFileParser.ReadNextDigit< int >();
        this->AddFaceToBlock( blockid, pos, faceid );
    }
}

void BlockMachine::AddLineToFace( int faceId, int position, int lineId )
{
    blkFaceSolver.AddLineToFace( faceId, position, lineId );
}

void BlockMachine::AddFaceToBlock( int blockId, int position, int faceId )
{
    blkFaceSolver.AddFace2Block( blockId, position, faceId );
}

void BlockMachine::GenerateGrid()
{
    blkFaceSolver.GenerateGrid();
}

EndNameSpace
