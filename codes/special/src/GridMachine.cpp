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
#include "GridLayout.h"
#include "GridLayoutParser.h"
#include "PointMachine.h"
#include "LineMachine.h"
#include "DomainMachine.h"
#include "BlockMachine.h"


BeginNameSpace( ONEFLOW )

GridMachine grid_Machine;

GridMachine::GridMachine()
{
}

GridMachine::~GridMachine()
{
}

void GridMachine::Run( const std::string & fileName )
{
    this->ResetState();
    try
    {
        const GridLayout layout = GridLayoutParser().Parse( fileName );
        this->ApplyLayout( layout );
        this->GenerateGrid();
    }
    catch ( ... )
    {
        this->ResetState();
        throw;
    }
    this->ResetState();
}

void GridMachine::ResetState()
{
    block_Machine.Reset();
    line_Machine.Reset();
    point_Machine.Reset();
    domain_Machine.Reset();
}

void GridMachine::ApplyLayout( const GridLayout & layout )
{
    for ( const auto & point : layout.points )
    {
        point_Machine.AddPoint( point.x, point.y, point.z, point.id );
    }

    for ( const auto & line : layout.lines )
    {
        line_Machine.AddLine( line.p1, line.p2, line.id );
    }

    for ( const auto & circle : layout.circles )
    {
        line_Machine.AddCircle( circle.p1, circle.pc, circle.p2, circle.id );
    }

    for ( const auto & dimension : layout.dimensions )
    {
        line_Machine.SetDimension( dimension.id, dimension.pointCount );
    }

    for ( const auto & distribution : layout.distributions )
    {
        line_Machine.SetDistribution( distribution );
    }

    for ( const auto & boundary : layout.boundaries )
    {
        domain_Machine.SetBcType( boundary.id, boundary.boundaryType );
    }

    for ( const auto & relation : layout.lineToFaces )
    {
        block_Machine.AddLineToFace(
            relation.faceId, relation.position, relation.lineId );
    }

    for ( const auto & relation : layout.faceToBlocks )
    {
        block_Machine.AddFaceToBlock(
            relation.blockId, relation.position, relation.faceId );
    }
}

void GridMachine::GenerateGrid()
{
    // BlockFaceSolver owns the complete mesh-generation sequence.
    block_Machine.GenerateGrid();
}

EndNameSpace
