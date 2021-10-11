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

#include "Visualize.h"

BeginNameSpace( ONEFLOW )

Visualize::Visualize()
{
    ;
}

Visualize::~Visualize()
{
    ;
}

int Plot::nWords = 5;
ostringstream * Plot::oss = 0;

Plot::Plot()
{
    ;
}

Plot::~Plot()
{
    ;
}

void Plot::DumpField( RealField & field )
{
    int nNodes = field.size();
    for ( int iNode = 0; iNode < nNodes; ++ iNode )
    {
        ( * Plot::oss ) << field[ iNode ] << " ";
        if ( ( iNode + 1 ) % Plot::nWords == 0 ) ( * Plot::oss ) << endl;
    }
    if ( nNodes % Plot::nWords != 0 ) ( * Plot::oss ) << endl;
}

void Plot::DumpField( IntField & l2g, RealField & field )
{
    int nNodes = l2g.size();
    for ( int iNode = 0; iNode < nNodes; ++ iNode )
    {
        int id = l2g[ iNode ];
        ( * Plot::oss ) << field[ id ] << " ";
        if ( ( iNode + 1 ) % Plot::nWords == 0 ) ( * Plot::oss ) << endl;
    }
    if ( nNodes % Plot::nWords != 0 ) ( * Plot::oss ) << endl;
}

void Plot::DumpFaceNodeLink( LinkField & f2n )
{
    int nFaces = f2n.size();
    int iCount = 0;
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        int nNodes = f2n[ iFace ].size();
        for ( int iNode = 0; iNode < nNodes; ++ iNode )
        {
            ( * Plot::oss ) << f2n[ iFace ][ iNode ] + 1 << " ";
            if ( ( iCount + 1 ) % Plot::nWords == 0 ) ( * Plot::oss ) << endl;
            iCount ++;
        }
    }
    if ( iCount % Plot::nWords != 0 ) ( * Plot::oss ) << endl;
}

void Plot::DumpFaceElementLink( IntField & elementId, int nElem )
{
    int nFaces = elementId.size();
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        int eId = elementId[ iFace ] + 1;
        if ( eId > nElem || eId < 0 ) eId = 0;

        ( * Plot::oss ) << eId << " ";
        if ( ( iFace + 1 ) % Plot::nWords == 0 ) ( * Plot::oss ) << endl;
    }
    if ( nFaces % Plot::nWords != 0 ) ( * Plot::oss ) << endl;
}

void Plot::DumpFaceNodeNumber( LinkField & f2n )
{
    int nFaces = f2n.size();
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        ( * Plot::oss ) << f2n[ iFace ].size() << " ";
        if ( ( iFace + 1 ) % Plot::nWords == 0 ) ( * Plot::oss ) << endl;
    }
    if ( nFaces % Plot::nWords == 0 ) ( * Plot::oss ) << endl;
}

int GetTotalNumFaceNodes( LinkField & f2n )
{
    int totalNumFaceNodes = 0;
    int nSize = f2n.size();
    for ( int i = 0; i < nSize; ++ i )
    {
        totalNumFaceNodes += f2n[ i ].size();
    }
    return totalNumFaceNodes;
}


EndNameSpace
