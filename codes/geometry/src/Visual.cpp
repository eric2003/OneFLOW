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
#include "Visual.h"
#include "Mesh.h"
#include "FaceMesh.h"
#include "CellMesh.h"
#include "NodeMesh.h"
#include "FaceTopo.h"
using namespace std;

BeginNameSpace( ONEFLOW )

int Visual::numberOfWords = 5;

Visual::Visual()
{
    ;
}

Visual::~Visual()
{
    ;
}

void Visual::SetNumberOfWords( int numberOfWords )
{
    Visual::numberOfWords = numberOfWords;
}

void Visual::Show( Mesh * mesh )
{
    fstream file;
    file.open( "tecplot.dat", ios_base::out );

    Visual::DumpTitle( file, mesh );
    Visual::DumpCoordinate( file, mesh );
    Visual::DumpFaceNodesLink( file, mesh );
    Visual::DumpFaceElementLink( file, mesh );

    file.close();
    file.clear();
}

void Visual::DumpTitle( fstream & file, Mesh * mesh )
{
    StringField titleOfTecplot;
    titleOfTecplot.push_back( "title=\"THE FLOW FIELD OF ONEFLOW\"" );
    titleOfTecplot.push_back( "variables=" );
    titleOfTecplot.push_back( "\"x\"" );
    titleOfTecplot.push_back( "\"y\"" );
    titleOfTecplot.push_back( "\"z\"" );

    for ( UInt i = 0; i < titleOfTecplot.size(); ++ i )
    {
        file << titleOfTecplot[ i ] << endl;
    }

    UInt totalNumFaceNodes = mesh->faceMesh->CalcTotalFaceNodes();
    UInt nFace = mesh->faceMesh->GetNFace();
    UInt numberOfNodes = mesh->nodeMesh->GetNumberOfNodes();
    UInt numberOfCells = mesh->cellMesh->GetNumberOfCells();

    UInt numberOfWordsInEachLine = 5;

    file << "ZONE\n";
    file << "ZoneType = FEPolygon\n";

    file << "Nodes    = " << numberOfNodes << endl;
    file << "Faces    = " << nFace << endl;
    file << "Elements = " << numberOfCells << endl;
    file << "TotalNumFaceNodes = " << totalNumFaceNodes << endl;
    file << "NumConnectedBoundaryFaces = 0\n";
    file << "TotalNumBoundaryConnections = 0\n";
}

void Visual::DumpCoordinate( fstream & file, Mesh * mesh )
{
    Visual::DumpCoordinate( file, mesh->nodeMesh->xN );
    Visual::DumpCoordinate( file, mesh->nodeMesh->yN );
    Visual::DumpCoordinate( file, mesh->nodeMesh->zN );
}

void Visual::DumpCoordinate( fstream & file, RealField & coordinate )
{
    UInt numberOfNodes = coordinate.size();
    for ( UInt iNode = 0; iNode < numberOfNodes; ++ iNode )
    {
        file << coordinate[ iNode ] << " ";
        if ( ( iNode + 1 ) % Visual::numberOfWords == 0 ) file << endl;
    }
    //If it's not full, the end line needs a line break
    if ( numberOfNodes % Visual::numberOfWords != 0 ) file << endl;
}

void Visual::DumpFaceNodesLink( fstream & file, Mesh * mesh )
{
    UInt nodeCount = 0;
    UInt nFace = mesh->faceMesh->GetNFace();
    for ( UInt iFace = 0; iFace < nFace; ++ iFace )
    {
        int numberOfNodesOnFace = mesh->faceMesh->faceTopo->f2n[ iFace ].size();
        for ( int iNodeOfFace = 0; iNodeOfFace < numberOfNodesOnFace; ++ iNodeOfFace )
        {
            file << mesh->faceMesh->faceTopo->f2n[ iFace ][ iNodeOfFace ] + 1 << " ";
            if ( ( nodeCount + 1 ) % Visual::numberOfWords == 0 ) file << endl;
            nodeCount ++;
        }
    }
    if ( nodeCount % Visual::numberOfWords != 0 ) file << endl;
}

void Visual::DumpFaceElementLink( fstream & file, Mesh * mesh )
{
    UInt nFace = mesh->faceMesh->GetNFace();
    UInt numberOfCells = mesh->cellMesh->GetNumberOfCells();

    Visual::DumpFaceElementLink( file, nFace, numberOfCells, mesh->faceMesh->faceTopo->lCell );
    Visual::DumpFaceElementLink( file, nFace, numberOfCells, mesh->faceMesh->faceTopo->rCell );
}

void Visual::DumpFaceElementLink( fstream & file, UInt nFace, UInt numberOfElements, IntField & faceElementIndex )
{
    for ( UInt iFace = 0; iFace < nFace; ++ iFace )
    {
        int elementIndex = faceElementIndex[ iFace ] + 1;
        if ( elementIndex > numberOfElements || elementIndex < 0 ) elementIndex = 0;

        file << elementIndex << " ";
        if ( ( iFace + 1 ) % Visual::numberOfWords == 0 ) file << endl;
    }
    if ( nFace % Visual::numberOfWords != 0 ) file << endl;
}

EndNameSpace