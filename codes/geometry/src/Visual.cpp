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
#include "Visual.h"
#include "Mesh.h"
#include "FaceMesh.h"
#include "CellMesh.h"
#include "NodeMesh.h"
#include "FaceTopo.h"


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
    std::fstream file;
    file.open( "tecplot.dat", std::ios_base::out );

    Visual::DumpTitle( file, mesh );
    Visual::DumpCoordinate( file, mesh );
    Visual::DumpFaceNodesLink( file, mesh );
    Visual::DumpFaceElementLink( file, mesh );

    file.close();
    file.clear();
}

void Visual::DumpTitle( std::fstream & file, Mesh * mesh )
{
    StringField titleOfTecplot;
    titleOfTecplot.push_back( "title=\"THE FLOW FIELD OF ONEFLOW\"" );
    titleOfTecplot.push_back( "variables=" );
    titleOfTecplot.push_back( "\"x\"" );
    titleOfTecplot.push_back( "\"y\"" );
    titleOfTecplot.push_back( "\"z\"" );

    for ( HXSize_t i = 0; i < titleOfTecplot.size(); ++ i )
    {
        file << titleOfTecplot[ i ] << std::endl;
    }

    HXSize_t totalNumFaceNodes = mesh->faceMesh->CalcTotalFaceNodes();
    HXSize_t nFaces = mesh->faceMesh->GetNFace();
    HXSize_t numberOfNodes = mesh->nodeMesh->GetNumberOfNodes();
    HXSize_t numberOfCells = mesh->cellMesh->GetNumberOfCells();

    HXSize_t numberOfWordsInEachLine = 5;

    file << "ZONE\n";
    file << "ZoneType = FEPolygon\n";

    file << "Nodes    = " << numberOfNodes << std::endl;
    file << "Faces    = " << nFaces << std::endl;
    file << "Elements = " << numberOfCells << std::endl;
    file << "TotalNumFaceNodes = " << totalNumFaceNodes << std::endl;
    file << "NumConnectedBoundaryFaces = 0\n";
    file << "TotalNumBoundaryConnections = 0\n";
}

void Visual::DumpCoordinate( std::fstream & file, Mesh * mesh )
{
    Visual::DumpCoordinate( file, mesh->nodeMesh->xN );
    Visual::DumpCoordinate( file, mesh->nodeMesh->yN );
    Visual::DumpCoordinate( file, mesh->nodeMesh->zN );
}

void Visual::DumpCoordinate( std::fstream & file, RealField & coordinate )
{
    HXSize_t numberOfNodes = coordinate.size();
    for ( HXSize_t iNode = 0; iNode < numberOfNodes; ++ iNode )
    {
        file << coordinate[ iNode ] << " ";
        if ( ( iNode + 1 ) % Visual::numberOfWords == 0 ) file << std::endl;
    }
    //If it's not full, the end line needs a line break
    if ( numberOfNodes % Visual::numberOfWords != 0 ) file << std::endl;
}

void Visual::DumpFaceNodesLink( std::fstream & file, Mesh * mesh )
{
    HXSize_t nodeCount = 0;
    HXSize_t nFaces = mesh->faceMesh->GetNFace();
    for ( HXSize_t iFace = 0; iFace < nFaces; ++ iFace )
    {
        int numberOfNodesOnFace = mesh->faceMesh->faceTopo->faces[ iFace ].size();
        for ( int iNodeOfFace = 0; iNodeOfFace < numberOfNodesOnFace; ++ iNodeOfFace )
        {
            file << mesh->faceMesh->faceTopo->faces[ iFace ][ iNodeOfFace ] + 1 << " ";
            if ( ( nodeCount + 1 ) % Visual::numberOfWords == 0 ) file << std::endl;
            nodeCount ++;
        }
    }
    if ( nodeCount % Visual::numberOfWords != 0 ) file << std::endl;
}

void Visual::DumpFaceElementLink( std::fstream & file, Mesh * mesh )
{
    HXSize_t nFaces = mesh->faceMesh->GetNFace();
    HXSize_t numberOfCells = mesh->cellMesh->GetNumberOfCells();

    Visual::DumpFaceElementLink( file, nFaces, numberOfCells, mesh->faceMesh->faceTopo->lCells );
    Visual::DumpFaceElementLink( file, nFaces, numberOfCells, mesh->faceMesh->faceTopo->rCells );
}

void Visual::DumpFaceElementLink( std::fstream & file, HXSize_t nFaces, HXSize_t numberOfElements, IntField & faceElementIndex )
{
    for ( HXSize_t iFace = 0; iFace < nFaces; ++ iFace )
    {
        int elementIndex = faceElementIndex[ iFace ] + 1;
        if ( elementIndex > numberOfElements || elementIndex < 0 ) elementIndex = 0;

        file << elementIndex << " ";
        if ( ( iFace + 1 ) % Visual::numberOfWords == 0 ) file << std::endl;
    }
    if ( nFaces % Visual::numberOfWords != 0 ) file << std::endl;
}

EndNameSpace
